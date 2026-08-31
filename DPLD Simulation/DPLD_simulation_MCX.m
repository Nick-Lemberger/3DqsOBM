% This is a simple Matlab script to simulate the detected-photon launch
% distribution (DPLD) required for phase or index retireval in qsOBM with
% the S11142-10 photodetector.
% Author: Nick Sidney Lemberger, 19.05.2026
%% Add the path to MXCLab
addpath(genpath('G:\User\Nick Lemberger\MCXStudio'));
%savepath
%%
mcxlab('gpuinfo') % Check if MXC finds the GPU

%% General description of the simulation. Scaling: Voxel = 0.04mm
% MCX Simulation for the S11142-10 detector taped to sample in
% epi-direction with a 170µm cover slip and small 518F oil filled gap.
clear cfg
cfg.nphoton=1e9; % Number of photons for simulation.
cfg.tstart=0;
cfg.tend=5e-9;
cfg.tstep=5e-9;
cfg.issrcfrom0=0; % Array start at 1 as customary in Matlab

% Define the material properties
n_518F = 1.5035; % Immersion oil 518F
n_glass = 1.507; % BK7 glass
n_si = 3.56; % Silicon of the photodiode

% Define the properties of your scattering medium

% 1% Intralipid in Agar scattering phantom
%n_scatter = 1.323; % Refractive index
%us_scatter = 2.35; % Full scattering coefficient  
%g_scatter = 0.6525; % Scattering anisotropy g

% Brain
%n_scatter =  1.36;
%us_scatter = 50;
%g_scatter = 0.9;

% Chicken muscle tissue
n_scatter =  1.37;
us_scatter = 10;
g_scatter = 0.93;

% Material list definition
% A weak dummy absoprtion of 0.001 is used for calculation of the jacobian
cfg.prop  =  [0, 0, 1, 1;           % Zero voxels, air   
              0.001, 0, 1, n_518F;  % Immersion oil 518F
              0.001, 0, 1, n_si;    % Detector sillicone reflection dummy 
              0.001, 0, 1, n_glass; % 170µm glass cover slip
              0.001, us_scatter, g_scatter, n_scatter]; % Scattering medium
% Set a virtual focal point in the sample.
focalPoint_z = 11; % Sample starts at 10, 40µm per voxel unit

% define the domain
dimxy = 375; % 15mm in x,y plane to cover the full detector
dimz = 250; % 10mm high
cfg.unitinmm = 0.04; % 40µm per voxel
cfg.vol=ones(dimxy,dimxy,dimz); % 15x15x10mm simulation volume

% For when the sample is the scattering medium
cfg.vol(:,:,1:2) = 2; % Detector silicon
cfg.vol(:,:,3:5) = 1; % Detector / glas gap filled with 518F oil
cfg.vol(:,:,6:9) = 3; % Cover slip (160µm as 4*40µm=160µm)
cfg.vol(:,:,10:250) = 4; % Sample 

% For beads in oil closed with cover slip and scattering medium on top 
%cfg.vol(:,:,1:2) = 2; % Detector silicon
%cfg.vol(:,:,3:5) = 1; % Detector / glas gap filled with 518F oil
%cfg.vol(:,:,6:9) = 3; % Cover slip (160µm as 4*40µm=160µm)
%cfg.vol(:,:,10:16) = 1; % Sample, immersed in 518F oil
%cfg.vol(:,:,17:20) = 3; % Cover slip on top (160µm as 4*40µm=160µm)
%cfg.vol(:,:,21:250) = 4; % Scattering medium

% Cut the 2mm hole in the silicone layer 
r_hole = 1 / cfg.unitinmm; % 1mm detector hole radius
[xi, yi, zi] = meshgrid(1:dimxy, 1:dimxy, 1:2);
dist = sqrt((xi - dimxy/2).^2 + (yi - dimxy/2).^2);
cfg.vol(dist < r_hole) = 1;

% Isotropic source at focal point to allow for all refraction angles to contribute 
cfg.srctype='isotropic';
cfg.srcpos=[dimxy/2,dimxy/2,focalPoint_z];
cfg.srcdir=[0,0,1];

% save d: det_id, p: ppath, x: exit position, v: exit direction
cfg.savedetflag='dpxv'; 
cfg.autopilot = 1;

% First, simulate photons reaching the detector from the focal point.
% Detect photons in epi direction (-z) leaving the simulation. 
% Select according to real detector dimentions in post selection.
% Set boundary conditions per-face (-x,-y,-z,+x,+y,+z) to absorption (a)
% and set face -z as detector (001000).
cfg.bc = 'aaaaaa001000';
cfg.isreflect = 1; % Calculate Frenel reflections
cfg.isspecular = 0;
cfg.maxdetphoton = 1e8; % Maximum number of photons to detect

% Preview of simulation volume
%mcxpreview(cfg);

% Start the simulation
[flux,detp,vol,seed] = mcxlab(cfg);
 
% Get the positions of the detected photons
px = (detp.p(:,1)-dimxy/2)*cfg.unitinmm;
py = (detp.p(:,2)-dimxy/2)*cfg.unitinmm;
pz = detp.p(:,3)*cfg.unitinmm;

% Select photons according to dimensions of real detector quadrants
pd_slit = 0.15; % Tiling slit of S11142-10 detector 
pd_hole = 2.3;  % Diameter of the hole in detector S11142-10
pd_aa = 14; % Length of the area in x and y
% Select photons indecies according to the geometry of the active area
s1 = px < pd_aa/2 & px > -pd_aa/2 & py < pd_aa/2 & py > -pd_aa/2; % Select photons in square active area
s2 = sqrt(px.^2 + py.^2) > pd_hole/2; % Remove photons from the center hole
s3 = (px > pd_slit/2 | px < -pd_slit/2) & (py > pd_slit/2 | py < -pd_slit/2); % Remove photons from tiling slits
sQ12 = py > 0; % Select the upper two quadrants
% Create an index selection mask via bitwise AND operation
mask = s1 & s2 & s3 & sQ12;  
clear('s1','s2','s3','sQ12','pd_slit','pd_hole','pd_aa');
% Select photos via index mask
px_sel = px(mask);
py_sel = py(mask);
pz_sel = pz(mask);


% Replay the simulation with photons reaching the actual detector geometry
% Modify to detect the replayed photons leaving the focal volume to get
% their directions
cfg_replay = cfg;
cfg_replay.seed = seed.data(:,mask);  % define cfg.seed using the returned seeds
cfg_replay.detphotons= detp.data(:,mask);  % define cfg.detphotons. using the returned detp data
%cfg_replay.outputtype='jacobian'; % tell mcx to output absorption (mu_a) Jacobian

% Sperical zero-shell around illumination point. Set everything 0 except
% for a small sphere
r_0shell = 4; % Radius of zero shell for detection to catch immediate direction of outgoing photons
[xi, yi, zi] = meshgrid(1:dimxy, 1:dimxy, 1:dimz);
dist = sqrt((xi - dimxy/2).^2 + (yi - dimxy/2).^2 + (zi - focalPoint_z).^2);
cfg_replay.vol(dist > r_0shell) = 0; % all zero outside of small sphere

% Introduce a new detector to detect the photon directions emitted by the
% source that fell within detector geometry in the previous simulation run.
cfg_replay.detpos=[dimxy/2,dimxy/2,focalPoint_z, r_0shell*2]; % New detector
cfg_replay.bc = 'aaaaaa000000'; % Disable old detector
cfg_replay.isreflect = 0; % Disable reflections
% Also turn of scattering and refraction by setting all µs=0, g=1, and n=1
% We only want the initial launch direction without any interaction.
cfg_replay.prop  =  [0, 0, 1, 1;      % Zero voxels, air   
                     0.001, 0, 1, 1;  % Air
                     0.001, 0, 1, 1;  % Detector sillicone reflection dummy 
                     0.001, 0, 1, 1;  % 170µm glass cover slip
                     0.001, 0, 1, 1]; % Scattring medium

% run replay
[flux_sphRep, detp_sphRep, vol_sphRep, seed_sphRep]=mcxlab(cfg_replay);

% Seperate photons with z+ and z- directions
vz = detp_sphRep.v(:,3); 
% Select photons going up as there corrospond to refraction
% The photons going down corrospond to reflected photons which we discard
% as weasume only weak phase objects with small index variations
s_up = vz > 0; 

% Copy original setup for replay of only upwards launched photons 
cfg_replay_up = cfg; 
cfg_replay_up.seed = seed_sphRep.data(:,s_up);        
cfg_replay_up.detphotons= detp_sphRep.data(:,s_up);

[flux_up, detp_up, vol_up, seed_up]=mcxlab(cfg_replay_up);

imgSize = round(dimxy/3);
imgMid = round(dimxy/2);
%mcxplotvol(log(flux_up.data(imgMid-imgSize:imgMid+imgSize,imgMid-imgSize:imgMid+imgSize,1:50)));
mcxplotvol(log(flux_up.data));
lighting none

%% Save the illumination direction data v as JSON
% Select photons going up as we assume forward scattering in the focal 
% point and no reflections for small delta n. 
v = detp_sphRep.v(detp_sphRep.v(:,3) > 0,:); 
MCXData.PhotonData.v = v;
savejson('MCXData',MCXData,'FileName','LaunchAngles_Oil518FChickenMuscle.jdat','Compression','gzip')

%% Visualization of DPLD
v = detp_sphRep.v(detp_sphRep.v(:,3) > 0,:); % Select photons going up
MCXData.PhotonData.v = v;

% In order to visualize the DPLD, it must first be cast into a 2D histogram
% spanning the discrete lateral frequency space of an image. Here we just
% use 256 bins for both axis and normalize the frequency space to [-1,1]
nb = 256;
% Corrospondence arrays in units 2/lambda (here rescaled to one)
% The corrospondence arrays uxQ, uyQ mark the respective bin each photon direction
% belongs to once discretized
kBins = linspace(-1,1-2/(nb),nb);
uxQ = discretize(v(:,1), kBins); 
uyQ = discretize(v(:,2), kBins);

DPLD = zeros(numel(kBins),numel(kBins)); % Source distribution grid
% Calculate direction distributions from corrospondence arrays
for i = 1:numel(uxQ)
    if(isnan(uxQ(i)) || isnan(uyQ(i)))
        %Skip
    else
        % Obliquity factor for sperical projection in inverse units of
        % 2/lambda (here rescaled to one)
        obFac = sqrt((1)^2 - kBins(uxQ(i))^2 - kBins(uyQ(i))^2);
        obFac = obFac*~imag(obFac); % remove imaginary entries caused by discreet binning
        DPLD(uxQ(i),uyQ(i)) = DPLD(uxQ(i),uyQ(i)) + obFac; % Increase count by +1*obFac
    end
end
% Exploit symetries and construct differential S
DPLD = DPLD + flipud(DPLD)-fliplr(DPLD+flipud(DPLD)); 
DPLD = DPLD/max(DPLD,[],'all'); % Normalize

img = DPLD;
figure;
hIm=imagesc(zeros(size(img)));
set(gca,'XTick',[], 'YTick',[], 'Position',[0,0,1,1]) %Fill the window with the image
set(hIm,'CData',img);
fprintf("Total contributed Photons: %d \n",size(v,1));
daspect([1 1 1])
colormap jet
colorbar

figure;
plot(kBins,DPLD(128,:));
figform;

