%% Load the image z-stack slice by slice and preprocess
% Modify according to your file and dataformat. 
% This script also includes an extra SHG channel. Comment out if not needed. 

clear('PGX_norm','PGY_norm','PGX_data','PGY_data','SUM_data' ... 
    ,'im1save','im2save','im3save','im4save','i','samplerate' ...
    ,'PGX','PGY');
load(sprintf('parameters.mat'));
% parameters.mat contains:
% Zpositions: 1xN double array containing sample Z positions in µm
% Zstack_size: size of Zstack
% galvo_amp: amplitude of galvometer, specific to your scanners. We always
%            pecify the lateral field of view to an amplitude of 1. 
% pixel_nr: number of lateral pixels, square images are assumed
load(sprintf('Zstack_IMG%d.mat',1));

% Used to cut out the image edges in order to remove scanning artifacts
% that could disturb the reconstruction. We cut 1024px images to 900px.
ROI_cutout = 900;
ogImsize = length(im1save);
if ROI_cutout < ogImsize
    imSize = ROI_cutout;
else
    imSize = ogImsize;
end
Zstack_nr = length(Zpositions);

% Read in the images
SUM = zeros(imSize, imSize, Zstack_nr);
PGX_raw = zeros(imSize, imSize, Zstack_nr);
PGY_raw = zeros(imSize, imSize, Zstack_nr);
%SHG = zeros(imSize, imSize, Zstack_nr);
for i = 1:Zstack_nr
    fprintf("Loading slice %d of %d!\n",i,Zstack_nr);
    load(sprintf('Zstack_IMG%d.mat',i));
    if ROI_cutout < ogImsize        
        PGX_raw(:,:,i) = imcrop(-im1save,[(ogImsize/2-imSize/2) (ogImsize/2-imSize/2) imSize-1 imSize-1]);
        PGY_raw(:,:,i) = imcrop(-im2save,[(ogImsize/2-imSize/2) (ogImsize/2-imSize/2) imSize-1 imSize-1]);
        SUM(:,:,i) = imcrop(abs(im3save),[(ogImsize/2-imSize/2) (ogImsize/2-imSize/2) imSize-1 imSize-1]);
        %SHG(:,:,i) = imcrop(abs(im4save),[(ogImsize/2-imSize/2) (ogImsize/2-imSize/2) imSize-1 imSize-1]);
    else       
        PGX_raw(:,:,i) = -im1save;
        PGY_raw(:,:,i) = -im2save;
        SUM(:,:,i) = abs(im3save);
        %SHG(:,:,i) = abs(im4save);
    end
end
    
G = 4.34;
% The SUM amplifier summs each photodiode channel with factor 1/4 to 
% prevent saturation. G compensates for this.
% However, we modulate and detect at 1MHz, which is almost at the bandwidth 
% limit of the SUM amplifier attenuating the signal somewhat and leading 
% to an effective factor of 1/4.34 
PGX(:,:,:) = PGX_raw(:,:,:) ./ (G.*SUM(:,:,:)); % Normalization of images
PGY(:,:,:) = PGY_raw(:,:,:) ./ (G.*SUM(:,:,:)); % Normalization of images

% Optional. Subtract a heavily smoothed copy to remove possible background gradients
for i = 1:Zstack_nr   
    PGX(:,:,i) = PGX(:,:,i) - imgaussfilt(PGX(:,:,i),imSize/5);
    PGY(:,:,i) = PGY(:,:,i) - imgaussfilt(PGY(:,:,i),imSize/5);
end

fprintf('Z-stack loaded!\n');
clear('PGX_raw','PGY_raw','im1save','im2save','im3save','im4save','i','samplerate','Z0');

figure;
orthosliceViewer(SUM);
colormap gray
figure;
orthosliceViewer(PGX);
colormap gray
figure;
orthosliceViewer(PGY);
colormap gray
%figure;
%orthosliceViewer(SHG);
%colormap gray
%% Optional rescaling for up/downsampling 
% to resize the lateral frequency space. Ths is needed if you subsample 
% your image, i.e, your resolution is below the resolution of your MOs NA.
% If this happens, the image frequency space would be smaller than the
% pupil function of your MO and might lead to errors in the following
% reconstruction. 
scaleFactor = 2;
PGX_new = zeros(imSize*scaleFactor, imSize*scaleFactor, Zstack_nr);
PGY_new = zeros(imSize*scaleFactor, imSize*scaleFactor, Zstack_nr);
SUM_new = zeros(imSize*scaleFactor, imSize*scaleFactor, Zstack_nr);
%SHG_new = zeros(imSize*scaleFactor, imSize*scaleFactor, Zstack_nr);
for i = 1:Zstack_nr       
    PGX_new(:,:,i) = imresize(PGX(:,:,i),scaleFactor);
    PGY_new(:,:,i) = imresize(PGY(:,:,i),scaleFactor);
    SUM_new(:,:,i) = imresize(SUM(:,:,i),scaleFactor);
    %SHG_new(:,:,i) = imresize(SHG(:,:,i),scaleFactor);
end
PGX = PGX_new;
PGY = PGY_new;
SUM = SUM_new;
%SHG = SHG_new;
imSize = imSize*scaleFactor;
clear('PGX_new','PGY_new','SUM_new','SHG_new');
figure;

orthosliceViewer(PGX);
colormap gray
%%
lambda = 1.050; % Wavelength in µm
NA = 0.8; % NA of microscope objective
FOVperGA = 276.93; 
% Field of view for a galvoamplitude of 1. 
% 281.4µm for SEWIA
% 387.7µm for Nikon 50x NA=0.6
% 276.93µm for Nikon 50x NA=0.8
% (These are FOV values for our MOs for our personal reference.
% You need to specify your own FOV for a certain galvo amplitude)

pxN = imSize; % Lateral pixel number
pzN = Zstack_nr; % Axial pixel number
galvo_amp_eff = galvo_amp * ROI_cutout/ogImsize; % Effective galvoamplitude if cropped by ROI
FOVx = FOVperGA*galvo_amp_eff; % Field of view in µm
FOVz = Zstack_size; % Axial field of view in µm

% Calculate lateral frequency grids
dx = FOVx/pxN; % Image pixel size in µm
kx_max = 0.5/dx; % Max spatial frequency in lines/µm after fft of image
dkx = 2*kx_max/pxN; % Frequency step per pixel in lines/µm of image fft

% CORRECTION removed factor 2*NA -> 1*NA
fc_NAp = NA/lambda; % Cut-off frequency from coherent illumination MO in lines/µm
mask_fc = fc_NAp/dkx;

% Calculate axial frequency axis
dz = FOVz/pzN; % Image pixel size in µm
kz_max = 0.5/dz; % Max spatial frequency in lines/µm after fft of image
dkz = 2*kz_max/pzN; % Frequency step per pixel in lines/µm of image fft

kxBins = linspace(-kx_max,kx_max-dkx,pxN); % lateral image frequency space k axis
kzBins = linspace(-kz_max,kz_max-dkz,pzN); % axial image frequency space k axis

% Add path to include the extra functions like calcPupil() etc.
addpath('G:\User\Nick Lemberger\Promotion\4Q-Epi');

% Calculate pupil function
% beta =  BFA / w, ratio between the 1/e^2 beam size w and the back focal
% apparture diameter BFA
% Use beta to model the Gaussian laser beam profile being clipped by the 
% back focal apparture of your MO. Use beta = 0 for a flat profile, i.e.
% plane waves. 
beta = 0.51; 
P = calcPupil(lambda,NA,FOVx,pxN,beta);
img = P;
figure;
hIm=imagesc(zeros(size(img)));
set(gca,'XTick',[], 'YTick',[], 'Position',[0,0,1,1]) %Fill the window with the image
set(hIm,'CData',img);
daspect([1 1 1])
colormap gray
colorbar

% Plot the profile to check if it fits your beam shape after the back focal
% apparture. 
%figure;
%plot(P(:,451));
%figform;

% Calculate the source function (detected photon launch distribution, DPLD)
% from a simulation. Simulate the distribution with the file
% "DPLD_MXC_simulation.m" using MXCLAB.
sim_data_path = 'G:\User\Nick Lemberger\Promotion\4Q-Epi\Simulation\MCX_Simulations\LaunchAngles_Oil518FGlassAgar.jdat';
%sim_data_path = 'G:\User\Nick Lemberger\Promotion\4Q-Epi\Simulation\MCX_Simulations\LaunchAngles_Oil518FGlassChickenMuscle.jdat';
S = calcSource(lambda,FOVx,pxN,P,sim_data_path);
S = imgaussfilt(S, pxN/200, 'Padding', 'symmetric'); % Smooth to reduce MC simulation noise

[KX,KY] = meshgrid(kxBins,kxBins);
bandMask = sqrt(KX.^2 + KY.^2) <= 1/lambda; % Incoherent DPC cutoff
S = S.*bandMask; % Cut possible frequency leakage caused by the smoothing

img = S;
figure;
hIm=imagesc(zeros(size(img)));
set(gca,'XTick',[], 'YTick',[], 'Position',[0,0,1,1]) %Fill the window with the image
set(hIm,'CData',img);
daspect([1 1 1])
colormap jet
colorbar
%%
% Calculate T2D and T3D

% Unoptimized old version. Very slow and not recommended for use,
% but easier to read and comprehend.
%gridRes = 1024;
%OTF = calc2D3DOTF(lambda,NA,FOVx,FOVz,pxN,pzN,gridRes,P,S);

% Much faster and highly optimized version
OTF = calc2D3DOTF_fast(lambda,NA,FOVx,FOVz,pxN,pzN,P,S);

T2D = OTF.T2D;
T3D = OTF.T3D;

% Directly save the OTF for later use, so you dont have to recalculate
% every time. 
save('OTF_justCalced.mat', 'T2D', 'T3D')

img = T2D;
figure;
hIm=imagesc(zeros(size(img)));
set(gca,'XTick',[], 'YTick',[], 'Position',[0,0,1,1]) %Fill the window with the image
set(hIm,'CData',img);
daspect([1 1 1])
colormap jet
colorbar 

img = sum(T3D,3);
figure;
hIm=imagesc(zeros(size(img)));
set(gca,'XTick',[], 'YTick',[], 'Position',[0,0,1,1]) %Fill the window with the image
set(hIm,'CData',img);
daspect([1 1 1])
colormap jet
colorbar

% Residual error between T2D and summed T3D. Should be 1 everywhere.
img = sum(T3D,3)./ T2D ./ pzN*2;
figure;
hIm=imagesc(zeros(size(img)));
set(gca,'XTick',[], 'YTick',[], 'Position',[0,0,1,1]) %Fill the window with the image
set(hIm,'CData',img);
daspect([1 1 1])
colormap jet
colorbar

figure;
orthosliceViewer(permute(-T3D,[2,1,3]));
colormap jet
colorbar

%% Optional prefiltering. Remove stuff outside of NA
[KX,KY] = meshgrid(kxBins,kxBins);
bandMask = sqrt(KX.^2 + KY.^2) <= 1*NA/lambda; % DPC cutoff
for slice = 1:Zstack_nr
   PGX_ = PGX(:,:,slice);
   PGY_ = PGY(:,:,slice);
   %SHG_ = SHG(:,:,slice);
   % Apply NA of illumination MO as mask to remove out of band noise
   F_PGX_ = fftshift(fft2(PGX_)).* bandMask;
   F_PGY_ = fftshift(fft2(PGY_)).* bandMask;
   %F_SHG_ = fftshift(fft2(SHG_)).* (P>0);
   
   % Remove some specific noise bands but be careful.
   %F_PGX_(:,451) = 0;
   %F_PGY_(:,451) = 0;
   %F_PGX_(451,:) = 0;
   %F_PGY_(451,:) = 0;
   %F_PGX_(:,612:613) = 0;
   %F_PGY_(:,612:613) = 0;
   %F_PGX_(:,289) = 0;
   %F_PGY_(:,289) = 0;
   
   PGX(:,:,slice) = real(ifft2(ifftshift(F_PGX_)));
   PGY(:,:,slice) = real(ifft2(ifftshift(F_PGY_)));
   %SHG(:,:,slice) = real(ifft2(ifftshift(F_SHG_)));
end
clear('PGX_','PGY_','SHG_','F_PGX_','F_PGY_','F_SHG_');

%% Reconstruct Z-Stack with T3D

% Applies the Tikhonov regularization for direct de-convolution and arrives
% automatically at a suitable alpha via generalized cross-validation.
[V_T3D alpha] = applyT3D_GCV3D(PGX,PGY, lambda/(4*pi) * T3D,'Mask','data','MaxPoints',1e12,'Plot',false);

% Here you can set alpha manually instead:
%alpha =1e-4;
%V_T3D = applyT3D(PGX,PGY, lambda/(4*pi)*T3D,alpha);

% You need to set n0 manually. For a true quantitative absolute
% reconstruction, n0 must be known precisely. For a relative 
% reconstruction, n0 can be set approximately as the relative error 
% is very small for small deviations of n0. 
n0 = 1.5037; % Refractive index 518F immersion oil
%n0 = 1.45; % Approximate average index of refraction chicken fat tissue
%n0 = 1.36 % Approximate average index of brain tissue

k = 2*pi*n0 / lambda;
n = n0 - n0/(2*k^2) * V_T3D;

figure;
orthosliceViewer(n);
caxis([1.45, 1.57])
colormap gray
colorbar

