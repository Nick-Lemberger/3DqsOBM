function P = calcPupil(lambda,NA,FOVx,pxN,beta)
% Calculate the pupil function for a gaussian beam arriving at the back
% focal apparture of an microscope objective.
% beta defines the ratio between the 1/e^2 beam size w and the back focal
% apparture diameter BFA of the microscope objective. 
% beta =  BFA / w
% beta = 0 defines a gaussian light source with infinite size,
% essentially plane waves at the back focal aparture. 

% Calculate lateral frequency grids
dx = FOVx/pxN; % Image pixel size in µm
kx_max = 0.5/dx; % Max spatial frequency in lines/µm after fft of image
dkx = 2*kx_max/pxN; % Frequency step per pixel in lines/µm of image fft

% CORRECTION: removed factor 2 
fc_NA = NA/lambda; % Cut-off frequency from MO in lines/µm 
mask_fc = fc_NA/dkx;


% Calculate pupil function
P = fspecial('disk', mask_fc) == 0;
% Check if the disk is larger then the grid
if floor(pxN-mask_fc)<0 % If yes, cut the disk to the grid size
    mask_pxR = size(P,1);
    P = imcrop(P,[(mask_pxR/2-pxN/2), ...
        (mask_pxR/2-pxN/2) pxN-1 pxN-1]);
    clear('mask_pxR');
    fprintf("Image frequency space smaller then NA! Use more pixels!");
else  % If not, padd array to match the mask to the grid size
    P = imresize(padarray(P, [floor((pxN/2)-mask_fc) ...
        floor((pxN/2)-mask_fc)], 1, 'both'), [pxN pxN]);
end
P = ~P; % invert

if(beta > 0)
    sigma = 0.5 * fc_NA / beta; % Size of gaussian distribution
    kxBins = linspace(-kx_max,kx_max-dkx,pxN);
    [kx,ky] = meshgrid(kxBins,kxBins);
    beam_pupil = exp( (-(kx.^2)-(ky.^2)) ./ (2.*sigma.^2));
    P = P.* beam_pupil;
end

end