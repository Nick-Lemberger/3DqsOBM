function OTF = calc2D3DOTF(lambda,NA,FOVx,FOVz,pxN,pzN,gridRes,P,S)
% Calculate lateral frequency grids
dx = FOVx/pxN; % Image pixel size in µm
kx_max = 0.5/dx; % Max spatial frequency in lines/µm after fft of image
dkx = 2*kx_max/pxN; % Frequency step per pixel in lines/µm of image fft

fc_NA = 2*NA/lambda; % Cut-off frequency from MO in lines/µm
mask_fc = fc_NA/dkx;

% Calculate axial frequency axis
dz = FOVz/pzN; % Image pixel size in µm
kz_max = 0.5/dz; % Max spatial frequency in lines/µm after fft of image
dkz = 2*kz_max/pzN; % Frequency step per pixel in lines/µm of image fft

kxBins = linspace(-kx_max,kx_max-dkx,pxN); % image frequency space k axis

r0 = 1/lambda;
% Resolution of the grid used to calculate the T3D in sperical coordinates
phiN = gridRes;
thetaN = gridRes;

T2D = zeros(pxN,pxN);
T3D = zeros(pxN,pxN,pzN);

for i = 1:pxN/2
    tic
    for j = 1:pxN/2
        qx = kxBins(i); % Select a qx
        qy = kxBins(j); % Select a qy
        % Only continue if shifted pupils would overlap as the contribution
        % is zero otherwise
        if(sqrt((qx/2)^2 + (qy/2)^2) <= (ceil(mask_fc))*dkx)
            P_plus  = imtranslate(P,  [qx/2/dkx, qy/2/dkx]); % Shift pupil
            P_minus = rot90(P_plus,2); % rotate 180° due to symmetry
            S_plus  = imtranslate(S,  [qx/2/dkx, qy/2/dkx]); % Shift illumination distribution
            S_minus = -rot90(S_plus,2); % rotate 180° and negate due to symmetry
            SP = (S_plus-S_minus).*P_plus.*P_minus; % Integrand without delta function (T2D)
            T2D(i,j) = sum(SP,'all'); % Integrate SP to calculate T2D at (qx,qz)
            % Calculate surface defined by the delta function of the T3D
            % Returned is a list with x,y,z indecies of the points on the surface.
            ew = EwSphereDiff([qx,qy],r0,phiN,thetaN,dkx,dkz,pxN,pzN);
            % Only a quarter of the T3D is calculated directly in the x,y
            % plane by looping only up to pxN/2.
            % The rest is constructed according to its symmetry.
            for m = 1:length(ew(:,3))
                if(ew(m,3) >= 1 && ew(m,3) <= pzN)
                    T3D(i,j,ew(m,3))             = T3D(i,j,ew(m,3))             +SP(ew(m,1),ew(m,2));
                    T3D(pxN+1-i,j,ew(m,3))       = T3D(pxN+1-i,j,ew(m,3))       -SP(ew(m,1),ew(m,2));
                    T3D(i,pxN+1-j,ew(m,3))       = T3D(i,pxN+1-j,ew(m,3))       +SP(ew(m,1),ew(m,2));
                    T3D(pxN+1-i,pxN+1-j,ew(m,3)) = T3D(pxN+1-i,pxN+1-j,ew(m,3)) -SP(ew(m,1),ew(m,2));
                end
            end
        end
    end
    toc
    fprintf("Done with sub-integration %d of %d!\n",i,pxN/2);
end
% Calculate the rest of T2D via mirroring according to its symmetry
T2D = T2D + fliplr(T2D);
T2D = T2D - flipud(T2D);
% Finally, multiply by integration step size.
OTF.T2D = T2D * dkx^2; 
OTF.T3D = T3D * dkx^2;
end