function p = EwSphereDiff(q,r0,phiN,thetaN,dkx,dkz,pxN,pzN)
% The delta function in the integral of the 3DOTF describes a surface in 3D 
% frequency space which is defined in x,y by the insersection of two
% shifted semi-spheres with hermitian symmetry by(qx/2,qy/2) and 
% (-qx/2,-qy/2), respectively, and in z by the subtraction of the 
% intersecting part. 
% A list of points is returned in positive pixel corrdinates that can be
% used as index to directly access all allowed u for integration. 
phiBins = linspace(0,2*pi,phiN);
thetaBins = linspace(0,pi/2,thetaN);
[phi, theta] = meshgrid(phiBins,thetaBins); % phi and theta meshgrid

% Convert sperical surface mesh to cartesian coordinates
xSphere = r0 * sin(theta) .* cos(phi);
ySphere = r0 * sin(theta) .* sin(phi);
xSphere = xSphere(:); % linear indexing
ySphere = ySphere(:); % linear indexing
% xSphere, ySphere: x and y coordinates of a sperical surface

% Calculate shifted ewald spheres
Zplus =  sqrt(r0^2 - (xSphere+q(1)/2).^2 - (ySphere+q(2)/2).^2);
Zminus = sqrt(r0^2 - (xSphere-q(1)/2).^2 - (ySphere-q(2)/2).^2);

% Calculate the difference, find imaginary contributions
% Zplus./abs(Zplus)+Zminus./abs(Zminus)-2 makes sure that contributions in
% the form of ia * ib are also removed by marking them as imaginary
Zdiffplane = (Zplus-Zminus) + Zplus./abs(Zplus) + Zminus./abs(Zminus) - 2;
Zdiffplane = Zdiffplane(:);
ind_realz = imag(Zdiffplane)==0; % Select only values that are truly real

% Select the corrosponding coordinates, convert to pixel coordinates 
x = xSphere(ind_realz)/dkx;
y = ySphere(ind_realz)/dkx;
z = Zdiffplane(ind_realz)/dkz;

% Round to nearest gridpoint in cartesian pixel coordinates
x = round(x);
y = round(y);
z = round(z);

% Make sure that all selected points are in bounds of the grid.
inBounds = (x>=-pxN/2+1 & x<=pxN/2) & (y>=-pxN/2+1 & y<=pxN/2) & (z>=-pzN/2+1 & z<=pzN/2);
x = x(inBounds);
y = y(inBounds);
z = z(inBounds);

% Convert to linear coordinates for easier removal of duplicate 
% coordinates after approximating to cartesian grid
% Important for energy conservation: 
% Each XYZ point must be unique and each point must also be unique in 
% the XY plane, meaning each point in the XY plane has only one Z value.
linCoords = (x + pxN/2) + (y + pxN/2)*(pxN+1) + (z + pzN/2)*(pxN+1)^2;
% Remove multiple apperances of the same pixel caused by the rounding
% Also remove possible NaNs
linCoords = rmmissing(unique(linCoords));
% Get back to 3D coordinates
x = round(mod(linCoords, pxN+1));
linCoords = floor(linCoords/(pxN+1));
y = round(mod(linCoords, pxN+1));
linCoords = floor(linCoords/(pxN+1));
z = round(linCoords);

% Set up buffers for fast Z averaging
% Buffer and count matrix used for summing and averaging of Z values of 
% duplicate points in the XY plane.
zSumBuff = zeros(pxN); 
zCount = zeros(pxN);
for i = 1:length(x)
    zSumBuff(x(i),y(i)) = zSumBuff(x(i),y(i)) + z(i); % Summ Z values
    zCount(x(i),y(i)) = zCount(x(i),y(i)) + 1; % Count entries per pixel
end
zSumBuff = round(zSumBuff./zCount); % Calculate mean Z value for each pixel
for i = 1:length(x)
    z(i) = zSumBuff(x(i),y(i)); % get the mean values out
end

% Convert again to linear coordinates for unique() operation as this is the
% fastest way to remove duplicates
linCoords = x + y*(pxN+1) + z*(pxN+1)^2;       
% Remove multiple apperances of the same pixel caused by the rounding
% Also remove possible NaNs
linCoords = rmmissing(unique(linCoords));
% Get back to 3D coordinates
x = round(mod(linCoords, pxN+1));
linCoords = floor(linCoords/(pxN+1));
y = round(mod(linCoords, pxN+1));
linCoords = floor(linCoords/(pxN+1));
z = round(linCoords);

p = [x,y,z];
end



% % New arrays to save the unique coordinates. Prefill with NaN for
% % simple removal of unsused space later with  rmmissing(). 
% x_uniqZ = zeros(length(linCoords),1) + NaN;
% y_uniqZ = zeros(length(linCoords),1) + NaN;
% z_uniqZ = zeros(length(linCoords),1) + NaN;
% % Slower-than-i moving index for following for-loop
% idc_uniqZ = 1;
% % Used to mark and skip duplicates already taken care of
% remainingDupes = ones(length(linCoords),1);
% 
% for i = 1:length(linCoords_zFlat)
%     if(remainingDupes(i))
%         k = linCoords_zFlat(i); % Look at a coordinate with removed Z
%         indc_mask = linCoords_zFlat == k; % And get index mask of all duplicates
%         remainingDupes(indc_mask) = 0; % Mark all duplicates of k as found
%         % Get the full XYZ linear coordinate of the identified XY duplicates
%         sel_linCoords = linCoords(indc_mask); 
%         
%         % Calculate the 3D components. Calculate the mean value for Z and
%         % round to nearest grid point.
%         x_uniqZ(idc_uniqZ) = round(mod(sel_linCoords(1), pxN+1));
%         
%         sel_linCoords = floor(sel_linCoords / (pxN+1));
%         y_uniqZ(idc_uniqZ) = round(mod(sel_linCoords(1), pxN+1));
%         
%         sel_linCoords = floor(sel_linCoords / (pxN+1));
%         z_uniqZ(idc_uniqZ) = round(sum(sel_linCoords)/length(sel_linCoords));
%         
%         idc_uniqZ = idc_uniqZ+1;
%     end
% end