% This code was optimized for speed using Claude AI, model Opus 5 with the 
% original self-written "calc2D3DOTF.m" used as working basis. 
% Minor errors regarding slight grid and bin alingment offsets were
% discovered in this process and fixed.

function OTF = calc2D3DOTF_fast(lambda,NA,FOVx,FOVz,pxN,pzN,P,S,varargin)
%CALC2D3DOTF_FAST  Drop-in replacement for calc2D3DOTF, without EwSphereDiff.
%
%   OTF = calc2D3DOTF_fast(lambda,NA,FOVx,FOVz,pxN,pzN,P,S)
%
%   Note the signature: there is no gridRes argument any more. The spherical
%   mesh is gone; the Ewald-difference surface is evaluated directly on the
%   Cartesian u-grid.
%
% -------------------------------------------------------------------------
% WHY THIS IS FASTER
%   Substituting w = u - q/2 in Eq. (1) turns every half-pixel shift into an
%   integer one:
%
%       integrand(w) = [S(w) - S(w+q)] * P(w) * P(w+q)
%       qz(w)        = Rz(w+q) - Rz(w),      Rz(w) = sqrt(1/lambda^2 - |w|^2)
%
%   Consequences:
%     * imtranslate is no longer needed - array slicing replaces it, and the
%       bilinear half-pixel interpolation (which blurred P and S) disappears.
%     * Rz is INDEPENDENT of q, so it is precomputed once. Both the integrand
%       and the qz surface then cost two array lookups each.
%     * EwSphereDiff is not needed. Its 2*gridRes^2 mesh points were rounded
%       onto this very grid and deduplicated to one z per pixel anyway - the
%       direct evaluation gives that result exactly, with 100 % pixel coverage
%       instead of 69 % at gridRes = 512.
%     * The scalar inner loop over the surface points becomes one accumarray.
%     * Only the pupil-overlap region is touched, not the full pxN x pxN array.
%
% DIFFERENCES FROM calc2D3DOTF (all deliberate, all corrections)
%   1. No interpolation: integer shifts only. T2D will differ slightly from the
%      old version because the old one blurred P and S by half a pixel.
%   2. qz = 0 is placed at index pzN/2+1, matching MATLAB's fftshift convention
%      for even pzN. The old code used pzN/2, one bin off, which corresponds to
%      a one-voxel axial shift of the reconstruction.
%   3. The mirror index is pxN+2-i, not pxN+1-i: on the grid
%      kxBins(a) = (a-1-pxN/2)*dkx, the partner of index a is pxN+2-a. The old
%      index broke the Hermitian symmetry of the OTF by one pixel.
%   Set 'Legacy',true to reproduce 2. and 3. as in the original.
%
% OPTIONS
%   'ZInterp'  'linear'  Splat qz linearly between the two neighbouring bins
%                        instead of rounding ('nearest'). Reduces axial
%                        discretisation noise at no measurable cost.
%   'Legacy'   false     Reproduce the old index conventions (see above).
%   'Parallel' false     Use parfor over the outer loop.
%   'Verbose'  true
% -------------------------------------------------------------------------

p = inputParser;
p.addParameter('ZInterp','linear');
p.addParameter('Legacy',false);
p.addParameter('Parallel',false);
p.addParameter('Verbose',true);
p.parse(varargin{:});
o = p.Results;

dkx = 1/FOVx;                 % = 2*(0.5/(FOVx/pxN))/pxN
dkz = 1/FOVz;
r0  = 1/lambda;
Rpup = (NA/lambda)/dkx;       % pupil radius in pixels (coherent cutoff)
qmax = 2*Rpup;                % pupil-overlap limit, = 2NA/lambda

assert(mod(pxN,2)==0 && mod(pzN,2)==0, 'pxN and pzN must be even.');
assert(Rpup <= r0/dkx, ...
   'Pupil radius %.1f px exceeds the Ewald radius %.1f px - check calcPupil.', ...
    Rpup, r0/dkx);

% ---- q-independent quantities, computed once ----------------------------
a  = (1:pxN)' - 1 - pxN/2;                 % integer offset, 0 at index pxN/2+1
[AX,AY] = ndgrid(a,a);                     % AX -> rows (qx), AY -> cols (qy)
Rad2 = (AX.^2 + AY.^2)*dkx^2;
Rz = sqrt(max(r0^2 - Rad2,0));
Rz(Rad2 > r0^2) = NaN;                     % evanescent: no propagating solution
clear AX AY Rad2

Pmask = P ~= 0;
% tight bounding box of the pupil, so slices stay small
rr = find(any(Pmask,2)); cc = find(any(Pmask,1));
if o.Verbose
    fprintf('calc2D3DOTF_fast: pupil %.1f px, |q|max %.1f px, pupil bbox %dx%d\n', ...
            Rpup, qmax, numel(rr), numel(cc));
end

zOff = pzN/2 + 1;  if o.Legacy, zOff = pzN/2; end
mOff = 2;          if o.Legacy, mOff = 1;     end   % mirror index pxN+mOff-i

half = pxN/2;
T2Dq = zeros(half+1,half+1);               % includes the q = 0 row/column
T3Dq = zeros(half+1,half+1,pzN);

loopBody = @(i) deal([]);  %#ok<NASGU>  (placeholder, real body below)

if o.Parallel
    parfor i = 1:half+1
        [t2row,t3row] = local_row(i,half,pxN,pzN,qmax,S,Pmask,P,Rz,dkz,zOff,o.ZInterp);
        T2Dq(i,:)     = t2row;
        T3Dq(i,:,:)   = t3row;
    end
else
    for i = 1:half+1
        [t2row,t3row] = local_row(i,half,pxN,pzN,qmax,S,Pmask,P,Rz,dkz,zOff,o.ZInterp);
        T2Dq(i,:)     = t2row;
        T3Dq(i,:,:)   = t3row;
        if o.Verbose && mod(i,5)==0
            fprintf('  row %d of %d\n', i, half+1);
        end
    end
end

% ---- place the quarter block and mirror ---------------------------------
T2D = zeros(pxN,pxN);
T3D = zeros(pxN,pxN,pzN);
T2D(1:half+1,1:half+1)   = T2Dq;
T3D(1:half+1,1:half+1,:) = T3Dq;
clear T2Dq T3Dq

% indices 1..half+1 are computed (q <= 0); 1 is Nyquist and has no partner,
% half+1 is q = 0 and is its own partner. Only 2..half get mirrored.
src = 2:half;
mir = pxN + mOff - src;
ok  = mir >= 1 & mir <= pxN;
sm  = src(ok);  mm = mir(ok);

% qy mirror (columns of T2D): even (+).   qx mirror (rows of T2D): odd (-).
T2D(1:half+1, mm)   =  T2D(1:half+1, sm);
T3D(1:half+1, mm,:) =  T3D(1:half+1, sm,:);
T2D(mm, :)          = -T2D(sm, :);
T3D(mm, :, :)       = -T3D(sm, :, :);

OTF.T2D = T2D * dkx^2;
OTF.T3D = T3D * dkx^2 / dkz;

if o.Verbose
    % T3D is returned as a density in qz (divided by dkz), so the discrete
    % form of  int T3D dqz = T2D  is   sum(T3D,3)*dkz = T2D.
    r = sum(OTF.T3D,3) * dkz ./ OTF.T2D;
    m = abs(OTF.T2D) > 1e-2*max(abs(OTF.T2D(:)));
    fprintf('  energy check  sum(T3D,3)*dkz/T2D :  median %.6f  min %.6f  max %.6f\n', ...
            median(r(m)), min(r(m)), max(r(m)));
    fprintf('  (plotting sum(T3D,3)./T2D without the dkz gives 1/dkz = FOVz = %.1f)\n', 1/dkz);
end
end

% =========================================================================
function [t2row, t3row] = local_row(i,half,pxN,pzN,qmax,S,Pmask,P,Rz,dkz,zOff,zinterp)
t2row = zeros(1,half+1);
t3row = zeros(1,half+1,pzN);
nqx = i - 1 - pxN/2;                 % T2D row index  -> shifts COLUMNS
for j = 1:half+1
    nqy = j - 1 - pxN/2;             % T2D col index  -> shifts ROWS
    if nqx^2 + nqy^2 > qmax^2, continue; end

    a1 = max(1,1-nqy); a2 = min(pxN,pxN-nqy);    % rows  shifted by nqy
    b1 = max(1,1-nqx); b2 = min(pxN,pxN-nqx);    % cols  shifted by nqx
    if a1>a2 || b1>b2, continue; end
    r0s = a1:a2;  c0s = b1:b2;
    r1s = r0s+nqy; c1s = c0s+nqx;

    m = Pmask(r0s,c0s) & Pmask(r1s,c1s);       % pupil overlap only
    if ~any(m,'all'), continue; end

    vals = (S(r0s,c0s) - S(r1s,c1s)) .* P(r0s,c0s) .* P(r1s,c1s);
    qz   = Rz(r1s,c1s) - Rz(r0s,c0s);

    sel  = m & isfinite(qz) & vals ~= 0;
    if ~any(sel,'all'), continue; end
    v = vals(sel);  z = qz(sel)/dkz;

    t2row(j) = sum(v);

    if strcmpi(zinterp,'linear')
        z0 = floor(z);  f = z - z0;
        k0 = z0 + zOff;  k1 = k0 + 1;
        g0 = k0>=1 & k0<=pzN;   g1 = k1>=1 & k1<=pzN;
        acc = accumarray(k0(g0), v(g0).*(1-f(g0)), [pzN 1]) + ...
              accumarray(k1(g1), v(g1).*   f(g1) , [pzN 1]);
    else
        k = round(z) + zOff;
        g = k>=1 & k<=pzN;
        acc = accumarray(k(g), v(g), [pzN 1]);
    end
    t3row(1,j,:) = acc;
end
end