% This code was checked for errors and optimized for speed using 
% Claude AI, model Opus 5.

function [V, alphaOpt, info] = applyT3D_GCV3D(PGX, PGY, T3D, varargin)
%   applyT3D with the Tikhonov parameter chosen by GCV.
%
%   Drop-in replacement for applyT3D: same argument order, same conventions,
%   same output, except that alpha is determined from the
%   data by generalized cross validation (Ledwig & Robles 2019, Eqs. 12-13).
%
% -------------------------------------------------------------------------
% NAME-VALUE OPTIONS
%   'Mask'       'support'  Which frequencies enter the GCV sums:
%                           'support' - where S > 0    (recommended)
%                           'all'     - the whole 3D volume
%                           'data'    - where S > 0 AND the data is non-zero;
%                                       use this if the PG stacks were
%                                       pre-filtered with (P>0) in qsOBM_3D.m
%                           or your own logical [Ny Nx Nz] array.
%   'AlphaRange' []         [aMin aMax]; default 10.^[-12 2] * max(S).
%   'NumAlpha'   200
%   'Refine'     true       Golden-section refinement in log10(alpha).
%   'MaxPoints'  2e7        Cap on the number of frequencies entering the GCV
%                           sums (random, fixed seed). Affects only the alpha
%                           search, never the reconstruction. Inf disables.
%   'n0'         []         With 'Lambda', converts V to refractive index.
%   'Lambda'     []
%   'Scale'      1          Calibration constant c multiplying V (see NOTE 1).
%   'Verbose'    true
%   'Plot'       false
%
% OUTPUTS
%   V         Scattering potential, or n if 'n0' and 'Lambda' are given.
%   alphaOpt  GCV-selected regularisation parameter.
%   info      .alphaGrid .gcv .alphaGridMin .atBoundary .gcvAtOpt
%             .nFreqTotal .nFreqUsed .subsampled .residual .traceTerm
%
%
p = inputParser;
p.addParameter('Mask',      'support');
p.addParameter('AlphaRange', [],   @(x)isempty(x)||numel(x)==2);
p.addParameter('NumAlpha',   200,  @(x)isscalar(x)&&x>=5);
p.addParameter('Refine',     true);
p.addParameter('MaxPoints',  2e7,  @(x)isscalar(x)&&x>=1e3);
p.addParameter('n0',         [],   @(x)isempty(x)||isscalar(x));
p.addParameter('Lambda',     [],   @(x)isempty(x)||isscalar(x));
p.addParameter('Scale',      1,    @isscalar);
p.addParameter('Verbose',    true);
p.addParameter('Plot',       false);
p.parse(varargin{:});
o = p.Results;

sz = size(PGX);
assert(numel(sz)==3,          'PGX must be a 3D z-stack.');
assert(isequal(size(PGY),sz), 'PGX and PGY must have identical size.');
assert(isequal(size(T3D),sz), 'T3D must have the same size as the PG stacks.');
assert(isreal(T3D), ['T3D is expected to be real, as produced by calc2D3DOTF. ' ...
                     'Do not pre-multiply by 1i - that is done internally.']);
assert(isempty(o.n0)==isempty(o.Lambda), ...
       'Supply both ''n0'' and ''Lambda'', or neither.');

T3Drot = rot90(T3D);                  % rotates in the first two dimensions
FPGX = fftshift(fftn(PGX));
FPGY = fftshift(fftn(PGY));

M = 2;
nFreqTotal = numel(T3D);

% accumulate the GCV statistics plane by plane
Sc = cell(sz(3),1);
Ec = cell(sz(3),1);
Nc = cell(sz(3),1);
maxS = 0;
for k = 1:sz(3)
    Tr = T3Drot(:,:,k);
    Tt = T3D(:,:,k);
    Fx = FPGX(:,:,k);
    Fy = FPGY(:,:,k);

    Sk = Tr.^2 + Tt.^2;            % |T_x|^2+|T_y|^2 (the 1i drops out)
    Nk = Tr.*Fx + Tt.*Fy;          % |N|^2 unaffected by the common -1i
    Ek = abs(Fx).^2 + abs(Fy).^2;

    if ischar(o.Mask) || isstring(o.Mask)
        switch lower(char(o.Mask))
            case 'all',     mk = true(size(Sk));
            case 'support', mk = Sk > 0;
            case 'data',    mk = (Sk > 0) & (Ek > 0);
            otherwise, error('Unknown ''Mask'' option ''%s''.', char(o.Mask));
        end
    else
        mk = o.Mask(:,:,k);
    end

    Sc{k} = double(Sk(mk));
    Ec{k} = double(Ek(mk));
    Nc{k} = double(abs(Nk(mk)).^2);
    maxS  = max(maxS, max(Sk(:)));
end
Sv = vertcat(Sc{:});
Ev = vertcat(Ec{:});
Nv = vertcat(Nc{:});
clear Sc Ec Nc Tr Tt Fx Fy Sk Nk Ek mk

Nq = numel(Sv);
assert(Nq > 0, 'No frequencies selected - check the ''Mask'' option.');

subsampled = false;
if Nq > o.MaxPoints
    rs  = RandStream('mt19937ar','Seed',42);
    idx = randperm(rs, Nq, round(o.MaxPoints));
    Sv = Sv(idx);  Ev = Ev(idx);  Nv = Nv(idx);
    Nq = numel(Sv);
    subsampled = true;
end
n = M*Nq;

% GCV sweep
if isempty(o.AlphaRange)
    assert(maxS > 0, 'T3D is identically zero.');
    aRange = maxS * [1e-12, 1e2];
else
    aRange = sort(double(o.AlphaRange(:)'));
end
aGrid = logspace(log10(aRange(1)), log10(aRange(2)), round(o.NumAlpha));

gcv = zeros(size(aGrid));
for k = 1:numel(aGrid)
    gcv(k) = local_G(aGrid(k), Sv, Ev, Nv, M, Nq, n);
end
[~, kMin]  = min(gcv);
atBoundary = (kMin==1) || (kMin==numel(aGrid));
alphaOpt   = aGrid(kMin);

if o.Refine && ~atBoundary
    lo = log10(aGrid(kMin-1));  hi = log10(aGrid(kMin+1));
    gr = (sqrt(5)-1)/2;
    c1 = hi-gr*(hi-lo);   c2 = lo+gr*(hi-lo);
    f1 = local_G(10^c1,Sv,Ev,Nv,M,Nq,n);
    f2 = local_G(10^c2,Sv,Ev,Nv,M,Nq,n);
    for it = 1:80
        if f1 < f2
            hi=c2; c2=c1; f2=f1;
            c1=hi-gr*(hi-lo);  f1=local_G(10^c1,Sv,Ev,Nv,M,Nq,n);
        else
            lo=c1; c1=c2; f1=f2;
            c2=lo+gr*(hi-lo);  f2=local_G(10^c2,Sv,Ev,Nv,M,Nq,n);
        end
        if (hi-lo) < 1e-8, break; end
    end
    alphaOpt = 10^((lo+hi)/2);
end

if atBoundary
    warning('qsOBM_GCV3D:Boundary', ...
        ['GCV minimum sits on the edge of [%.3g %.3g]. Widen ''AlphaRange''. ' ...
         'A minimum pinned at the LOWER bound usually means the data were ' ...
         'pre-filtered or smoothed, so GCV sees no noise left to suppress.'], ...
         aRange(1), aRange(2));
end

% reconstruct, identical to applyT3D.m
denom = T3Drot.^2 + T3D.^2 + alphaOpt;
Fphi  = (conj(1i*T3Drot).*FPGX + conj(1i*T3D).*FPGY) ./ denom;
clear FPGX FPGY denom T3Drot
V = real(ifftn(ifftshift(Fphi))) * o.Scale;
clear Fphi

if ~isempty(o.n0)
    kk = 2*pi*o.n0 / o.Lambda;
    V  = o.n0 - o.n0*V./(2*kk^2);         % manuscript Eq. (3); see NOTE 3
end

% -------------------------------- report -----------------------------------
info = struct('alphaGrid',aGrid, 'gcv',gcv, 'alphaGridMin',aGrid(kMin), ...
              'atBoundary',atBoundary, 'nFreqTotal',nFreqTotal, ...
              'nFreqUsed',Nq, 'subsampled',subsampled);
[info.residual, info.traceTerm] = local_parts(alphaOpt,Sv,Ev,Nv,M,Nq);
info.gcvAtOpt = local_G(alphaOpt,Sv,Ev,Nv,M,Nq,n);

if o.Verbose
    if subsampled, tag = ' (random subset)'; else, tag = ''; end
    fprintf('qsOBM_GCV3D: %dx%dx%d, %d of %d frequencies used%s\n', ...
            sz(1),sz(2),sz(3), Nq, nFreqTotal, tag);
    fprintf('             alpha_GCV = %.6e   (grid minimum %.6e)\n', ...
            alphaOpt, aGrid(kMin));
end

if o.Plot
    figure; loglog(aGrid,gcv,'-','LineWidth',1.4); hold on
    loglog(alphaOpt,info.gcvAtOpt,'o','MarkerSize',8,'LineWidth',1.4); grid on
    xlabel('$\alpha$','Interpreter','latex');
    ylabel('$V(\alpha)$','Interpreter','latex');
    legend({'GCV functional','selected $\alpha$'}, 'Interpreter','latex', ...
           'Location','northwest');
end
end

% ==========================================================================
function [num, tr] = local_parts(a, Sv, Ev, Nv, M, Nq)
den = Sv + a;
num = sum( Ev - Nv .* (Sv + 2*a) ./ den.^2 );
tr  = M*Nq - sum( Sv ./ den );
end

function g = local_G(a, Sv, Ev, Nv, M, Nq, n)
[num, tr] = local_parts(a, Sv, Ev, Nv, M, Nq);
g = (num/n) / (tr/n)^2;
end