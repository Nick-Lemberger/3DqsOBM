function phi = applyT3D(PGX,PGY,T3D,alpha)
FPGX = fftshift(fftn(PGX));
FPGY = fftshift(fftn(PGY));
% Wiener filtering for deconvolution with T3Dx and T3Dy
denom = (abs(rot90(T3D)).^2 + abs(T3D).^2 + alpha);
Fphi = (conj(1i*rot90(T3D)).*FPGX + conj(1i*T3D).*FPGY)./denom;
phi = real(ifftn(ifftshift(Fphi)));
end