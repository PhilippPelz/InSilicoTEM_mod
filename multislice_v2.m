function psi_multislice  = multislice_v2(psi_t,Nm,n,lambda,q_m,dzprop,show,use_gpu)
%multislice Performs multislice calculations 
% SYNOPSIS:
% psi_multislice  = multislice(psi_t,Nm,n,lambda,q_m,dzprop)
%
% PARAMETERS:
%   psi_t: the phase grating 
%      Nm: the sieze of the image
%       n: number of slices
%  lambda: wavelength of the incoming electrons
%      qm: frequencies in the Fourier domain
%  dzprop: propagation distance of the Fresnel propagator
%
% OUTPUT:
%  psi_multislice: electron wave exiting the specimen

% (C) Copyright 2013
%  Quantitative Imaging Group      Leiden University Medical Center
%  Faculty of Applied Sciences     Department of Molecular Cell Biology
%  Delft University of Technology  Section Electron Microscopy
%  Lorentzweg 1                    2300 RC Leiden
%  2628 CJ Delft
%  The Netherlands
%
%  Milos Vulovic
%
% modified in 2017 by Robert Bücker, MPSD Hamburg

if nargin < 7, show = false; end
if nargin < 8, use_gpu = false; end

% Free Fresnel propagator
P = exp(-1i*pi * lambda * q_m.^2 * dzprop); % Fresnel propagator

% multislice wave function
psi_multislice = ones(Nm, Nm, 'like', 1+1i);

% overload dip_image definition
ft = @(I) fftshift(fftn(ifftshift(I))) ./ sqrt(numel(I));
ift = @(I) ifftshift(ifftn(fftshift(I))) .* sqrt(numel(I));

if show
    figure(803), clf
end

if use_gpu
    P = gpuArray(P);
    psi_multislice = gpuArray(psi_multislice);
end

for ii = 1:n
    fprintf('Slice %d out of %d\n', ii, n)    
    psi_multislice = ift( ft(psi_multislice .* psi_t(:,:,ii)) .* P);
    if show
        imagesc(angle(psi_multislice))
        drawnow
    end
    %psi_multislice = ift(ft(psi_multislice*squeeze(psi_t(:,:,ii-1)))*P(dzprop));
end

if use_gpu
    psi_multislice = gather(psi_multislice);
end
