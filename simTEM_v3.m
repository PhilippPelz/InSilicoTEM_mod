
%simTEM simulates the full image formation given the input interaction
% potential and acquisition parameters (dose, optics and detector)
%
% SYNOPSIS:
% [imStructOut,paramsout]= simTEM(InputVol, params2)
%
% PARAMETERS:
%  InputVol: Interaction potenital specimen volume
%   params2: Structure containing various input simulation paramters
%
% OUTPUT:
% imStructOut: Structure containg output images (or stacks): noisy, noiseless and exit wave image
%              Optional: ctf and mtf (if params.disp.ctf or params.disp.mtfdqe are set)

% (C) Copyright 2013
%  Quantitative Imaging Group      Leiden University Medical Center
%  Faculty of Applied Sciences     Department of Molecular Cell Biology
%  Delft University of Technology  Section Electron Microscopy
%  Lorentzweg 1                    2300 RC Leiden
%  2628 CJ Delft
%  The Netherlands
%
%  Milos Vulovic

[psi_exit, t, dzprop] = getExitWave(InputVol, params2);

[intensity, ctf, perf_intensity] = imageExitWave(psi_exit, params2);

[series, noiseless_series] = detectElectrons(intensity,params2);

%% ---------------------------------Output structure
imStructOut.series           = series;
imStructOut.noiseless_series = noiseless_series;
imStructOut.exit             = psi_exit;
imStructOut.ctf = ctf_out;
if params2.disp.mtfdqe
    imStructOut.mtf = params2.cam.mtf;
    imStructOut.dqe = params2.cam.dqe;
end