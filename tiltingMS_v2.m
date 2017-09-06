function [OutVol, n] = tiltingMS_v2(InVol,tiltang, params2)
%tiltingMS Tilts (rotates) the input volume
% SYNOPSIS:
% [TiltVol, n] = tiltingMS(InVol,tiltang, params2)
%
% PARAMETERS:
%   InVol: input volume
% titlang: tilt angle
% params2: structure containing various input physical and processing parameters
%
% OUTPUT:
%  TiltVol: Tilted volume
%  n: number of slices for multislice

% (C) Copyright 2013
%  Quantitative Imaging Group      Leiden University Medical Center
%  Faculty of Applied Sciences     Department of Molecular Cell Biology
%  Delft University of Technology  Section Electron Microscopy
%  Lorentzweg 1                    2300 RC Leiden
%  2628 CJ Delft
%  The Netherlands
%
%  Milos Vulovic
voxSz = params2.acquis.pixsize;

if tiltang ~= 0
    TiltVol = rotation(dip_image(InVol), tiltang, 2, 'bspline','zero');
    %fix bug :rotation.m function doesn't work with complex input volume. Rotate 2 times (real and imag) or ask Bernd?
    TiltVol(TiltVol < min(InVol)) = min(InVol);
    TiltVol = dip_array(cut(TiltVol, [size(InVol,1), size(InVol,2), size(TiltVol,3)]));
else
    TiltVol = InVol;
end

thicknessfull = voxSz * size(InVol,3) / cos(tiltang); % thickness of tilted volume
n = ceil(thicknessfull / params2.inter.msdz); % to ensure the integer number of slices

if params2.inter.msdz > thicknessfull
    n = 1;
    fprintf('The slice thickness can not be larger than specimen thickness. Use only one slice...\n')
end

% the next part makes the voxel number and slice number commensurate
if mod(single(thicknessfull/voxSz), n) > 0
    intslices = thicknessfull/voxSz - mod(single(thicknessfull/voxSz), n) + n;
else
    intslices = thicknessfull/voxSz - mod(single(thicknessfull/voxSz), n);
end

% pad the output Volume such that its thickness becomes commensurate with
% the number of slices, if necessary (disabled, as this is taken care of
% now in the slicing outside this function)
if false %size(TiltVol,3) ~= intslices
    fprintf('The potential array z size is incommensurate with the slice thickness. Consider changing.')
    OutVol = zeros(size(InVol,1), size(InVol,1), intslices,'like',TiltVol);
    TV_pos = floor((size(OutVol,3) - size(TiltVol,3))/2) + 1;
    
    OutVol(:,:,TV_pos:(TV_pos + size(TiltVol,3)-1)) = TiltVol;
    
else
    
    OutVol = TiltVol;
    
end

