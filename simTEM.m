function imStructOut= simTEM(InputVol, params2)
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


voxSz = params2.acquis.pixsize;% the voxel size
% prealocate memory for the stack (series)
nTiltAngles = length(params2.acquis.tilt);

if strcmp(params2.seriesout,'tilt')
    Nseries= nTiltAngles;
elseif strcmp(params2.seriesout,'defocus')
    Nseries= length(params2.acquis.df);
elseif strcmp(params2.seriesout,'dose')
    Nseries= length(params2.acquis.dose_on_sample);
else 
    Nseries=1;
end
series        = newim(params2.proc.N,params2.proc.N,Nseries); 
btot_i        = newim(params2.proc.N,params2.proc.N,Nseries);
extprojstack  = newim(params2.proc.N,params2.proc.N,Nseries, 'complex');

poten0 = InputVol;
thickness = voxSz*size(poten0,3);

switch params2.inter.type
    case{'pa+wpoa', 'pa', 'wpoa', 'tpga'}
        if params2.proc.cores>1
           matlabpool(params2.proc.cores)
        end
        poten0 = permute(poten0,[1 3 2]);
    % calculate the projected potential (neglects the thickness of the specimen)
    switch params2.inter.type
        case {'pa+wpoa','pa'}
            poten0_stack = FreqDomRadon_2d(poten0,params2.acquis.tilt);
    % in case of wpoa and tpga the Ewald sphere (parabola) is sampled and propagation through the specimen thickness is taken into account 
        case {'wpoa','tpga'}
            poten0_stack = FreqDomRadon_Ewald_2d(poten0,params2.acquis.tilt, params2);
    end
    if params2.proc.cores>1
        matlabpool close
    end
    poten0_stack = permute(poten0_stack,[1 3 2]);
    poten0_stack = cut(poten0_stack,[size(InputVol,1) size(InputVol,2) size(poten0_stack,3)]);
    
    % imaginary part of the potential (amplitude contrast)
    if params2.spec.imagpot == 1
       poten_stack = poten0_stack + 1i*params2.spec.potenampl*poten0_stack;
    elseif params2.spec.imagpot == 2 || params2.spec.imagpot == 3
       imagPotProj = double((newim(poten0_stack)+1)*dip_image(reshape(params2.spec.potenampl*round(thickness/(voxSz))./cos(params2.acquis.tilt), [1 1 nTiltAngles])));
       poten_stack = poten0_stack + 1i*imagPotProj;
    else
       poten_stack = poten0_stack;
    end
    
    % transmission function
    if strcmp(params2.inter.type, 'pa')|| strcmp(params2.inter.type, 'tpga') 
        psi_exit   = exp(1i*params2.inter.sig_transfer*poten_stack*voxSz);
    % the first term (corresponds to the single scattering)    
    elseif strcmp(params2.inter.type, 'pa+wpoa')|| strcmp(params2.inter.type, 'wpoa')
        psi_exit   = 1+1i*params2.inter.sig_transfer*poten_stack*voxSz;
    end
    extprojstack = psi_exit;
    
% if multislice then the exit wave is calculated separately for each tilt angle.  
case 'ms' 
    psi_exit = newim(params2.proc.N, params2.proc.N, nTiltAngles, 'dcomplex');
    %if strcmp(params2.seriesout,'tilt')
    for ll=1:nTiltAngles
        tiltang = params2.acquis.tilt(ll);
        fprintf('Simulate tilt angle %3.0f\n', tiltang*180/pi)     
        [potenext, n] = tiltingMS(poten0,tiltang, params2); 
        sizepot = size(potenext);
        Nm      = max(sizepot(1), sizepot(2));
        thicknessfull = sizepot(3)*voxSz;

        % construct imaginary part 
        if params2.spec.imagpot == 0
            potenfull = potenext;
        elseif params2.spec.imagpot == 1
            potenfull = potenext +1i*params2.spec.potenampl*potenext;
        elseif params2.spec.imagpot == 2 
            potenfull = potenext +1i*params2.spec.potenampl*(newim(potenext)+1);
        elseif params2.spec.imagpot == 3
            imagPot = params2.spec.potenampl*(potenext==0) + params2.spec.proteinampl*(potenext>0);
            potenfull = potenext +1i*imagPot; 
        else
            error('This option for the params.spec.potenampl is not valid. Please choose between 0-3');
        end
        
        % preallocation for multislice
        t  = newim([Nm, Nm, n],'dcomplex');

        % projected potential within slices (phase grating)
        for ii = 0:n-1
            t(:,:,ii) = mean(potenfull(:,:,ii*(sizepot(3)/n):(ii+1)*(sizepot(3)/n)-1),[],3);
        end

        xwm = (voxSz)*Nm;%pixelsize for multislice * size sample
        %Fourier domain
        q_true_pix_m = 1/xwm;
        q_m = rr([Nm Nm])*q_true_pix_m; % frequencies in Fourier domain

        % propagator function Fresnel Propagation 
        dzprop = thicknessfull/n;
        % transmission functions; the slice thickness is constant
        psi_t = exp(1i*params2.inter.sig_transfer*t*dzprop); 
        
        psi_tr = dip_array(real(t));
        psi_ti = dip_array(imag(t));
        h5create(params2.fn,'/vreal',size(psi_tr),'Deflate',9,'ChunkSize',[10 10 10]);
        h5create(params2.fn,'/vimag',size(psi_ti),'Deflate',9,'ChunkSize',[10 10 10]);
        h5write(params2.fn,'/vreal',psi_tr);
        h5write(params2.fn,'/vimag',psi_ti);
%         h5write(params2.fn,'/vimag',psi_ti);
        h5write(params2.fn,'/dz',dzprop);
        
%         dipisosurface(dip_image());
%         dipshow(dip_image(imag(psi_t(:,:,30))));
        % multislice
        PsiExit = multislice(psi_t, Nm, n, params2.inter.lambda, q_m, dzprop); 
        psi_exit(:,:,ll-1) = PsiExit;
        if strcmp(params2.spec.source, 'amorph')
            thicknessfull = params2.spec.thick/cos(tiltang);
            psi_exit(:,:,ll-1) = psi_exit(:,:,ll-1)*exp(-params2.inter.sig_transfer*params2.spec.potenampl*thicknessfull);
        end

        psi = psi_tr + 1j* psi_ti;
        proj_psi_t = exp(1i*params2.inter.sig_transfer*sum(psi,3)*dzprop); 
        psi_exit_angle = angle(proj_psi_t);
        fprintf('max ang %2.3f\n', max(max(psi_exit_angle)));
        imagesc(psi_exit_angle);
        colorbar;
        
        %projpot_r = resample(squeeze(real(psi_exit(:,:,ll-1))),[(voxSz*1e-10/params2.acquis.pixsize) (voxSz*1e-10/params2.acquis.pixsize)]);
        %projpot_i = resample(squeeze(imag(psi_exit(:,:,ll-1))),[(voxSz*1e-10/params2.acquis.pixsize) (voxSz*1e-10/params2.acquis.pixsize)]);
        %proj_pot  = projpot_r+1i*projpot_i; % projected potential
        %%extproj   = extend_exit(squeeze(proj_pot),params2);
        extprojstack(:,:,ll-1) = squeeze(psi_exit(:,:,ll-1));
        
    end
end
if strcmp(params2.seriesout,'defocus') || strcmp(params2.seriesout,'dose')
   psi_exit = repmat(psi_exit,[1 1 Nseries]);
end

%% ---------------------------------- CTF with df, ast, envelopes and optionally phase plate

psi_exit=double(psi_exit);   
btot_i=double(btot_i);
perfbtot_i=double(btot_i);
for jjj= 1:Nseries
    if strcmp(params2.seriesout,'defocus')
        params2.acquis.df_run=params2.acquis.df(jjj);
    else
        params2.acquis.df_run=params2.acquis.df(1);
    end
    if ~mod(jjj,5)||~mod(jjj,Nseries)
        fprintf(['Calculate the CTF for the ' params2.seriesout sprintf(' series. Image number %3d of %3d\n',  jjj, Nseries)]);
    end
  
  [ctf] = simulateCTF(params2);
  if params2.mic.PPflag 
    ctf = 0.9 * ctf;
  end 
  params2.mic.Cs                  = 0;   % Spherical aberration  [m]
    params2.mic.C_c               = 0;   % Chromatic aberration  [m]
    params2.mic.deltaE            = 0.0;      % Energy spread of the source [eV]
    params2.mic.diam_obj          = 200e-6;   % Diameter of objective aperture [m]
    params2.mic.foc               = 4.7e-3;   % Focal distance [m]
    params2.acquis.df             = [0]*1e-9; % Defocus [m]. Note: undefocus df>0; overfocus df <0. 
    params2.acquis.df_run = params2.acquis.df;
    params2.acquis.ast            = [0]*1e-9;    % Astigmatism [m]. 
    params2.acquis.astangle       = 0*pi/180;    % Astigmatism angle [rad]
    params2.mic.PPflag             = 1;
    [perfectctf] = simulateCTF(params2);
  if params2.mic.PPflag 
    perfectctf = 0.9 * perfectctf;
  end 
  if params2.disp.ctf 
    ctf_out(:,:,jjj)=double(ctf);
%     perfect_ctf_out(:,:,jjj)=double(perfectctf);
  end
    figure;
    imagesc(double(abs(ctf)));
    colorbar;   
  btot = ctf*dip_fouriertransform(dip_image(psi_exit(:,:,jjj)),'forward',[1 1 ]); %0]);
  figure;
    imagesc(log10(double(abs(btot))));
    colorbar;    
  btot_i(:,:,jjj) = double(abs(dip_fouriertransform(btot,'inverse',[1 1])).^2); % intensity in the image without camera influence              
  
  perfbtot = perfectctf*dip_fouriertransform(dip_image(psi_exit(:,:,jjj)),'forward',[1 1 ]); %0]);
  perfbtot_i(:,:,jjj) = double(abs(dip_fouriertransform(perfbtot,'inverse',[1 1])).^2); % intensity in the image without camera influence              

end
noiseless_tilt_series = dip_image(btot_i); 
noiseless_tilt_series_no_aberration = perfbtot_i; 

            
%%  --------------------------------- Camera influence          
noiseless_tilt_series=double(noiseless_tilt_series);
series = double(series);
for iii= 1:Nseries
    if strcmp(params2.seriesout,'dose')
        params2.influx=params2.acquis.dose(iii);
    else
        params2.influx=params2.acquis.dose(1);
    end
    if ~mod(iii,5)||~mod(iii,Nseries)
        fprintf(['Calculate the DQE for the ' params2.seriesout sprintf(' series. Image number %3d of %3d\n', iii, Nseries)]);
    end
    [IntenDetect1,~] = DetectSim(squeeze(dip_image(noiseless_tilt_series(:,:,iii))), params2); 
    [~,IntenNoiseless2] = DetectSim(squeeze(dip_image(noiseless_tilt_series_no_aberration(:,:,iii))), params2); 
    series(:,:,iii) = double(IntenDetect1);

    h5create(params2.fn,strcat('/img_noiseless_',num2str(iii)),size(noiseless_tilt_series(:,:,iii)),'Deflate',9,'ChunkSize',[20 20]);
    h5create(params2.fn,strcat('/img_noisy_',num2str(iii)),size(series(:,:,iii)),'Deflate',9,'ChunkSize',[20 20]);
    
    h5write(params2.fn,strcat('/img_noiseless_',num2str(iii)),double(IntenNoiseless2));    
    h5write(params2.fn,strcat('/img_noisy_',num2str(iii)),double(IntenDetect1));
end
series = dip_image(series);

 %% ---------------------------------Output structure
imStructOut.series           = series;
imStructOut.noiseless_series = IntenNoiseless2;
imStructOut.exit             = extprojstack;
if params2.disp.ctf 
   imStructOut.exit = ctf_out;
end
if params2.disp.mtfdqe
   imStructOut.mtf = params2.cam.mtf;
   imStructOut.dqe = params2.cam.dqe;
end