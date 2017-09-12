function [psi_exit, t, dzprop] = getExitWave(InputVol, params2)

voxSz = params2.acquis.pixsize; % the voxel size

% prealocate memory for the stack (series)
nTiltAngles = length(params2.acquis.tilt);

thickness = voxSz*size(InputVol,3);

switch params2.inter.type
    case{'pa+wpoa', 'pa', 'wpoa', 'tpga'}
        error('In this version, only multislice is supported (for now)')
        if params2.proc.cores>1
            parpool(params2.proc.cores)
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
            parpool close
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
        
    case 'ms'
        % if multislice then the exit wave is calculated separately for each tilt angle.
        
        psi_exit = zeros(params2.proc.N, params2.proc.N, nTiltAngles, 'like', 1+1i);
        
        for ll = 1:nTiltAngles
            tiltang = params2.acquis.tilt(ll);
            fprintf('Simulate tilt angle %3.0f\n', tiltang * 180/pi)
            
            % potential post-processing for MS: commensurability with slice
            % number, and tilt
            [potential, n] = tiltingMS_v2(dip_array(InputVol), tiltang, params2);
            
            sizepot = size(potential);
            Nm      = max(sizepot(1), sizepot(2));
            thicknessfull = sizepot(3)*voxSz;
            dzprop = thicknessfull / n;
            
            % modify potential by constructing imaginary part
            % ...could this be moved to t instead of potential?
            if params2.spec.imagpot == 0
                % potential = potential; % do nothing
            elseif params2.spec.imagpot == 1
                potential = potential +1i*params2.spec.potenampl*potential;
            elseif params2.spec.imagpot == 2
                potential = potential +1i*params2.spec.potenampl*(newim(potential)+1);
            elseif params2.spec.imagpot == 3
                potential = potential + ...
                    1i * (params2.spec.potenampl*(potential==0) + params2.spec.proteinampl*(potential>0));
            else
                error('This option for the params.spec.potenampl is not valid. Please choose between 0-3');
            end
            
            % preallocation of potential array for multislice
            t = zeros(Nm,Nm,n,'like',1+1i);
            
            % projected potential within slices (phase grating)
            zpix_per_slice = ceil(sizepot(3)/n);
            zpix_off = round((n*zpix_per_slice - sizepot(3)) / 2);
            fprintf('Building slice potentials (%d slices, %d voxels thickness)\n', n, zpix_per_slice)
            for ii = 1:n
                zrng = ((ii-1)*zpix_per_slice+1 : ii*zpix_per_slice) - zpix_off;
                zrng(zrng < 1) = []; zrng(zrng > sizepot(3)) = [];
                t(:,:,ii) = sum(potential(:,:,zrng),3);
            end
            t = t / zpix_per_slice;
            
            xwm = voxSz * Nm;           % pixelsize for multislice * size sample
            q_true_pix_m = 1 / xwm;     % Fourier domain
            q_m = dip_array(rr([Nm Nm])) * q_true_pix_m; % 2D frequencies in Fourier domain
            
            % dipshow(dip_image(imag(psi_t(:,:,30))));
            
            % run actual multislice.
            % as first argument, it requires the potential propagator, not
            % the slice potential itself!
            psi_exit(:,:,ll) = multislice_v2(exp(1i * params2.inter.sig_transfer * t * dzprop), ...
                Nm, n, params2.inter.lambda, q_m, dzprop, ...
                false, params2.proc.gpu);
            % original version for troubleshooting
            %psi_exit(:,:,ll) = dip_array(multislice(dip_image(psi_t), Nm, n, params2.inter.lambda, q_m, dzprop));
            
            if strcmp(params2.spec.source, 'amorph')
                thicknessfull = params2.spec.thick/cos(tiltang);
                psi_exit(:,:,ll) = psi_exit(:,:,ll) * exp(-params2.inter.sig_transfer * params2.spec.potenampl * thicknessfull);
            end
            
            proj_psi_t = exp(1i*params2.inter.sig_transfer*sum(t,3)*dzprop);
            %psi_exit_angle = angle(proj_psi_t); %WTF?!
            psi_exit_angle = angle(psi_exit(:,:,ll));
            fprintf('max ang %2.3f\n', max(max(psi_exit_angle)));
            figure
            imagesc(psi_exit_angle);
            colorbar;
            title('Exit angle (multislice)')
            figure
            imagesc(angle(proj_psi_t));
            colorbar;
            title('Projected potential')
            figure
            imagesc(psi_exit_angle - angle(proj_psi_t));
            colorbar;
            title('Difference')
            
            % export stuff if desired...
            if params2.proc.write_hdf
                h5create(params2.proc.h5file,'/vreal',size(t),'Deflate',9,'ChunkSize',[10 10 10]);
                h5create(params2.proc.h5file,'/vimag',size(t),'Deflate',9,'ChunkSize',[10 10 10]);
                h5write(params2.proc.h5file,'/vreal',real(t));
                h5write(params2.proc.h5file,'/vimag',imag(t));
                h5write(params2.proc.h5file,'/dz',dzprop);
            end
        end
end