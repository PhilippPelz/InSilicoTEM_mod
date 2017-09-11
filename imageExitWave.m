
function [intensity, ctf, perf_intensity] = imageExitWave(psi_exit, params2)

nFocSeries = length(params2.acquis.df);
nPsiExit = size(psi_exit,3);
nRuns = nFocSeries * nPsiExit;

intensity = zeros(params2.proc.N,params2.proc.N,nRuns);
ctf = dip_image(zeros(params2.proc.N,params2.proc.N,nFocSeries));

for k = 1:nFocSeries
    params2.df_run = params2.df(k);
    fprintf('Calculating CTF for df = %g\n', params2.df_run)
    ctf(:,:,k) = simulateCTF(params2);
end
if params2.mic.PPflag
    ctf = 0.9 * ctf;
end

if nargout > 2
    
    perf_intensity = zeros(params2.proc.N,params2.proc.N,nRuns);
    
    perf = params2;
    perf.mic.Cs                = 0;   % Spherical aberration  [m]
    perf.mic.C_c               = 0;   % Chromatic aberration  [m]
    perf.mic.deltaE            = 0.0;      % Energy spread of the source [eV]
    perf.mic.diam_obj          = 200e-6;   % Diameter of objective aperture [m]
    perf.mic.foc               = 4.7e-3;   % Focal distance [m]
    perf.acquis.df             = [0]*1e-9; % Defocus [m]. Note: undefocus df>0; overfocus df <0.
    perf.acquis.df_run         = perf.acquis.df;
    perf.acquis.ast            = [0]*1e-9;    % Astigmatism [m].
    perf.acquis.astangle       = 0*pi/180;    % Astigmatism angle [rad]
    perf.mic.PPflag            = 1;
    perfectctf = simulateCTF(perf);
    
    if perf.mic.PPflag
        perfectctf = 0.9 * perfectctf;
    end
    
end

for k = 1:nFocSeries
    for jjj = 1:nPsiExit
        
        set = (jjj-1) * nFocSeries + k;
        
        fprintf('Applying CTF to set %d out of %d', set, nRuns);
        
        % this could be made much more elegant (blockwise with bsxfun)... but for now...
        btot = ctf(:,:,k)*dip_fouriertransform(dip_image(psi_exit(:,:,jjj)),'forward',[1 1 ]); %0]);
        
        intensity(:,:,set) = dip_array(abs(dip_fouriertransform(btot,'inverse',[1 1])).^2); % intensity in the image without camera influence
        
        if nargout > 2
            perfbtot = perfectctf*dip_fouriertransform(dip_image(psi_exit(:,:,jjj)),'forward',[1 1 ]); %0]);
            perf_intensity(:,:,jjj) = dip_array(abs(dip_fouriertransform(perfbtot,'inverse',[1 1])).^2); % intensity in the image without camera influence
        end
        
    end
end

end

