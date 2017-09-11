
function [series, noiseless_series] = detectElectrons(intensity)

nDose = length(params2.acquis.dose);
nIntensity = size(intensity,3);
nRuns = nDose * nIntensity;

series = zeros([size(intensity,1), size(intensity,2), nRuns]);
noiseless_series = zeros(size(series));

for iii= 1:nIntensity
    for k = 1:nDose
        
        set = (jjj-1) * nDose + k;
        
        fprintf('Applying dose setting %d (out of %d) to input image %d (out of %d)',k,iii,nDose,nIntensity);
        
        params2.influx = params2.acquis.dose(k);
        [IntenDetect1, IntenNoiseless2] = DetectSim(squeeze(dip_image(intensity(:,:,iii))), params2);
        
        series(:,:,set) = dip_array(IntenDetect1);
        noiseless_series(:,:,set) = dip_array(IntenNoiseless2);
        
        if params2.proc.write_hdf
            h5create(params2.proc.h5file,strcat('/img_noiseless_',num2str(iii)),size(noiseless_series(:,:,iii)),'Deflate',9,'ChunkSize',[20 20]);
            h5create(params2.proc.h5file,strcat('/img_noisy_',num2str(iii)),size(series(:,:,iii)),'Deflate',9,'ChunkSize',[20 20]);
            
            h5write(params2.proc.h5file,strcat('/img_noiseless_',num2str(iii)),double(IntenNoiseless2));
            h5write(params2.proc.h5file,strcat('/img_noisy_',num2str(iii)),double(IntenDetect1));
        end
        
    end
end

end