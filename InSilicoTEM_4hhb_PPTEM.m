% 4V5F eukaryotic 80S ribosomes <--------------------------
% 4.7 MDa

% 1RYP - 20S PROTEASOME FROM YEAST  <--------------------------
% Total no. of amino acids	 = 	6386
% Molecular weight of 1RYP.pdb	 = 	704.711 kDa
% Volume of 1RYP.pdb	 = 	853257.000 Å3

% 3J9I - acidophilum 20S proteasome
% Total no. of amino acids	 = 	5978
% Molecular weight of 3j9i	 = 	657.768 kDa
% Volume of 3j9i	 = 	804858.600 Å3

% 2WRI - ribosome
% Total no. of amino acids	 = 	3229
% Molecular weight of 2WRI.pdb	 = 	365.365 kDa
% Volume of 2WRI.pdb	 = 	444916.000 Å3

% 2GTL - LUMBRICUS ERYTHROCRUORIN
% Total no. of amino acids	 = 	2419
% Molecular weight of 2GTL.pdb	 = 	275.247 kDa
% Volume of 2GTL.pdb	 = 	331225.600 Å3


% --------------------------------- InSilicoTEM
% The software simulates TEM image formation of
% (macro)molecules by taking into account the interaction potential,
% microscope aberrations and detector characteristics.
% Reference: M.Vulovic et al., Journal of Structural Biology, Volume 183, 2013, Pp 19�32
%
% Milos Vulovic 2013


close all; clear all;
% 4V5F eukaryotic 80S ribosomes <--------------------------
% 4.7 MDa

% 1RYP - 20S PROTEASOME FROM YEAST  <--------------------------
% Total no. of amino acids	 = 	6386
% Molecular weight of 1RYP.pdb	 = 	704.711 kDa
% Volume of 1RYP.pdb	 = 	853257.000 Å3

% 3J9I - acidophilum 20S proteasome
% Total no. of amino acids	 = 	5978
% Molecular weight of 3j9i	 = 	657.768 kDa
% Volume of 3j9i	 = 	804858.600 Å3

% 2WRI - ribosome
% Total no. of amino acids	 = 	3229
% Molecular weight of 2WRI.pdb	 = 	365.365 kDa
% Volume of 2WRI.pdb	 = 	444916.000 Å3

% 2GTL - LUMBRICUS ERYTHROCRUORIN
% Total no. of amino acids	 = 	2419
% Molecular weight of 2GTL.pdb	 = 	275.247 kDa
% Volume of 2GTL.pdb	 = 	331225.600 Å3

% 1SA0 - TUBULIN-COLCHICINE
% Total no. of amino acids	 = 	1934
% Molecular weight of 1SA0.pdb	 = 	216.551 kDa
% Volume of 1SA0.pdb	 = 	259272.800 Å3

% 4hhb - hemoglobin     <--------------------------------------
% Total no. of amino acids	 = 	574
% Molecular weight of 4hhb	 = 	64 kDa
% Volume of 4hhb	 = 	75960.400 Å3


% ----------------------- General processing parameters (proc field)
params.proc.N               = 1024; % Image size (field of view)
params.proc.partNum         = 50;   % Number of particles.
params.proc.geom            = 0;  % Specify orientation and translation of particles in 'PartList.m' (=1) or generate them randomly (=0)
params.proc.cores           = 40;  % the numer of matlab pools to be open for parfor loops
params.proc.gpu             = true; % use GPU acceleration
params.proc.write_hdf       = true;
params.proc.scratch_dir     = '/scratch';
params.proc.output_dir      = '/home/bueckerr/local_data';
params.proc.fn              = '4hhb_300_phase-plate_TEM_20_4.h5';

% ----------------------- Specimen (spec field)
params.spec.source          = 'pdb';        % Options: 'map', 'pdb', or 'amorph'
params.spec.pdbin        = '4hhb' ;      % Required if params.spec.source = 'map' or 'pdb' ( 2GTL, 2WRJ, 1SA0, 1RYP or any other pdb entree)
params.spec.mapsample    = '2WRJ_VoxSize1.0A.mrc'; % if input is a map the name should contain the info about voxel size with following convention '_VoxSize'%02.2f'
params.spec.potcontribution = 'iasa';    % Potential type. Options: 'iasa' or 'iasa+pb'(note: for 'iasa+pb' you first need to calcualate the pb potential via apbs. See the manual)
params.spec.motblur         = 0;            % Motion blur in [A]
params.spec.thick           = 50e-09;       % Thickness of the specimen[m].
params.spec.imagpot         = 3;            % Amplitude contrast flag. Options: (=0, none) (=1 constant Q) (=2 ice "plasmons") (=3 'plasmons" of ice and protein)
params.spec.overlap         = false;         % true does not work (yet)        

% ---------------------- Electron-specimen interaction (inter field)
params.inter.type           = 'ms'; % Options: 'pa' (projection), 'wpoa' (weak-phase), 'pa+wpoa', 'tpga'(thick-phase grating), 'ms'(multislice)
params.inter.msdz       = 1.0e-9; % Approximate thickness of the slice for multislice [m] Required if params.inter.type = 'ms'


% ---------------------- Microscope (mic field)
params.mic.Cs          = 2.7e-3;   % Spherical aberration  [m]
params.mic.a_i         = 0.050e-3; % Illumination aperture [rad]
params.mic.C_c         = 2.7e-3;   % Chromatic aberration  [m]
params.mic.deltaE      = 0.7;      % Energy spread of the source [eV]
% ----------------------- aperture
params.mic.diam_obj          = 200e-6;   % Diameter of objective aperture [m]
params.mic.foc               = 4.7e-3;   % Focal distance [m]
% ----------------------- ideal phase plate (optional)
params.mic.PPflag             = 1;        % Phase plate flag
params.mic.PP_Phase       = 0.5*pi;   % Phase plate phase shift [rad]
params.mic.PP_qcuton      = 1/300*1e10;% Cut-on frequency [1/m]. Typical values: 1/25, 1/75, 1/125, 1/250, 1/250


%----------------------- Acquisition settings (acquis field)
params.acquis.pixsize        = 1.7e-10;   % Pixel size in the specimen plane [m] 3748
params.acquis.df             = [50]*1e-9; % Defocus [m]. Note: undefocus df>0; overfocus df <0.
params.acquis.df_run = params.acquis.df;
params.acquis.ast            = [0]*1e-9;    % Astigmatism [m].
params.acquis.astangle       = 0*pi/180;    % Astigmatism angle [rad]
params.acquis.Voltage        = 300e3;       % Acceleration voltage
params.acquis.tilt           = [0]/180*pi;  % Tilt geometry -20:5:20
params.acquis.dose_on_sample = [20]/length(params.acquis.tilt); % Integrated flux (for the full tilt series) [e-/A2]

% ----------------------- Detector-Camera (cam field)
params.cam.type              = 'K2'; % Options: 'custom', 'Eagle4k', 'US4000', 'US1000GIF', 'FalconI', 'perfect' (64% at Nq -counting mode), 'ideal' (100% at Nq)
params.cam.bin               = 1; % hardware binning
% if params.cam.type = 'custom' please characterize your camera (provide MTF and DQE files and if needed readnoise and dark current).
% The characterization can be done via e.g. the tools and methods described in Vulovic et al. 2010 Acta Crist. D
% If you want to simulate MTF without accurate characterization set params.cam.GenMTFasEMG
params.cam.GenMTFasEMG = 1;
params.cam.DQEflag          = 1; % Flag for dqe (=0 means NTF = MTF)


% ---------------------- Display what? (disp field)
params.disp.generateWhat   = 'im'; % Options: 'im', 'exitw', 'imNoiseless'
params.disp.ctf            = 1; % Flag to display CTF
params.disp.mtfdqe         = 1; % Flag to display MTF and DQE

dip_randomseed(randi(9999))

% ---------------------- Parse parameters
params2 = parsePar(params);


% ---------------------- Generate and/or load 3D potential of a particle
[PartPot, params2] = loadSamples(params2);


% ---------------------- Padding and placing the particles within the volume
[InputVol, PosOrient] = generateFullVolume(PartPot,params2);

if params2.proc.write_hdf
    if exist(params2.proc.h5file,'file')
        delete(params2.proc.h5file)
    end
    
    h5create(params2.proc.h5file,'/dose',1);
    h5write(params2.proc.h5file,'/dose',params.acquis.dose_on_sample(1));
    
    h5create(params2.proc.h5file,'/dx',1);
    h5create(params2.proc.h5file,'/dz',1);
    h5write(params2.proc.h5file,'/dx',params.acquis.pixsize);
end

%% ---------------------- Image simulation
[imStructOut] = simTEM_v2(InputVol, params2);

%% ---------------------- Display
params.disp.generateWhat = 'im';
switch params.disp.generateWhat
    case 'im'
        dipshow(imStructOut.series, 'percentile')
    case 'exitw'
        dipshow(imStructOut.exit, 'lin')
    case 'imNoiseless'
        dipshow(imStructOut.noiseless_series, 'lin')
end

% --------------------- Writing the image (stack in case of a series).
WriteMRC(dip_array(imStructOut.series),params2.acquis.pixsize,[params2.proc.h5file(1:end-3) '_noisy.mrc'])
WriteMRC(dip_array(imStructOut.noiseless_series),params2.acquis.pixsize,[params2.proc.h5file(1:end-3) '_noiseless.mrc'])
