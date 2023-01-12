% light sheet simulation and stripe pattern deconvolution simulation
% 
% Author: Xiongtao Ruan 


codeRT = fileparts(which(mfilename));
addpath(genpath(codeRT));


%% generate light sheets

% user can define the root path for the results, if it is not defined, the
% results will be saved in LLS_simulation outside the code directory.

rt = '/clusterfs/fiona/Data/LLS_simulation_202208/';
if ~exist(rt, 'dir')
    rt = sprintf('%s/../LLS_simulation_results/', codeRT);
    mkdir(rt);
end

b = load('plane wave sets for common lattices.mat');

figure_names = {'Fig_MB_Sq45_NA0p3', 'Fig_MB_Hex_NA0p43', 'Fig_AC_Sq45_NA0p3_NAsigma0p09', 'Fig_AC_Hex_NA0p4', ...
                'Fig_MB_Sq45_NA0p31', 'Fig_AC_Sq45_NA0p3_NAsigma0p07', 'Fig_AC_AxialSW_NA0p3', 'Fig_MB_Hex_S2_NA0p43', ...
                'Fig_HB_Hex_NA0p5', 'Fig_HB_HexRect_NA0p5', 'Fig_Sinc_lateralSW_NA0p32', 'Fig_AC_rectangular_NA0p25', ...
                'Fig_Confined_AxialSW_NA0p45', 'Fig_MB_Sq45_NP0p41', 'Fig_FST_Hex_NA0p4', 'Fig_FST_Sq_NA0p3', ...
                'Fig_Gaussian_NA0p21_NAsimga0p21', 'Fig_Gaussian_NA0p21_NAsimga0p42', 'Fig_Sinc_lateralSW_NA0p4', ...
                'Fig_MB_Sq45_NA0p445_Crop0p05', 'Fig_MB_Sq45_NA0p445_Crop0p1', 'Fig_MB_Sq45_NA0p495_Crop0p2', ...
                };

% user can define a valid index between 1 - numel(figure_names), if not, it
% will popromt to ask for a valid lattice light sheet index. 
lls_ind = -1;

if lls_ind <= 0
    disp('Light sheet names and indices:')
    
    nLLS = numel(figure_names);

    for i = 1 : nLLS
        fprintf('%2d  %s\n', i, figure_names{i})
    end
    fprintf('\n')

    x = input(sprintf('Please choose a light sheet to simulate with the index (1 - %d): ', nLLS));
    if ~isa(x, 'double') || (x < 1 || x > nLLS)
        error('Input must be a number from 1 to %d!', nLLS);
    end
    lls_ind = x;
end

figure_name = figure_names{lls_ind}

% default algorithm
algorithm = 'original';
xyPol = [1, 1, 0];
crop_factor = 0.1;  
switch figure_name
    case 'Fig_MB_Sq45_NA0p3'
        lightsheet_descrip = 'multi-Bessel';
        PW =  b.PW_Square;
        algorithm = 'multi-Bessel'; 
        NAlattice = .30;
        NAsigma = 5.0;
        NAannulus = [0.38, 0.22];
    case 'Fig_MB_Hex_NA0p43'
        lightsheet_descrip = 'multi-Bessel';        
        PW = b.PW_Hex;
        algorithm = 'multi-Bessel';
        NAlattice = .43;
        NAsigma = 3.0;
        NAannulus = [0.47, 0.40];                
        crop_factor = 0.08;        
    case 'Fig_AC_Sq45_NA0p3_NAsigma0p09'
        lightsheet_descrip = 'Sq45';
        PW =  b.PW_Square;
        NAlattice = .30;
        NAsigma = 0.09;
        NAannulus = [0.4, 0.2];
    case 'Fig_AC_Hex_NA0p4'
        lightsheet_descrip = 'Hex';
        PW = b.PW_Hex; 
        NAlattice = .40;
        NAsigma = 0.075;
        NAannulus = [0.6, 0.2];        
    case 'Fig_MB_Sq45_NA0p31'
        lightsheet_descrip = 'multi-Bessel';
        PW =  b.PW_Square;
        algorithm = 'multi-Bessel'; 
        NAlattice = .31;
        NAsigma = 3.0;
        NAannulus = [0.4, 0.3];
    case 'Fig_AC_Sq45_NA0p3_NAsigma0p07'
        lightsheet_descrip = 'Sq45';
        PW =  b.PW_Square;
        NAlattice = .30;
        NAsigma = 0.07;
        NAannulus = [0.4, 0.2]; 
    case 'Fig_AC_AxialSW_NA0p3'
        lightsheet_descrip = 'AxialSW';
        PW =  b.PW_AxialStandingWave;
        NAlattice = .30;
        NAsigma = 0.10;
        NAannulus = [0.4, 0.2];
        crop_factor = 0.02;        
    case 'Fig_MB_Hex_S2_NA0p43'
        lightsheet_descrip = 'multi-Bessel';        
        PW =  b.PW_Hex_sqrt2;   
        algorithm = 'multi-Bessel';
        NAlattice = .43;
        NAsigma = 3.0;
        NAannulus = [0.47, 0.4];        
    case 'Fig_HB_Hex_NA0p5'
        lightsheet_descrip = 'Hex';        
        PW = b.PW_Hex;
        algorithm = 'harmonic balanced';    
        NAlattice = .50;
        NAsigma = 0.075;
        NAannulus = [0.6, 0.4]; 
        crop_factor = 0.01;         
    case 'Fig_HB_HexRect_NA0p5'
        lightsheet_descrip = 'HexRect';
        PW =  b.PW_HexRect;
        algorithm = 'harmonic balanced';            
        NAlattice = .50;
        NAsigma = 0.15;
        NAannulus = [0.6, 0.4];
        crop_factor = 0.01;
    case 'Fig_Sinc_lateralSW_NA0p32'
        lightsheet_descrip = 'ExpSinc';        
        PW =  b.PW_LateralStandingWave; 
        algorithm = 'new';
        NAlattice = .32;
        NAsigma = 5.0;
        NAannulus = [0.4, 0.2];        
    case 'Fig_AC_rectangular_NA0p25'
        lightsheet_descrip = 'ExpAxialSW';        
        PW =  b.PW_AxialStandingWave;
        NAlattice = .25;
        NAsigma = 0.13;
        NAannulus = [0.35, 0.20];
        crop_factor = 0.02;        
    case 'Fig_Confined_AxialSW_NA0p45'
        lightsheet_descrip = 'AxialSW';                
        PW =  b.PW_AxialStandingWave;
        NAlattice = .45;
        NAsigma = 0.065;
        NAannulus = [0.6, 0.2];  
        crop_factor = 0.05;
    case 'Fig_MB_Sq45_NP0p41'
        lightsheet_descrip = 'multi-Bessel';                
        PW =  b.PW_Square;
        algorithm = 'multi-Bessel';
        NAlattice = .41;
        NAsigma = 3.;
        NAannulus = [0.6, 0.4];  
        crop_factor = 0.1;
    case 'Fig_FST_Hex_NA0p4'
        lightsheet_descrip = 'FST';                
        PW =  b.PW_Hex;
        algorithm = 'cosine-sinc';
        NAlattice = .40;
        NAsigma = 0.1;
        NAannulus = [0.435, 0.365];  
        crop_factor = 0.1;        
    case 'Fig_FST_Sq_NA0p3'
        lightsheet_descrip = 'FST';                
        PW =  b.PW_Square;
        algorithm = 'cosine-sinc';
        NAlattice = .30;
        NAsigma = 0.1;
        NAannulus = [0.375, 0.225];  
        crop_factor = 0.1;                
    case 'Fig_Gaussian_NA0p21_NAsimga0p21'
        lightsheet_descrip = 'ExpGaussian';                
        PW =  b.PW_Gaussian;
        algorithm = 'new';
        NAlattice = .21;
        NAsigma = 0.21;
        NAannulus = [0.4, 0.2];  
        crop_factor = 0.1;
    case 'Fig_Gaussian_NA0p21_NAsimga0p42'
        lightsheet_descrip = 'ExpGaussian';                
        PW =  b.PW_Gaussian;
        algorithm = 'new';
        NAlattice = .21;
        NAsigma = 0.42;
        NAannulus = [0.6, 0.2];  
        crop_factor = 0.1;
    case 'Fig_Sinc_lateralSW_NA0p4'
        lightsheet_descrip = 'ExpSinc';        
        PW =  b.PW_LateralStandingWave; 
        algorithm = 'new';
        NAlattice = .4;
        NAsigma = 5.0;
        NAannulus = [0.6, 0.2];        
    case 'Fig_MB_Sq45_NA0p445_Crop0p05'
        lightsheet_descrip = 'multi-Bessel';                
        PW =  b.PW_Square;
        algorithm = 'multi-Bessel';
        NAlattice = .445;
        NAsigma = 5.;
        NAannulus = [0.55, 0.44];  
        crop_factor = 0.05;
    case 'Fig_MB_Sq45_NA0p445_Crop0p1'
        lightsheet_descrip = 'multi-Bessel';                        
        PW =  b.PW_Square;
        algorithm = 'multi-Bessel';
        NAlattice = .445;
        NAsigma = 5.;
        NAannulus = [0.55, 0.44];  
        crop_factor = 0.1;
    case 'Fig_MB_Sq45_NA0p495_Crop0p2'
        lightsheet_descrip = 'multi-Bessel';                                
        PW =  b.PW_Square;
        algorithm = 'multi-Bessel';
        NAlattice = .495;
        NAsigma = 5.;
        NAannulus = [0.55, 0.44];  
        crop_factor = 0.2;
end

NAdet = 1.0;
xz_det_PSF = b.xz_PSF_RW_510nm_NA1p0;
xz_det_OTF = b.xz_OTF_RW_510nm_NA1p0;
gamma = .5;
yplot_range = 99;
ystepsize_in = 3;
visualize =~false;
pixsize = 0.1;

outputdir = [rt, 'simulations/'];
mkdir(outputdir);
subfolder = ['NAlattice', num2str(NAlattice, '%1.2f'), filesep, lightsheet_descrip, filesep, 'NAAnnulusMax', num2str(NAannulus(1), '%1.2f'), filesep, 'NAsigma', num2str(NAsigma, '%1.2f'), filesep];

% figure MB_Hex_NA0p43 and fig MB_Hex_S2_NA0p43 have the same subfolder name 
% based on their parameters the subfolder defined above is used for figure MB_Hex_NA0p43; 
% figure MB_Hex_S2_NA0p43 use a different subfolder name.
if strcmp(figure_name, 'Fig_MB_Hex_S2_NA0p43')
    subfolder = 'Fig_MB_Hex_PW_Hex_Fund_MaxComp_sqrt2_NAexc0p43_annulus0p47-0p40_crop0p10/';
end

resultPath = sprintf('%s/%s/decon_simulation/', outputdir, subfolder)
mkdir(resultPath);

llsFn = sprintf('%s/PSF_OTF_simulation.mat', resultPath);

PSF_folder = sprintf('%sPSFs/', resultPath);
mkdir(PSF_folder);

figurePath = sprintf('%sfigures/', resultPath);
mkdir(figurePath);


%% run LLS simulation

if ~false
    fprintf('Start lattice light sheet simulation for %s...\n', figure_name)
    tic
    [yFWHM, PWb, OnAxisIntensity, CumulativeFWHM, Overall_Dithered_xz_OTF2, DitheredxzPSFCrossSection, ...
        PlotSLMPhase, MaskIntensity, SampleIntensityAtFocus, DitheredIntensity, OverallCrossSection, ...
        SheetCrossSection, DetPSFCrossSection] = Calc_and_Plot_3D_LLS_PSFs_and_OTFs_v9e_cleaned(algorithm, ...
        lightsheet_descrip, xyPol, PW, NAlattice, NAsigma, NAannulus, crop_factor, NAdet, xz_det_PSF, ...
        xz_det_OTF, gamma, yplot_range, ystepsize_in, outputdir, subfolder);
    toc
    save('-v7.3', llsFn, 'yFWHM', 'PWb', 'OnAxisIntensity', 'CumulativeFWHM', ...
        'Overall_Dithered_xz_OTF2', 'DitheredxzPSFCrossSection', ...
        'PlotSLMPhase', 'MaskIntensity', 'SampleIntensityAtFocus', 'DitheredIntensity', ...
        'OverallCrossSection', 'SheetCrossSection', 'DetPSFCrossSection', 'yplot_range', ...
        'ystepsize_in'); 
    fprintf('Lattice light sheet simulation for %s is done!\n', figure_name)
end



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% theoretical 3D PSF generation, stripe pattern generation and decon simulation

% decon simulation parameters
% snr for the stripe patterns
snrs = [20];
snr_ref = 20;

generatePSF = true;
generateData = true;


%% load 3D detection PSF and DitheredxzPSFCrossSection 

if generatePSF
    fn = [codeRT, '/RW_PSFs/', 'Det_PSF_OTF_3D_510_NA1p0_px_38p346nm_RichardsWolf.mat'];
    if ~exist(fn, 'file')
        error(['The 3D detection PSF file %s does not exist, please follow the instruction', ...
            'in github repo to download the file or check if it is in the right location!'], fn);
    end
    a = load(fn, 'PSF_RW_510nm_NA1p0');
    PSF_RW_510nm_NA1p0 = a.PSF_RW_510nm_NA1p0;
    
    a = load(llsFn, 'DitheredxzPSFCrossSection', 'yFWHM');
    DitheredxzPSFCrossSection = a.DitheredxzPSFCrossSection;
    yFWHM = a.yFWHM
    clear a
else
    a = load(llsFn, 'yFWHM');
    yFWHM = a.yFWHM
    clear a
end

% center psf index 
cind = 1;
cme = [cind];


%% generate center/middle/edge PSFs

if generatePSF
    fprintf('Generate 3D theoretical PSF and OTF for %s ...\n', figure_name)
    thr_eff_psz = 0.1 * 488 / 1.33 / 1000;
    xyPixelSize = 0.108;
    dz = 0.108;
    
    psf_mat = zeros([size(PSF_RW_510nm_NA1p0), numel(cme)]);
    for i = 1 : numel(cme)
        ind = cme(i);
        overall_PSF_3d = PSF_RW_510nm_NA1p0 .* permute(DitheredxzPSFCrossSection(ind, 501:1501), [1, 3, 2]);
        
        lambda = (ind - 1) * ystepsize_in;
    
        psf_mat(:, :, :, i) = overall_PSF_3d;
    
        % resize psfs to 0.108um
        overall_PSF_3d = imresize3(overall_PSF_3d, round([thr_eff_psz / xyPixelSize, thr_eff_psz / xyPixelSize, thr_eff_psz / dz] .* size(overall_PSF_3d)), 'linear');
        fnout = sprintf('%s/PSF_3d_px_108nm_%dlambda.tif', PSF_folder, lambda);
        writetiff(single(overall_PSF_3d), fnout);
    end
    clear PSF_RW_510nm_NA1p0 
    fprintf('3D theoretical PSF and OTF for %s is generated!\n', figure_name)
end


%% generate data

if generateData
    fprintf('Simulate stripe pattern images for %s ...\n', figure_name)
    %% simulate line image groundtruth    
    rng(10);
    im = zeros(799, 799); 
    
    n = 1 : 25;
    inds = (n - 1) .* n / 2 + 1;
    for i = 1 : 1
        if rem(i, 2) == 0
            inds_i = max(inds) + 1 - flip(inds) + (i - 1) * 500;
        else
            inds_i = inds + (i - 1) * 500;
        end
        inds_i(inds_i > size(im, 1)) = [];
        im(inds_i, :) = 1;
    end
    
    im([400, 650], :) = 1;
    
    im = padarray(im, [101, 101], 0, 'both');
    
    
    %% convolve with OTF
    
    psf_bbox = [1, 1, 1001, 1001];
    
    thr_eff_psz = 0.1 * 488 / 1.33 / 1000;
    xyPixelSize = 0.108;
    dz = 0.108;
    ysteps  = [1];
    
    num_ysteps = 3;
    
    % im raw 
    im_raw_mat = zeros(num_ysteps, 340, 340, 340);
    rng(1);

    % change to constant lines
    im_gd = im;

    im_gd_3d = padarray(permute(im_gd, [3, 2, 1]), [700, 200, 200], 0, 'both');

    % write out groundtruth before convolution (with crop)
    writetiff(single(im_gd_3d(201 : end - 200, 201 : end - 200, 201 : end - 200)), sprintf('%s/simulated_ground_truth.tif', resultPath));

    im_gd_3d_fft = fftn(im_gd_3d);
    
    for j = 1
        ind = cme(j);    
        lambda = (ind - 1) * ystepsize_in;
    
        psf = psf_mat(:, :, :, j);
        psf = single(psf);
        psf = psf ./ sum(psf(:));
    
        otf = psf2otf(psf, size(im_gd_3d_fft));
        
        im_raw = ifftn(im_gd_3d_fft .* otf);
        im_raw = im_raw(201 : end - 200, 201 : end - 200, 201 : end - 200);
        im_raw = im_raw .* (im_raw > 0);
        im_raw = imresize3(im_raw, round([thr_eff_psz / xyPixelSize, thr_eff_psz / xyPixelSize, thr_eff_psz / dz] .* size(im_raw)), 'linear');
        im_raw_gt = im_raw;
        writetiff(single(im_raw_gt), sprintf('%s/simulated_ground_truth_raw_%dlambda.tif', resultPath, lambda));
    end
    clear psf_mat;

    
    %% add noise for different snrs
    noise_sigma = 4;
    bg = 100;
    
    for i = 1
        ind = cme(i);    
        lambda = (ind - 1) * ystepsize_in;
        
        rawFn = sprintf('%s/simulated_ground_truth_raw_%dlambda.tif', resultPath, lambda);
        im_raw_gt = readtiff(rawFn);
        
        for j = 1 : numel(snrs)
            mint = snrs(j) .^ 2;        
            im_raw = im_raw_gt ./ prctile(im_raw_gt(:), 99.99) * mint;
            rng(1);
            im_gauss = bg + randn(size(im_raw)) * noise_sigma;
            im_poisson = randn(size(im_raw)) .* sqrt(im_raw);
            
            im_i = im_raw + im_gauss + im_poisson;
    
            fnout = sprintf('%s/simulated_raw_%dlambda_snr_%0.1f.tif', resultPath, lambda, snrs(j));
            writetiff(single(im_i), fnout)
        end
    end
    fprintf('Stripe pattern images for %s is generated!\n', figure_name)    
end


%% deconvolution 

fprintf('Deconvolve stripe pattern image with 3D theoretical PSF for %s and visualize results ...\n', figure_name)    

dataPath_exps = {resultPath};
for i = 1
    ind = cme(i);    
    lambda = (ind - 1) * ystepsize_in;

    disp(dataPath_exps);
    
    xyPixelSize = 0.108;
    dz = 0.108;
    Reverse = true;
    dzPSF = 0.108;
    parseSettingFile = ~true;
    
    ChannelPatterns = {
                        'raw_0lambda_snr', ...
                       };

    psfFn = sprintf('%s/PSF_3d_px_108nm_%dlambda.tif', PSF_folder, lambda);
    PSFFullpaths = {
                    psfFn, ...
                    }; 

    Background = 100;
    DeconIter = 80;
    Save16bit = true;
    fixIter = true;
    zarrFile = false;
    cpusPerTask = 24;
    largeFile = false;
    GPUJob = false;
    debug = true;
    psfGen = false;

    EdgeErosion = 0;

    deconPathstr = sprintf('matlab_decon_psf_%dlambda', lambda);

    XR_decon_data_wrapper(dataPath_exps, 'deconPathstr', deconPathstr, 'xyPixelSize', xyPixelSize, 'dz', dz, 'Reverse', Reverse, ...
        'ChannelPatterns', ChannelPatterns, 'PSFFullpaths', PSFFullpaths, 'dzPSF', dzPSF, ...
        'parseSettingFile', parseSettingFile, 'Background', Background, 'EdgeErosion', EdgeErosion, 'CPPdecon', false, ...
        'CudaDecon', ~true, 'DeconIter', DeconIter,  'fixIter', fixIter, 'zarrFile', zarrFile, 'Save16bit', Save16bit, ...
        'largeFile', largeFile, 'GPUJob', GPUJob, 'debug', debug, 'psfGen', psfGen, 'cpusPerTask', cpusPerTask, 'parseCluster', ~true);
end


%% FSC for raw data 

dataPath_exps = {resultPath};

disp(dataPath_exps);

xyPixelSize = 0.108;
dz = 0.108;
dr = 7.5;
dtheta = pi / 6;
resThreshMethod = 'one-bit';
resThresh = 0.2;
resAxis = 'xz';
% half size of region for FSC (typically leave some extra space for border)
N = [151, 151, 151];
bbox = [];
ChannelPatterns = {'simulated_raw_0lambda_snr'};
Channels = [488];
suffix = 'raw';
iterInterval = 5;
skipConeRegion = ~false;

XR_FSC_analysis_wrapper(resultPath, 'xyPixelSize', xyPixelSize, 'dz', dz, 'dr', dr, 'dtheta', dtheta, ...
    'resThreshMethod', resThreshMethod, 'resThresh', resThresh, 'resAxis', resAxis, 'N', N, ...
    'bbox', bbox, 'ChannelPatterns', ChannelPatterns, 'Channels', Channels, 'suffix', suffix, ...
    'iterInterval', iterInterval, 'skipConeRegion', skipConeRegion)


%% find the optimal number of iterations

dataPaths = arrayfun(@(x) sprintf('%s/matlab_decon_psf_%dlambda/', resultPath, x), (cme(1) - 1) * ystepsize_in, 'UniformOutput', 0);

dataPath_exps = cell(numel(dataPaths), 1);
for i = 1 : numel(dataPaths)
    dataPath_i = dataPaths{i};    
    dir_info = dir([dataPath_i, '*_debug']);
    dataPath_exps{i} = cellfun(@(x) [dataPath_i, x, '/'], {dir_info.name}', 'unif', 0);
end
dataPath_exps = cat(1, dataPath_exps{:});

disp(dataPath_exps);

xyPixelSize = 0.108;
dz = 0.108;
dr = 7.5;
dtheta = pi / 6;
resThreshMethod = 'one-bit';
resThresh = 0.2;
resAxis = 'xz';
% half size of region for FSC (typically leave some extra space for border)
N = [151, 151, 151];
bbox = [];
ChannelPatterns = {'lambda'};
Channels = [488];
suffix = 'decon';
iterInterval = 5;
skipConeRegion = ~false;

XR_FSC_analysis_wrapper(dataPath_exps, 'xyPixelSize', xyPixelSize, 'dz', dz, 'dr', dr, 'dtheta', dtheta, ...
    'resThreshMethod', resThreshMethod, 'resThresh', resThresh, 'resAxis', resAxis, 'N', N, ...
    'bbox', bbox, 'ChannelPatterns', ChannelPatterns, 'Channels', Channels, 'suffix', suffix, ...
    'iterInterval', iterInterval, 'skipConeRegion', skipConeRegion)


%% collect optimal number of iterations

iters = 5 : 5 : 80;
oiter_mat = zeros(numel(snrs), 1);
for i = 1
    lambda_psf = (cme(i) - 1) * ystepsize_in;    
    for j = 1
        for k = 1 : numel(snrs)
            lambda_j = (cme(j) - 1) * ystepsize_in;
            resFn = sprintf('%s/matlab_decon_psf_%dlambda/simulated_raw_%dlambda_snr_%0.1f_debug/FSCs/fsc_res_info_c1.mat', resultPath, lambda_psf, lambda_j, snrs(k));
            
            a = load(resFn);
            [~, mind] = min(a.res_mu_mat);
    
            oiter = iters(mind);
            oiter_mat(k) = oiter;
        end
    end
end

oiter = median(oiter_mat(:));

oiterFn = sprintf('%s/figures/fsc_optimal_iterations_snrs.mat', resultPath);
        
save('-v7.3', oiterFn, 'iters', 'oiter_mat', 'oiter', 'cme', 'snrs');

oiterFn = sprintf('%s/figures/fsc_optimal_iterations_snrs.txt', resultPath);
        
writematrix([snrs; oiter_mat'], oiterFn);


%% run deconvolution for different SNRs with optimal number of iterations

dataPath_exps = {resultPath};
for i = 1
    ind = cme(i);    
    lambda = (ind - 1) * ystepsize_in;
       
    disp(dataPath_exps);

    for k = 1 : numel(snrs)
        
        snr = snrs(k);

        xyPixelSize = 0.108;
        dz = 0.108;
        Reverse = true;
        dzPSF = 0.108;
        parseSettingFile = ~true;
        
        ChannelPatterns = {
                            sprintf('lambda_snr_%0.1f', snr), ...
                           };
                       
        psfFn = sprintf('%s/PSF_3d_px_108nm_%dlambda.tif', PSF_folder, lambda);
        PSFFullpaths = {
                        psfFn, ...
                        }; 
                    
                    
        Background = 100;
        DeconIter = oiter_mat(k);
        Save16bit = true;
        fixIter = true;
        zarrFile = false;
        cpusPerTask = 24;
        largeFile = false;
        GPUJob = false;
        debug = ~true;
        psfGen = false;
    
        EdgeErosion = 0;
    
        deconPathstr = sprintf('matlab_decon_optimal_iteration_psf_%dlambda', lambda);
    
        XR_decon_data_wrapper(dataPath_exps,  'deconPathstr', deconPathstr, 'xyPixelSize', xyPixelSize, 'dz', dz, 'Reverse', Reverse, ...
            'ChannelPatterns', ChannelPatterns, 'PSFFullpaths', PSFFullpaths, 'dzPSF', dzPSF, ...
            'parseSettingFile', parseSettingFile, 'Background', Background, 'EdgeErosion', EdgeErosion, 'CPPdecon', false, ...
            'CudaDecon', ~true, 'DeconIter', DeconIter,  'fixIter', fixIter, 'zarrFile', zarrFile, 'Save16bit', Save16bit, ...
            'largeFile', largeFile, 'GPUJob', GPUJob, 'debug', debug, 'psfGen', psfGen, 'cpusPerTask', cpusPerTask, 'parseCluster', ~true);
    end
end


%% save raw center slice (xz, yz), optimal center slice, and line cuts to mat file 

zc = 171;

line_bbox_y = [171, 171, 1, 171, 171, 340];

for i = 1
    lambda_psf = (cme(i) - 1) * ystepsize_in;    
    for j = 1
        lambda_j = (cme(j) - 1) * ystepsize_in;    

        % load optimal decon image
        deconPath = sprintf('%s/matlab_decon_psf_%dlambda/simulated_raw_%dlambda_snr_%0.1f_debug/', resultPath, lambda_psf, lambda_j, snr_ref);
        k = find(snrs == snr_ref);
        optmFn = sprintf('%s/Iter_%04d.tif', deconPath, oiter_mat(k));
        im_decon = readtiff(optmFn);
        im_decon_xz = squeeze(im_decon(zc, :, :));
        im_decon_yz = squeeze(im_decon(:, zc, :));
        
        bb = line_bbox_y;
        decon_linecut = squeeze(im_decon(bb(1) : bb(4), bb(2) : bb(5), bb(3) : bb(6)));
        
        % load raw decon image
        rawFn = sprintf('%s/simulated_raw_%dlambda_snr_%0.1f.tif', resultPath, lambda_j, snr_ref);
        im_raw = readtiff(rawFn);
        
        im_raw_xz = squeeze(im_raw(zc, :, :));
        im_raw_yz = squeeze(im_raw(:, zc, :));
        
        raw_linecut = squeeze(im_decon(bb(1) : bb(4), bb(2) : bb(5), bb(3) : bb(6)));
        
        % save central slices and linecuts
        fnout = sprintf('%s/raw_%dlambda_snr_%0.1f_decon_%dlambda_optimal_iteration_central_slices.mat', resultPath, lambda_j, snr_ref, lambda_psf);
        save('-v7.3', fnout, 'im_decon_xz', 'im_decon_yz', 'decon_linecut', 'im_raw_xz', 'im_raw_yz', 'raw_linecut')
    end
end


%% run deconvolution and save intermediate results every iteration

dataPath_exps = {resultPath};
for i = 1
    ind = cme(i);    
    lambda = (ind - 1) * ystepsize_in;

    disp(dataPath_exps);
    
    xyPixelSize = 0.108;
    dz = 0.108;
    Reverse = true;
    dzPSF = 0.108;
    parseSettingFile = ~true;
    
    ChannelPatterns = {
                        'raw_0lambda_snr', ...
                       };

    psfFn = sprintf('%s/PSF_3d_px_108nm_%dlambda.tif', PSF_folder, lambda);
    PSFFullpaths = {
                    psfFn, ...
                    }; 

    Background = 100;
    switch figure_name
        case {'Fig_Confined_AxialSW_NA0p45'}
            DeconIter = 100;
        otherwise
            DeconIter = 50;
    end

    Save16bit = true;
    fixIter = true;
    zarrFile = false;
    cpusPerTask = 24;
    largeFile = false;
    GPUJob = false;
    debug = true;
    saveStep = 1;
    psfGen = false;

    EdgeErosion = 0;

    deconPathstr = sprintf('matlab_decon_psf_%dlambda_save_every_iteration', lambda);

    XR_decon_data_wrapper(dataPath_exps,  'deconPathstr', deconPathstr, 'xyPixelSize', xyPixelSize, ...
        'dz', dz, 'Reverse', Reverse, 'ChannelPatterns', ChannelPatterns, 'PSFFullpaths', PSFFullpaths, ...
        'dzPSF', dzPSF, 'parseSettingFile', parseSettingFile, 'Background', Background, ...
        'EdgeErosion', EdgeErosion, 'CPPdecon', false, 'CudaDecon', ~true, 'DeconIter', DeconIter, ...
        'fixIter', fixIter, 'zarrFile', zarrFile, 'Save16bit', Save16bit, 'largeFile', largeFile, ...
        'GPUJob', GPUJob, 'debug', debug, 'saveStep', saveStep, 'psfGen', psfGen, 'cpusPerTask', cpusPerTask, ...
        'parseCluster', ~true);
end


%% run FSC every iteration

dataPaths = arrayfun(@(x) sprintf('%s/matlab_decon_psf_%dlambda_save_every_iteration/', resultPath, x), (cme(1) - 1) * ystepsize_in, 'UniformOutput', 0);

dataPath_exps = cell(numel(dataPaths), 1);
for i = 1 : numel(dataPaths)
    dataPath_i = dataPaths{i};    
    dir_info = dir([dataPath_i, '*_debug']);
    dataPath_exps{i} = cellfun(@(x) [dataPath_i, x, '/'], {dir_info.name}', 'unif', 0);
end
dataPath_exps = cat(1, dataPath_exps{:});

disp(dataPath_exps);

xyPixelSize = 0.108;
dz = 0.108;
dr = 7.5;
dtheta = pi / 6;
resThreshMethod = 'one-bit';
resThresh = 0.2;
resAxis = 'xz';
% half size of region for FSC (typically leave some extra space for border)
N = [151, 151, 151];
bbox = [];
ChannelPatterns = {'lambda'};
Channels = [488];
suffix = 'decon';
iterInterval = 1;
skipConeRegion = ~false;

XR_FSC_analysis_wrapper(dataPath_exps, 'xyPixelSize', xyPixelSize, 'dz', dz, 'dr', dr, 'dtheta', dtheta, ...
    'resThreshMethod', resThreshMethod, 'resThresh', resThresh, 'resAxis', resAxis, 'N', N, ...
    'bbox', bbox, 'ChannelPatterns', ChannelPatterns, 'Channels', Channels, 'suffix', suffix, ...
    'iterInterval', iterInterval, 'skipConeRegion', skipConeRegion)



%% make movies for decon across iterations

if ~false

lls_str = sprintf('%s, Annulus NA %1.2f/%1.2f, lattice NA=%1.2f, NAsigma=%1.2f', lightsheet_descrip, NAannulus(1), NAannulus(2), NAlattice, NAsigma);
fig_str = 'Simulated stripe pattern';

switch figure_name
    case 'FigS15_Confined_AxialSW'
        iters = 1 : 100;
    otherwise
        iters = 1 : 50;
end

opt_iter = oiter_mat(snr_ref == snrs);

frameRate = 24;
% required frame rate for each frame (variant frame rate)
propFR =  [linspace(2, 10, numel(iters) + 1), 1/2];
frameNums = round(frameRate ./ propFR);
% central slice
zc = 171;

cmap = jet(4096);
cmap(1 : 512, :) = [zeros(512, 2), linspace(0, 1, 512)'];
[ceR, ceB, ceG, ceP, ceO, cfR, cfB, cfG, cfP, cfO,cfK,ceK] = GU_colorCodes;

fontName = 'Helvetica';
fontSize = 14;

line_y = 1 : 340;

aps = [0.015, 0.075, 0.1, 0.85; 
       0.11, 0.075, 0.5, 0.85;
       0.61, 0.075, 0.125, 0.85; 
       0.78, 0.515, 0.2, 0.405;
       0.8, 0.05, 0.15, 0.4
       ];

aps = [0.035, 0.095, 0.1, 0.8; 
       0.1125, 0.095, 0.5, 0.8;
       0.6, 0.095, 0.125, 0.8; 
       0.775, 0.515, 0.195, 0.375;
       0.79, 0.05, 0.16, 0.42
       ];

% ground truth
gtFn = sprintf('%s/simulated_ground_truth.tif', resultPath);
gt = single(readtiff(gtFn));
% take central slice
gt = squeeze(gt(501, :, :))';
gt_line = gt(:, 501);
gt_line = gt_line ./ max(gt_line);
gt_line(gt_line < 0.001) = 0.001;


for i = 1
    lambda_psf = (cme(i) - 1) * ystepsize_in;    
    for j = 1
        lambda_j = (cme(j) - 1) * ystepsize_in;    
        deconPath = sprintf('%s/matlab_decon_psf_%dlambda_save_every_iteration/simulated_raw_%dlambda_snr_%0.1f_debug/', resultPath, lambda_psf, lambda_j, snr_ref);
        
        lambda_j = (cme(j) - 1) * ystepsize_in;

        % load FSC for decon data
        resFn = sprintf('%s/matlab_decon_psf_%dlambda_save_every_iteration/simulated_raw_%dlambda_snr_%0.1f_debug/FSCs/fsc_res_info_c1.mat', resultPath, lambda_psf, lambda_j, snr_ref);
        
        a = load(resFn);
        thetas = a.thetas;
        res_mat = a.res_mat;

        % replace infinity/nan with fillmissing
        res_mat(isinf(res_mat)) = nan;
        res_mat = fillmissing(res_mat, 'linear', 2);
        
        % load FSC for raw data
        resFn = sprintf('%s/FSCs/simulated_raw_%dlambda_snr_%0.1f.mat', resultPath, lambda_j, snr_ref);
        a = load(resFn);

        res_mat = [a.res_mu.res'; res_mat];
        res_mu_mat = mean(res_mat(:, 1 : end - 1), 2);
      
        vFn = sprintf('%s/simulated_raw_%dlambda_snr_%0.1f_decon_psf_%dlambda_with_linecut_FSC.avi', figurePath, lambda_j, snr_ref, lambda_psf);
        vPath = sprintf('%s/simulated_raw_%dlambda_snr_%0.1f_decon_psf_%dlambda_with_linecut_FSC', figurePath, lambda_j, snr_ref, lambda_psf);
        mkdir(vPath);
        v = VideoWriter(vFn, 'Uncompressed AVI');
        v.FrameRate = frameRate;
        open(v);
        
        for k = 1 : numel(iters) + 2                
            if k == 1
                fn = sprintf('%s/simulated_raw_%dlambda_snr_%0.1f.tif', resultPath, lambda_j, snr_ref);
            elseif k == numel(iters) + 2
                kopt = find(iters == opt_iter) + 1;
                fn = sprintf('%s/Iter_%04d.tif', deconPath, opt_iter); 
            else
                fn = sprintf('%s/Iter_%04d.tif', deconPath, iters(k-1));                
            end
            im = readtiff(fn);
            im_c = single(squeeze(im(zc, :, :))');
            im_cx = single(squeeze(im(:, zc, :))');

            figure('visible', 'off');
            
            set(gcf, 'Renderer', 'opengl', 'Position', [2100 10 1920 1080]);
            set(gcf, 'color', 'k')
            
            % yz central slice            
            h_1 = axes('position', aps(1, :));
                        
            imagesc(im_cx(:, 140 : 200))
            colormap(cmap)
            axis equal
            xlim([1, 61])
            ylim([1, 340])
            set(gca, 'YDir', 'normal')
            
            xticks([1, 57])
            xticklabels({'15', '21'})

            yticks([1, 168, 334])
            yticklabels({'0', '18', '36'})

            h_1.XAxis.FontSize = fontSize;
            h_1.XAxis.FontName = fontName;
            h_1.YAxis.FontSize = fontSize;
            h_1.YAxis.FontName = fontName;
            xlabel('y (\mum)', 'FontName', fontName, 'FontSize', fontSize)
            ylabel('z (\mum)', 'FontName', fontName, 'FontSize', fontSize)

            h_1.XColor = 'w';
            h_1.YColor = 'w';
            h_1.TickDir = 'none';


            hold on,                        
            % xz central slice
            h_2 = axes('position', aps(2, :));
            imagesc(im_c)
            colormap(cmap)            
            
            hold on
            plot(zc * ones(size(line_y)), line_y, 'color', cfO, 'linewidth', 1.5);
            axis equal
            xlim([1, 340])
            ylim([1, 340])
            set(h_2, 'YDir', 'normal')
            
            xticks([1, 167, 333])
            xticklabels({'0', '18', '36'})
            h_2.XAxis.FontSize = fontSize;
            h_2.XAxis.FontName = fontName;

            xlabel('x (\mum)', 'FontName', fontName, 'FontSize', fontSize)            
            set(h_2, 'ytick', [])
            h_2.XColor = 'w';
            h_2.YColor = 'w';
            h_2.TickDir = 'none';


            % line cut plot
            hold on
            h_3 = axes('position', aps(3, :));

            line_yf = squeeze(im_c(line_y, zc));
            if k == 1
                line_yf = line_yf - median(im_c(:));
                line_yf = line_yf .* (line_yf > 0);
            end
            line_yf = line_yf ./ max(line_yf);
            line_yf(line_yf < 0.001) = 0.001;
            plot(line_yf,  'color', cfO, 'linewidth', 1.5);
            xlim([1, line_y(end) - line_y(1) + 1]);
            ylim([0, 1]);

            % plot groundtruth line cut
            hold on
            plot(h_3, linspace(1, 340, numel(gt_line)), gt_line, 'color', cfB, 'linewidth', 1.0)

            set(h_3, 'xtick', [])
            set(h_3, 'view', [90, -90]);  
            h_3.YAxis.FontSize = fontSize;
            h_3.YAxis.FontName = fontName;
            set(h_3, 'yscale', 'log')
            ylim([0.01, 1])
            yticks([0.01, 0.1, 1])
            ylabel('Relative Intensity', 'FontName', fontName, 'FontSize', fontSize)

            h_3.Color = 'k';
            h_3.XColor = 'w';
            h_3.YColor = 'w';
            h_3.TickDir = 'none';

            % FSC plot
            hold on
            h_4 = axes('position', aps(4, :));
            if k < numel(iters) + 2
                plot([0, iters(1 : k - 1)]', res_mu_mat(1 : k), 'color', cfB, 'LineWidth', 1.5);
            else
                % only plot good smooth region
                [~, idx] = max(diff(res_mu_mat));
                plot([0, iters(1 : idx - 1)]', res_mu_mat(1 : idx), 'color', cfB, 'LineWidth', 1.5);
                hold on

                plot(opt_iter .* ones(1, 100), linspace(res_mu_mat(kopt) + 0.01, 1.0, 100)', '-', 'color', 'w', 'LineWidth', 1);
                
                % add tex for optimal iterations
                text(h_4, opt_iter - 5, 1.05, sprintf('Optimal iteration: %d', opt_iter), 'FontSize', fontSize, ...
                    'FontWeight', 'bold', 'FontName', fontName, 'color', 'w');
                ha = annotation('arrow', 'position', [opt_iter, res_mu_mat(kopt)+ 0.5, 0, -0.49], 'color', 'w');
                ha.Parent = h_4;
                
            end
            xlim([0, max(iters)]);
            ylim([0, 2]);
            xlabel('RL iteration', 'FontName', fontName, 'FontSize', fontSize)
            ylabel('Avg. xz resolution (a.u.)', 'FontName', fontName, 'FontSize', fontSize)

            hold on, 
            if k == 1
                plot(0, res_mu_mat(k), '.', 'color', ceB, 'MarkerSize', 30);
            elseif k == numel(iters) + 2                
            else
                plot(iters(k - 1), res_mu_mat(k), '.', 'color', ceB, 'MarkerSize', 30);
            end
            xlim([0, max(iters)]);
            ylim([0.25, 2]);
            yticks([0.25 : 0.5 : 2])
            grid on
            box off
            h_4.XAxis.FontSize = fontSize;
            h_4.XAxis.FontName = fontName;
            h_4.YAxis.FontSize = fontSize;
            h_4.YAxis.FontName = fontName;
            
            h_4.Color = 'k';
            h_4.XColor = 'w';
            h_4.YColor = 'w';
            h_4.TickDir = 'out';
            
            h_4.XAxis.LineWidth = 1;
            h_4.YAxis.LineWidth = 1;

            h_4.GridAlpha = 0.5;


            % FSC res polar plot
            hold on
            h_5 = polaraxes('position', aps(5, :));
            if k < numel(iters) + 2                        
                res_k = res_mat(k, :);
            else
                res_k = res_mat(kopt, :);
            end
            res_k(res_k > 2) = 2.4;
            p = polarplot(thetas, res_k, 'color', ceG, 'LineWidth', 1.5);
                
            rlim([0, 2.4])
            rticks([0 : 0.5 : 2, 2.4])
            rticklabels({'0', '0.5', '1', '1.5', '2', '2+'})
            h_5.FontSize = fontSize;
            h_5.FontName = fontName;
            h_5.RAxis.FontSize = 10;
           
            h_5.Color = 'k';
            h_5.RColor = 'w';
            h_5.ThetaColor = 'w';
            h_5.GridColor = 'w';
            h_5.GridAlpha = 0.5;

            % lattice description text
            text(h_2, 0.66, 1.07, fig_str, 'units', 'normalized', 'FontSize', 16, ...
                'FontWeight', 'bold', 'FontName', fontName, 'color', 'w');

            if k == 1
                title_str = sprintf('Raw Image');
            elseif k == numel(iters) + 2
                title_str = sprintf('Optimal deconvolution');                        
            else
                title_str = sprintf('Iter %02d', iters(k-1));
            end
            text(h_2, 0.835, 1.03, title_str, 'units', 'normalized', 'FontSize', 16, ...
                'FontWeight', 'bold', 'FontName', fontName, 'color', 'w', 'HorizontalAlignment', 'center');
            
            % also save high resolution png files
            figFn = sprintf('%s/frame%04d.png', vPath, k);
            set(gcf, 'InvertHardcopy', 'off');
            print(figFn, '-dpng', '-r0');
        
            cdata = print('-RGBImage', '-r0');
            fid = im2frame(cdata);
            for fidx = 1 : frameNums(k)
                writeVideo(v, fid);
            end
            close(gcf);
        end
        close(v);
    end 
end

end


%% psf analysis

% parameters for psf image
dz = 0.108;
angle = 32.45;
% parameters for psf analysis
ChannelPatterns = {'PSF_3d_px_108nm'};
% channel wavelength
Channels = [488];
Deskew = ~true;
ObjectiveScan = true;
bgfactor = 0;
save16bit = true;
% RW path for corresponding channels
RWFn = {
        [codeRT, '/RW_PSFs/PSF_RW_515em_128_128_101_100nmSteps.tif']
        };

% source description
sourceStr = 'PSFs';
masterCompute = ~true;

XR_psf_analysis_wrapper(PSF_folder, 'dz', dz, 'angle', angle, 'ChannelPatterns', ChannelPatterns, ...
    'Channels', Channels, 'Deskew', Deskew, 'ObjectiveScan', ObjectiveScan, 'sourceStr', sourceStr, ...
    'RWFn', RWFn, 'bgfactor', bgfactor, 'save16bit', save16bit, 'masterCompute', masterCompute);


%%

fprintf('Deconvlution and visualization is done!\n')

