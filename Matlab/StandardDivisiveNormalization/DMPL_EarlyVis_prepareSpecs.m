function  M = DMPL_EarlyVis_prepareSpecs(M, levelIn, levelOut)


if (nargin<2); levelIn  = 0;   end
if (nargin<3); levelOut = Inf; end
assert(levelIn <= levelOut)


%% Check seed struct for consistency
assert(M.model_level >= 2)
assert(M.prepare_level >= 1)


%% Steps
%- Level 0: Basic preparation
M.prepare_level = max(M.prepare_level,2);
if (levelIn  <= 1); M.EarlyVis_spec = local_level_0(M); end
if (levelOut <= 0); return; end

%- Level 1: Image filters for deriving responses of Classical-Receptive-Field and Suppressive-Channels
if (levelIn <= 1)
    M.EarlyVis_filters = local_level_1(M);
    M.EarlyVis_spec.prepare_level = 1;
end
if (levelOut <= 1); return; end

%- Level 2: Kernels for xy-spatial pooling of the Suppressive-Channels
if (levelIn <= 2)
    M.EarlyVis_poolKernels = local_level_2(M);
    M.EarlyVis_spec.prepare_level = 2;
end
if (levelOut <= 2); return; end

%- Level 3: Kernels for orientation/frequency pooling of the Suppressive-Channels
if (levelIn <= 3)
    M.EarlyVis_divNorm = local_level_3(M);
    M.EarlyVis_spec.prepare_level = 3;
end
if (levelOut <= 3); return; end

%- Level 4: Calibration
if (levelIn <= 4)
    M.EarlyVis_calibration = local_level_4(M);
    M.EarlyVis_spec.prepare_level = 4;
end
if (levelOut <= 4); return; end

%- Level 5: Noise
% if (levelIn <= 5)
%     M.EarlyVis_calibration = local_level_5(M);
%     M.EarlyVis_spec.prepare_level = 5;
% end
% if (levelOut <= 5); return; end

%%% Return M
end  %%% of main function DMPL_EarlyVis_prepareSpecs


%%
%%
%% L O C A L     F U N C T I O N S
%%
%%


%% Level 0: Basic preparation
function  evS = local_level_0(M)

M.EarlyVis_spec.prepare_level = 0;

evS = M.EarlyVis_spec;

evS.num_rfLoc = size(M.EarlyVis_spec.domain_rfLoc_deg,1);

%- Add new fields to evS
evS.preparer = sprintf('%s, %s', mfilename(), datestr(now));

%- Check fields of evS
if(1~=strcmp(evS.phaseInvariance_method,{'quadrature','cos-only','sin-only'}))
    warning('Invalid M.EarlyVis_spec.phaseInvariance_method. ''quadrature'' assumed.');
    evS.phaseInvariance_method = 'quadrature';
end

if(1~=strcmp(evS.pool_xy_method,{'constant','linear'}))
    warning('Invalid M.EarlyVis_spec.pool_xy_method. ''constant'' assumed.');
    evS.pool_xy_method = 'constant';
end

if(1~=strcmp(evS.filters_type,{'gabor'}))
    warning('Invalid M.EarlyVis_spec.filters_type. ''gabor'' assumed.');
    evS.filters_type = 'gabor';
end

%%% Return evS
end  % of local function local_level_0


%%
%%


%% Level 1: Image filters for receptive-field convolution
function  EarlyVis_filters = local_level_1(M)

assert(M.EarlyVis_spec.prepare_level>=0);

stS = M.stim_spec;
evS = M.EarlyVis_spec;
num_chan_SpFreq = evS.num_channel_spFreq;
num_orient = evS.num_orient;
num_spFreq = evS.num_spFreq;
num_rfLoc  = evS.num_rfLoc;
num_pixels = prod(M.stim_spec.imageSize_pix);

%- Prepare to pad the stimulus with zeros
size_imgPadded  = pow2(floor(log2(stS.imageSize_pix))+1); % == 2^n where n is an arbitrary integer
size_leftBtmPad = floor((size_imgPadded - stS.imageSize_pix)./2);        % left  and bottom regions
size_rghtTopPad =  size_imgPadded - stS.imageSize_pix - size_leftBtmPad; % right and   top  regions

EarlyVis_filters.size_imgPadded  = size_imgPadded;
EarlyVis_filters.size_leftBtmPad = size_leftBtmPad;
EarlyVis_filters.size_rghtTopPad = size_rghtTopPad;

xaxis = 2; % in 'xy' image format, the x coordinate grows with the second index
yaxis = 1;
% iaxis = yaxis;
% jaxis = xaxis;

padLeft   = size_leftBtmPad(xaxis);
padBottom = size_leftBtmPad(yaxis);


%% Simple/Complex stimulus signal (classical receptive field)
%- Calculate the bandwidth parameters of the image filters
fwhh_clrf_bandwidth_spFreq_cpd = Log_to_Linear_FWHH(evS.fwhh_bandwidth_spFreq_oct, evS.domain_spFreq_cpd, 2); % [1,neurSpf]
EarlyVis_filters.fwhh_clrf_bandwidth_spFreq_cpd = fwhh_clrf_bandwidth_spFreq_cpd; % [1,neurSpf]

sigma_clrf_spFreq_cpd = NormalDistrib_FWHH_to_StdDev(fwhh_clrf_bandwidth_spFreq_cpd); % [1,neurSpf]
EarlyVis_filters.sigma_clrf_spFreq_cpd = sigma_clrf_spFreq_cpd; % [1,neurSpf]

sigma_clrf_orient_rad = NormalDistrib_FWHH_to_StdDev(evS.fwhh_bandwidth_orient_deg.*(pi/180)); % [1,1]
EarlyVis_filters.sigma_clrf_orient_rad = sigma_clrf_orient_rad; % [1,1]

%- Invert to the space domain (see Table 2.2 in Graham, 1989)
% The stdev width  determines the  frequency  bandwidth: (|) -> (||) -> (||...||)
% The stdev length determines the orientation bandwidth: (:) -> (==) -> (=.....=)
sigma_clrf_width  =(1/(2*pi)) ./  sigma_clrf_spFreq_cpd; % [1,neurSpf]
sigma_clrf_length =(1/(2*pi)) ./ (sigma_clrf_orient_rad.*evS.domain_spFreq_cpd); % [1,neurSpf]
EarlyVis_filters.sigma_clrf_width  = sigma_clrf_width;   % [1,neurSpf]
EarlyVis_filters.sigma_clrf_length = sigma_clrf_length;  % [1,neurSpf]

%- Allocate memory for the image filters
EarlyVis_filters.filters_clrf_CosImg = NaN(num_pixels,num_orient,num_spFreq,num_rfLoc); % [yx(ij),neurOri,neurSpf,neurRfl]
EarlyVis_filters.filters_clrf_SinImg = NaN(num_pixels,num_orient,num_spFreq,num_rfLoc); % [yx(ij),neurOri,neurSpf,neurRfl]

specs_gabor_clrf.gauss_peakHeight='normal'; % sum(gauss)==1
specs_gabor_clrf.phase_rad=0;
specs_gabor_clrf.amplitude=1;
for f = 1:num_spFreq
    specs_gabor_clrf.frequency    = evS.domain_spFreq_cpd(f);
    specs_gabor_clrf.sigma_width  = sigma_clrf_width(f);
    specs_gabor_clrf.sigma_length = sigma_clrf_length(f);
    
    for o = 1:num_orient
        specs_gabor_clrf.orientation_deg = evS.domain_orient_deg(o);
        
        for l = 1:num_rfLoc
            specs_gabor_clrf.center=evS.domain_rfLoc_deg(l,:);
            
            %-- Cos-phase
            specs_gabor_clrf.type='cos';
            gaborCos = Gabor2D(specs_gabor_clrf,stS.gridX_deg,stS.gridY_deg);
            EarlyVis_filters.filters_clrf_CosImg(:,o,f,l) = gaborCos(:);
            
            %-- Sin-phase
            specs_gabor_clrf.type='sin';
            gaborSin = Gabor2D(specs_gabor_clrf,stS.gridX_deg,stS.gridY_deg);
            EarlyVis_filters.filters_clrf_SinImg(:,o,f,l) = gaborSin(:);
            
        end % for l = 1:num_rfLoc
    end  % for o = 1:num_orient
end  % for f = 1:num_spFreq


%% Surround suppressive signal
%- Calculate the bandwidth parameters of the image filters
fwhh_chan_bandwidth_spFreq_cpd = Log_to_Linear_FWHH(evS.fwhh_bandwidth_spFreq_oct, evS.domain_channel_spFreq_cpd, 2); % [1,chanSpf]
EarlyVis_filters.fwhh_chan_bandwidth_spFreq_cpd = fwhh_chan_bandwidth_spFreq_cpd; % [1,chanSpf]

sigma_chan_spFreq_cpd = NormalDistrib_FWHH_to_StdDev(fwhh_chan_bandwidth_spFreq_cpd); % [1,chanSpf]
EarlyVis_filters.sigma_chan_spFreq_cpd = sigma_chan_spFreq_cpd; % [1,chanSpf]

sigma_chan_orient_rad = NormalDistrib_FWHH_to_StdDev(evS.fwhh_bandwidth_orient_deg.*(pi/180)); % [1,1]
EarlyVis_filters.sigma_chan_orient_rad = sigma_chan_orient_rad; % [1,1]

%- Invert to the space domain (see Table 2.2 in Graham, 1989)
% The stdev width  determines the  frequency  bandwidth: (|) -> (||) -> (||...||)
% The stdev length determines the orientation bandwidth: (:) -> (==) -> (=.....=)
sigma_chan_width  =(1/(2*pi)) ./  sigma_chan_spFreq_cpd; % [1,chanSpf]
sigma_chan_length =(1/(2*pi)) ./ (sigma_chan_orient_rad.*evS.domain_channel_spFreq_cpd); % [1,chanSpf]
EarlyVis_filters.sigma_chan_width  = sigma_chan_width;   % [1,chanSpf]
EarlyVis_filters.sigma_chan_length = sigma_chan_length;  % [1,chanSpf]

%- Allocate memory for the image filters
EarlyVis_filters.filters_chan_CosImg = cell(num_orient,num_chan_SpFreq); % {chanOri,chanSpf}
EarlyVis_filters.filters_chan_SinImg = cell(num_orient,num_chan_SpFreq); % {chanOri,chanSpf}

%- Allocate memory for the Fourier spectrums of the image filters
EarlyVis_filters.filters_chan_CosFFT = cell(num_orient,num_chan_SpFreq); % {chanOri,chanSpf}
EarlyVis_filters.filters_chan_SinFFT = cell(num_orient,num_chan_SpFreq); % {chanOri,chanSpf}

%- Compute the filters
padded_rf = zeros(size_imgPadded); % background of receptive field is always 0
specs_gabor_chan.gauss_peakHeight='normal'; % sum(gauss)==1
specs_gabor_chan.phase_rad=0;
specs_gabor_chan.amplitude=1;
for f = 1:num_chan_SpFreq
    specs_gabor_chan.frequency    = evS.domain_channel_spFreq_cpd(f);
    specs_gabor_chan.sigma_width  = sigma_chan_width(f);
    specs_gabor_chan.sigma_length = sigma_chan_length(f);
    
    for o = 1:num_orient
        specs_gabor_chan.orientation_deg = evS.domain_orient_deg(o);
        
        specs_gabor_chan.center=[0,0];
        specs_gabor_chan.type='cos';
        gaborCos = Gabor2D(specs_gabor_chan,stS.gridX_deg,stS.gridY_deg);
        EarlyVis_filters.filters_Surr_CosImg{o,f} = gaborCos;
        specs_gabor_chan.type='sin';
        gaborSin = Gabor2D(specs_gabor_chan,stS.gridX_deg,stS.gridY_deg);
        EarlyVis_filters.filters_Surr_SinImg{o,f} = gaborSin;

        %-- Fourier spectrum of filters are used in Convolution theorem
        padded_rf((1:stS.imageSize_pix(yaxis))+padBottom,...
                  (1:stS.imageSize_pix(xaxis))+padLeft) = gaborCos;
        EarlyVis_filters.filters_Surr_CosFFT{o,f} = conj( fft2(fftshift(padded_rf)) );

        padded_rf((1:stS.imageSize_pix(yaxis))+padBottom,...
                  (1:stS.imageSize_pix(xaxis))+padLeft) = gaborSin;
        EarlyVis_filters.filters_Surr_SinFFT{o,f} = conj( fft2(fftshift(padded_rf)) );
        
    end  % for o = 1:num_orient
end  % for f = 1:num_chan_SpFreq

%%% Return EarlyVis_filters
end  % of local function local_level_1


%%
%%


%% Level 2: Kernels for xy-spatial pooling of the suppressive-drive
function  EarlyVis_poolKernels = local_level_2(M)

assert(M.EarlyVis_spec.prepare_level>=1)

stS = M.stim_spec;
evS = M.EarlyVis_spec;

%- Convert the kernel bandwidth parameters (sigma = standard deviation)
sigmaSurr_deg = NaN(evS.num_spFreq,1); % [neurSpf,1]
switch (evS.pool_xy_method)
    case 'constant'
        if isfield(evS,'fwhh_pool_surr_xy_deg')
            sigmaSurr_deg(:) = NormalDistrib_FWHH_to_StdDev(evS.fwhh_pool_surr_xy_deg);
        elseif isfield(evS,'fwhh_pool_surr_xy_cyc')
            sigmaSurr_deg(:) = NormalDistrib_FWHH_to_StdDev(evS.fwhh_pool_surr_xy_cyc/M.EarlyVis_calibration.calibr_spFreq_cpd);
        end
    case 'linear'
        if isfield(evS,'fwhh_pool_surr_xy_cyc')
            sigmaSurr_deg(:) = NormalDistrib_FWHH_to_StdDev(evS.fwhh_pool_surr_xy_cyc./evS.domain_spFreq_cpd(:));
        elseif isfield(evS,'fwhh_pool_surr_xy_deg')
            sigmaSurr_deg(:) = NormalDistrib_FWHH_to_StdDev(evS.fwhh_pool_surr_xy_deg*M.EarlyVis_calibration.calibr_spFreq_cpd./evS.domain_spFreq_cpd(:));
        end
    otherwise
        error('Invalid M.EarlyVis_spec.pool_xy_method.');
end
EarlyVis_poolKernels.sigmaSurr_deg = sigmaSurr_deg; % [neurSpf,1]

%- Gaussian kernels for spatial pooling
EarlyVis_poolKernels.filters_pool_surr_yx = NaN(evS.num_rfLoc,evS.num_spFreq,stS.num_pixels); % [neurRfl,neurSpf,yx]

xaxis = 2;  % in 'xy' image format, the x coordinate grows with the second index
yaxis = 1;

specs_gauss_surr.peakHeight = 'normal'; % sum(gauss)==1
for sf = 1:evS.num_spFreq
    specs_gauss_surr.sigmaX = sigmaSurr_deg(sf);
    specs_gauss_surr.sigmaY = sigmaSurr_deg(sf);
    for rf = 1:evS.num_rfLoc
        specs_gauss_surr.centerX = evS.domain_rfLoc_deg(rf,xaxis);
        specs_gauss_surr.centerY = evS.domain_rfLoc_deg(rf,yaxis);
        filterGaussSurr = Gauss2D(specs_gauss_surr,stS.gridX_deg,stS.gridY_deg);
        EarlyVis_poolKernels.filters_pool_surr_yx(rf,sf,:) = filterGaussSurr(:); % [neurRfl,neurSpf,yx]
    end % rf = 1:evS.num_rfLoc
end % sf = 1:evS.num_spFreq

%%% Return EarlyVis_poolKernels
end  % of local function local_level_2


%%
%%


%% Level 3: Kernels for orientation/frequency pooling of the suppressive-drive
function  EarlyVis_divNorm = local_level_3(M)

assert(M.EarlyVis_spec.prepare_level>=2);

evS = M.EarlyVis_spec;

%% Kernel for orientation pooling: von Mises distribution
%- Convert the orientation pooling parameters
%   von Mises distribution use the direction domain that has a period of 360 degrees
%   whereas the orientation domain has a period of 180 degrees
temp_fwmh_pool_surr_direction_deg = 2 * evS.fwmh_pool_surr_orient_deg; % [0~180] -> [0~+360]
temp_domain_orient_rad = (2*pi/180).*evS.domain_orient_deg;            % [0~180] -> [0~2pi]

%- Convert full width at middle height (FWMH) to Kappa
orient_kernel_kappa = VonMisesDistrib_FWMH_to_Kappa(temp_fwmh_pool_surr_direction_deg);
EarlyVis_divNorm.orient_kernel_kappa = orient_kernel_kappa;

%- Make a matrix representing the orientation kernel
orient_kernel_weightMatrix = NaN(evS.num_orient,evS.num_orient); % [neurOri,chanOri]
if (evS.fwmh_pool_surr_orient_deg <= 3) % von Mises distribution -> delta function
    orient_kernel_weightMatrix(:) = 0;
    for o = 1:evS.num_orient
        orient_kernel_weightMatrix(o,o) = 1; % delta function
    end
else
    for o = 1:evS.num_orient
        temp_vector = temp_domain_orient_rad-temp_domain_orient_rad(o);
        temp_vector = exp(orient_kernel_kappa.*cos(temp_vector)); % von Mises
        temp_vector = temp_vector/sum(temp_vector); % sum(temp_vector)==1
        orient_kernel_weightMatrix(o,:) = temp_vector; % [neurOri,chanOri]
    end % o = 1:evS.num_orient
end
EarlyVis_divNorm.orient_kernel_weightMatrix = orient_kernel_weightMatrix; % [neurOri,chanOri]


%% Kernel for frequency pooling: Gaussian distribution
%- Convert full width at half height (FWHH) to Standard deviation (sigma)
spFreq_kernel_sigma_oct = NormalDistrib_FWHH_to_StdDev(evS.fwhh_pool_surr_spFreq_oct); % [1,1]
EarlyVis_divNorm.spFreq_kernel_sigma_oct = spFreq_kernel_sigma_oct; % [1,1]

%- Make a matrix representing the frequency kernel
spFreq_kernel_weightMatrix = NaN(evS.num_spFreq,evS.num_channel_spFreq); % [neurSpf,chanSpf]
for f = 1:evS.num_spFreq
    temp_vector = evS.domain_channel_spFreq_l2cpd - evS.domain_spFreq_l2cpd(f) - evS.shift_pool_surr_spFreq_oct;
    temp_vector = exp( -(temp_vector.^2)./(2*spFreq_kernel_sigma_oct.^2) ); % Gauss
    temp_vector = temp_vector/sum(temp_vector); % sum(temp_vector)==1
    spFreq_kernel_weightMatrix(f,:) = temp_vector;
end % f = 1:evS.num_spFreq
EarlyVis_divNorm.spFreq_kernel_weightMatrix = spFreq_kernel_weightMatrix; % [neurSpf,chanSpf]

%%% Return EarlyVis_divNorm
end  % of local function local_level_3


%%
%%


%% Level 4: Calibration
function  EarlyVis_calibration = local_level_4(M)

assert(M.EarlyVis_spec.prepare_level>=3);

EarlyVis_calibration = M.EarlyVis_calibration;
evS = M.EarlyVis_spec;
stS = M.stim_spec;

%% Make a 2D grating image for the calibration
calibr_orient_deg = EarlyVis_calibration.calibr_orient_deg;
calibr_spFreq_cpd = EarlyVis_calibration.calibr_spFreq_cpd;

calibrationGrating_stim = NaN([stS.imageSize_pix]);
grating_spec_calibr.amplitude = stS.rangeLum/2;
grating_spec_calibr.frequency = calibr_spFreq_cpd;
grating_spec_calibr.orientation_deg = calibr_orient_deg;
grating_spec_calibr.centerX   = 0;
grating_spec_calibr.centerY   = 0;
grating_spec_calibr.phase_deg = 0;
grating_spec_calibr.type = 'Cos';
calibrationGrating_stim(:,:) = Grating2D(grating_spec_calibr,stS.gridX_deg,stS.gridY_deg);
calibrationGrating_stim = calibrationGrating_stim + stS.bgLum;
EarlyVis_calibration.calibrGrating_stim = calibrationGrating_stim;

calibrationGrating_supp = NaN([stS.imageSize_pix]);
grating_spec_calibr.frequency = calibr_spFreq_cpd + evS.shift_pool_surr_spFreq_oct;
calibrationGrating_supp(:,:) = Grating2D(grating_spec_calibr,stS.gridX_deg,stS.gridY_deg);
calibrationGrating_supp = calibrationGrating_supp + stS.bgLum;
EarlyVis_calibration.calibrGrating_supp = calibrationGrating_supp;

%% Derive the calibration factors
EarlyVis_calibration.calibrFactor_numerator_simp = 1;
EarlyVis_calibration.calibrFactor_numerator_comp = 1;
EarlyVis_calibration.calibrFactor_denominator    = 1;

[calibr1_Comp,...
 calibr1_Simp,...
 calibr1_Supp] = DMPL_EarlyVis_PoolRectifiedLinChannels(M,calibrationGrating_stim); 

%- The 1st calibration for stimulus drive
values_simp = calibr1_Simp(... % [neurOri, neurSpf, neurRfl, neurPhs, img]
                          EarlyVis_calibration.calibr_orient_cell_idx, ... % neurOri
                          EarlyVis_calibration.calibr_spFreq_cell_idx, ... % neurSpf
                          EarlyVis_calibration.calibr_rfLoc_cell_idx,  ... % neurRfl
                          :,:);                                            % neurPhs,img
values_comp = calibr1_Comp(... % [neurOri, neurSpf, neurRfl, img]
                          EarlyVis_calibration.calibr_orient_cell_idx, ... % neurOri
                          EarlyVis_calibration.calibr_spFreq_cell_idx, ... % neurSpf
                          EarlyVis_calibration.calibr_rfLoc_cell_idx,  ... % neurRfl
                          :);                                              % img
calibrFactor_simp = max(values_simp(:));
assert(calibrFactor_simp > 0);
calibrFactor_comp = max(values_comp(:));
assert(calibrFactor_comp > 0);

EarlyVis_calibration.calibrFactor_numerator_simp = calibrFactor_simp;
EarlyVis_calibration.calibrFactor_numerator_comp = calibrFactor_comp;
M.EarlyVis_calibration = EarlyVis_calibration;

%- The 2nd calibration for suppressive drive
[calibr2_Comp,...
 calibr2_Simp,...
 calibr2_Supp] = DMPL_EarlyVis_PoolRectifiedLinChannels(M,calibrationGrating_supp); 

[dummy1,calibr_suppressDrive,dummy2] = DMPL_EarlyVis_MeanChannelResp(M, calibr2_Comp, ...
                                                                        calibr2_Simp, ...
                                                                        calibr2_Supp);
                                           
calibrFactor_denominator = calibr_suppressDrive(EarlyVis_calibration.calibr_orient_cell_idx, ... % neurOri
                                                EarlyVis_calibration.calibr_spFreq_cell_idx, ... % neurSpf
                                                EarlyVis_calibration.calibr_rfLoc_cell_idx,  ... % neurRfl
                                                :);
calibrFactor_denominator = max(calibrFactor_denominator(:));
assert(calibrFactor_denominator > 0);

EarlyVis_calibration.calibrFactor_denominator = calibrFactor_denominator;

%%% Return EarlyVis_calibration
end  % of local function local_level_4


%%
%%


%% Level 5: Noise
% function  EarlyVis_calibration = local_level_5(M)
% 
% %%% Return EarlyVis_calibration
% end  % of local function local_level_5

%%%%%%%% End of file