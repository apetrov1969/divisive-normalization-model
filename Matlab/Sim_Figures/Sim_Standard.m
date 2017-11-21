%% Simulation Experiments using the standardized Divisive Normalization model
%  The standard parameter set

% (c) Laboratory for Cognitive Modeling and Computational Cognitive
% Neuroscience at the Ohio State University, http://cogmod.osu.edu
%
% Please read the LICENSE and NO WARRANTY statement in:
% SawadaPetrov_License.txt


%% Clean up for debugging
cd(fileparts(mfilename('fullpath'))); % Change directory of this script
clear;
close all;


%% Generate a "model" structure M with default specifications
% Small stimulus size
model_level = 2;
M_Small = DMPL_defaultSpecs([],0,model_level);
M_Small = DMPL_prepareSpecs(M_Small);

% Large stimulus size
M_Large = M_Small;
M_Large.stim_spec.imageSize_pix = [128,128];
M_Large = DMPL_prepareSpecs(M_Large);

%%
complex = 1;
simple = @(ph) (mod(round(ph/90),4)+2);

neuron_orient_deg = 0;
neuron_spFreq_cpd = 2;

maxSize_Small = max(abs(M_Small.stim_spec.gridX_deg(:)))*2;
maxSize_Large = max(abs(M_Large.stim_spec.gridX_deg(:)))*2;


%% Subsection: Size tuning

    %% Size tuning
    specs_Size.contrast_pcent = 100;
    specs_Size.neuron_orient_deg = neuron_orient_deg;
    specs_Size.neuron_spFreq_cpd = neuron_spFreq_cpd;
    data_Size = Explore_Size(M_Large,specs_Size);

    data_Size.plot('complex');
    mCRF = data_Size.peak_x(complex);
    fprintf('measured CRF = %f\n',mCRF);

    %% Size tuning (contrast)
    specs_SizeByCon.neuron_orient_deg = neuron_orient_deg;
    specs_SizeByCon.neuron_spFreq_cpd = neuron_spFreq_cpd;
    data_SizeByCon = Explore_Size_by_Contrast(M_Large,specs_SizeByCon);

    data_SizeByCon.plot('complex');

    %% Size tuning (orientation)
        % Finding sub-optimal orientation
        subopt_specs_Orient.contrast_pcent = 100;
        subopt_specs_Orient.diameter_deg = mCRF;
        subopt_specs_Orient.neuron_orient_deg = neuron_orient_deg;
        subopt_specs_Orient.neuron_spFreq_cpd = neuron_spFreq_cpd;
        subopt_data_Orient = Explore_Orient(M_Small,subopt_specs_Orient);

    specs_SizeByOri.neuron_orient_deg = neuron_orient_deg;
    specs_SizeByOri.neuron_spFreq_cpd = neuron_spFreq_cpd;
    specs_SizeByOri.fwhh_orient_deg = subopt_data_Orient.bandwidth_deg(complex);
    data_SizByOri = Explore_Size_by_Orient(M_Large,specs_SizeByOri);

    data_SizByOri.plot('complex');

    %% Size tuning (spatial-frequency)
        % Finding sub-optimal spatial-frequency
        subopt_specs_SpFreq.contrast_pcent = 100;
        subopt_specs_SpFreq.diameter_deg = mCRF;
        subopt_specs_SpFreq.neuron_orient_deg = neuron_orient_deg;
        subopt_specs_SpFreq.neuron_spFreq_cpd = neuron_spFreq_cpd;
        subopt_data_SpFreq = Explore_SpFreq(M_Small,subopt_specs_SpFreq);

    specs_SizeBySpf.neuron_orient_deg = neuron_orient_deg;
    specs_SizeBySpf.neuron_spFreq_cpd = neuron_spFreq_cpd;
    specs_SizeBySpf.neuron_spFreq_idx = FindClosestIdx(M_Large.EarlyVis_spec.domain_spFreq_cpd,neuron_spFreq_cpd);
    specs_SizeBySpf.spFreq_low_oct  = subopt_data_SpFreq.marking_x(complex,2,specs_SizeBySpf.neuron_spFreq_idx);
    specs_SizeBySpf.spFreq_high_oct = subopt_data_SpFreq.marking_x(complex,3,specs_SizeBySpf.neuron_spFreq_idx);
    data_SizBySpf = Explore_Size_by_SpFreq(M_Large,specs_SizeBySpf);

    data_SizBySpf.plot('complex');

    %% Hole-size tuning (contrast)
    specs_HoleSizeByCon.neuron_orient_deg = neuron_orient_deg;
    specs_HoleSizeByCon.neuron_spFreq_cpd = neuron_spFreq_cpd;
    data_HoleSizeByCon = Explore_SizeHole_by_Contrast(M_Large,specs_HoleSizeByCon);

    data_HoleSizeByCon.plot('complex');


%% Subsection: Contrast sensitivity

    %%  Contrast sensitivity curve
    specs_Contrast.neuron_orient_deg = neuron_orient_deg;
    specs_Contrast.neuron_spFreq_cpd = neuron_spFreq_cpd;
    specs_Contrast.diameter_deg = maxSize_Small;
    data_Contrast = Explore_Contrast(M_Small,specs_Contrast);

    data_Contrast.plot('complex');

    %% Contrast sensitivity by noise
    specs_ConByNoise.neuron_orient_deg = neuron_orient_deg;
    specs_ConByNoise.neuron_spFreq_cpd = neuron_spFreq_cpd;
    specs_ConByNoise.diameter_deg = mCRF;
    data_ConByNoise = Explore_Contrast_by_Noise(M_Small,specs_ConByNoise);

    data_ConByNoise.plot('complex');

    %% Contrast sensitivity by orientation
    specs_ConByOri.neuron_orient_deg = neuron_orient_deg;
    specs_ConByOri.neuron_spFreq_cpd = neuron_spFreq_cpd;
    specs_ConByOri.diameter_deg = mCRF;
    data_ConByOri = Explore_Contrast_by_Orient(M_Small,specs_ConByOri);

    data_ConByOri.plot('complex');

    %% Contrast sensitivity by spatial-frequency
    specs_ConBySpf.neuron_orient_deg = neuron_orient_deg;
    specs_ConBySpf.neuron_spFreq_cpd = neuron_spFreq_cpd;
    specs_ConBySpf.diameter_deg = mCRF;
    data_ConBySpf = Explore_Contrast_by_SpFreq(M_Small,specs_ConBySpf);

    data_ConBySpf.plot('complex');

    %% Contrast sensitivity by size
    specs_ContBySize.neuron_orient_deg = neuron_orient_deg;
    specs_ContBySize.neuron_spFreq_cpd = neuron_spFreq_cpd;
    data_ContBySize = Explore_Contrast_by_Size(M_Large,specs_ContBySize);

    data_ContBySize.plot('complex');


%% Subsection: Orientation tuning

    %% Orientation tuning
    specs_Orient.neuron_orient_deg = neuron_orient_deg;
    specs_Orient.neuron_spFreq_cpd = neuron_spFreq_cpd;
    specs_Orient.contrast_pcent = 100;
    specs_Orient.diameter_deg = maxSize_Small;
    data_Orient = Explore_Orient(M_Small,specs_Orient);

    data_Orient.plot('complex');

    %% Orientation tuning by contrast
    specs_OriByCon.neuron_orient_deg = neuron_orient_deg;
    specs_OriByCon.neuron_spFreq_cpd = neuron_spFreq_cpd;
    specs_OriByCon.diameter_deg = maxSize_Small;
    data_OriByCon = Explore_Orient_by_Contrast(M_Small,specs_OriByCon);

    data_OriByCon.plot('complex');

    %% Orientation tuning by size
    specs_OriBySize.neuron_orient_deg = neuron_orient_deg;
    specs_OriBySize.neuron_spFreq_cpd = neuron_spFreq_cpd;
    specs_OriBySize.contrast_pcent = 100;
    specs_OriBySize.diameter_deg = mCRF;
    data_OriBySize = Explore_Orient_by_Size(M_Large,specs_OriBySize);

    data_OriBySize.plot('complex');


%% Subsection: Spatial-frequency tuning

    %% Spatial-frequency tuning
    specs_SpFreq.neuron_orient_deg = neuron_orient_deg;
    specs_SpFreq.neuron_spFreq_cpd = neuron_spFreq_cpd;
    specs_SpFreq.contrast_pcent = 100;
    specs_SpFreq.diameter_deg = maxSize_Small;
    data_SpFreq = Explore_SpFreq(M_Small,specs_SpFreq);

    data_SpFreq.plot('complex');

    %% Spatial-frequency tuning by contrast
    specs_SpfByCon.neuron_orient_deg = neuron_orient_deg;
    specs_SpfByCon.neuron_spFreq_cpd = neuron_spFreq_cpd;
    specs_SpfByCon.diameter_deg = maxSize_Small;
    data_SpfByCon = Explore_SpFreq_by_Contrast(M_Small,specs_SpfByCon);

    data_SpfByCon.plot('complex');

    %% Spatial-frequency tuning by size
    specs_SpfBySize.neuron_orient_deg = neuron_orient_deg;
    specs_SpfBySize.neuron_spFreq_cpd = neuron_spFreq_cpd;
    specs_SpfBySize.contrast_pcent = 100;
    specs_SpfBySize.diameter_deg = mCRF;
    data_SpfBySize = Explore_SpFreq_by_Size(M_Large,specs_SpfBySize);

    data_SpfBySize.plot('complex');

    
%% Subsection: X-orientation-suppression

    %% Orientation tuning of X-suppression
    specs_XOri.neuron_orient_deg = neuron_orient_deg;
    specs_XOri.neuron_spFreq_cpd = neuron_spFreq_cpd;
    specs_XOri.contrast1_pcent = 15;
    specs_XOri.contrast2_pcent = 25;
    specs_XOri.diameter_deg = maxSize_Small;
    data_XOri = Explore_XOrient(M_Small,specs_XOri);

    data_XOri.plot('complex');

    %% Spatial-frequency of X-suppression
    specs_XSpf.neuron_orient_deg = neuron_orient_deg;
    specs_XSpf.neuron_spFreq_cpd = neuron_spFreq_cpd;
    specs_XSpf.contrast1_pcent = 10;
    specs_XSpf.contrast2_pcent = 25;
    specs_XSpf.diameter_deg = maxSize_Small;
    data_XSpf = Explore_XSpFreq(M_Small,specs_XSpf);

    data_XSpf.plot('complex');

    %% Contrast sensitivity with X-suppression
    specs_XCont.neuron_orient_deg = neuron_orient_deg;
    specs_XCont.neuron_spFreq_cpd = neuron_spFreq_cpd;
    specs_XCont.diameter_deg = mCRF;
    data_XCont = Explore_XContrast_on_Contrast(M_Small,specs_XCont);

    data_XCont.plot('complex');


%% Subsection: Surround-suppression

    %% Orientation tuning with Surround-suppression
    specs_SurrOri.neuron_orient_deg = neuron_orient_deg;
    specs_SurrOri.neuron_spFreq_cpd = neuron_spFreq_cpd;
    specs_SurrOri.contrast_pcent = 100;
    specs_SurrOri.center_diameter_deg   = mCRF;
    specs_SurrOri.surround_outer_diameter_out_deg = maxSize_Large;
    specs_SurrOri.surround_inner_diameter_out_deg = mCRF;
    data_SurrOri = Explore_SurrOrient(M_Large,specs_SurrOri);

    data_SurrOri.plot('complex');

    %% Spatial-frequency with Surround-suppression
    specs_SurrSpf.neuron_orient_deg = neuron_orient_deg;
    specs_SurrSpf.neuron_spFreq_cpd = neuron_spFreq_cpd;
    specs_SurrSpf.contrast_pcent = 100;
    specs_SurrSpf.center_diameter_deg   = mCRF;
    specs_SurrSpf.surround_outer_diameter_out_deg = maxSize_Large;
    specs_SurrSpf.surround_inner_diameter_out_deg = mCRF;
    data_SurrSpf = Explore_SurrSpFreq(M_Large,specs_SurrSpf);
    
    data_SurrSpf.plot('complex');

    %% Contrast sensitivity with Surround-suppression
    specs_SurrCon.neuron_orient_deg = neuron_orient_deg;
    specs_SurrCon.neuron_spFreq_cpd = neuron_spFreq_cpd;
    specs_SurrCon.center_diameter_deg   = mCRF;
    specs_SurrCon.surround_outer_diameter_out_deg = maxSize_Large;
    specs_SurrCon.surround_inner_diameter_out_deg = mCRF;
    data_SurrCon = Explore_SurrContrast_on_Contrast(M_Large,specs_SurrCon);

    data_SurrCon.plot('complex');

    %% Contrast sensitivity with Surround-suppression(orientation)
    specs_SurrOri_Con.neuron_orient_deg = neuron_orient_deg;
    specs_SurrOri_Con.neuron_spFreq_cpd = neuron_spFreq_cpd;
    specs_SurrOri_Con.surr_contrast_pcent = 50;
    specs_SurrOri_Con.center_diameter_deg   = mCRF;
    specs_SurrOri_Con.surround_outer_diameter_out_deg = maxSize_Large;
    specs_SurrOri_Con.surround_inner_diameter_out_deg = mCRF;
    data_SurrOri_Con = Explore_SurrOrient_on_Contrast(M_Large,specs_SurrOri_Con);

    data_SurrOri_Con.plot('complex');

    
%% Subsection: Comparing Cross/Surround-suppression

    %% Orientation tuning of the Suppressive drive with Cross/Surround-suppression
    specs_SurrXOri.neuron_orient_deg = neuron_orient_deg;
    specs_SurrXOri.neuron_spFreq_cpd = neuron_spFreq_cpd;
    specs_SurrXOri.contrast_pcent = 100;
    specs_SurrXOri.center_diameter_deg   = mCRF;
    specs_SurrXOri.surround_outer_diameter_out_deg = maxSize_Large;
    specs_SurrXOri.surround_inner_diameter_out_deg = mCRF;
    data_SurrXOri = Explore_SurrXOrient_Supp(M_Large,specs_SurrXOri);

    data_SurrXOri.plot('complex');

    %% Spatial-frequency of the Suppressive drive with Cross/Surround-suppression
    specs_SurrXSpf.neuron_orient_deg = neuron_orient_deg;
    specs_SurrXSpf.neuron_spFreq_cpd = neuron_spFreq_cpd;
    specs_SurrXSpf.contrast_pcent = 100;
    specs_SurrXSpf.center_diameter_deg   = mCRF;
    specs_SurrXSpf.surround_outer_diameter_out_deg = maxSize_Large;
    specs_SurrXSpf.surround_inner_diameter_out_deg = mCRF;
    data_SurrXSpf = Explore_SurrXSpFreq_Supp(M_Large,specs_SurrXSpf);
    
    data_SurrXSpf.plot('complex');


%% Subsection: Measuring simple-cell receptive field

    %% Measuring Simple-cell with a light spot
    specs_Simple_LightSpot.neuron_orient_deg = neuron_orient_deg;
    specs_Simple_LightSpot.neuron_spFreq_cpd = neuron_spFreq_cpd;
    data_Simple_LightSpot = Explore_Simple_LightSpot(M_Small,specs_Simple_LightSpot);

    data_Simple_LightSpot.plot();

    %% Measuring Simple-cell with reverse-correlation method
    specs_Simple_RevCorr.neuron_orient_deg = neuron_orient_deg;
    specs_Simple_RevCorr.neuron_spFreq_cpd = neuron_spFreq_cpd;
    specs_Simple_RevCorr.size_block_pix = 2;
    specs_Simple_RevCorr.num_images = 50000;
    data_Simple_RevCorr = Explore_Simple_ReverseCorrelation(M_Small,specs_Simple_RevCorr);

    data_Simple_RevCorr.plot();

    %% Measuring Simple-cell with light/dark bars
    specs_Simple_Bar.neuron_orient_deg = neuron_orient_deg;
    specs_Simple_Bar.neuron_spFreq_cpd = neuron_spFreq_cpd;
    specs_Simple_Bar.bar_width = M_Small.stim_spec.degPerPixel;
    specs_Simple_Bar.contrast_pcent = 100;
    specs_Simple_Bar.diameter_deg = maxSize_Small;
    data_Simple_Bar = Explore_Simple_Bar(M_Small,specs_Simple_Bar);

    data_Simple_Bar.plotA('simple(0)');
    data_Simple_Bar.plotB('simple(0)');

      % Measured receptive-field by bar-width
      specs_Simple_Bar_Width.neuron_orient_deg = neuron_orient_deg;
      specs_Simple_Bar_Width.neuron_spFreq_cpd = neuron_spFreq_cpd;
      data_Simple_Bar_Width = Explore_Simple_Bar_Width(M_Small,specs_Simple_Bar_Width);

      data_Simple_Bar_Width.plot('simple(0)');

      % Measured receptive-field by bar-contrast
      specs_Simple_Bar_Contrast.neuron_orient_deg = neuron_orient_deg;
      specs_Simple_Bar_Contrast.neuron_spFreq_cpd = neuron_spFreq_cpd;
      data_Simple_Bar_Contrast = Explore_Simple_Bar_Contrast(M_Small,specs_Simple_Bar_Contrast);

      data_Simple_Bar_Contrast.plot('simple(0)');




