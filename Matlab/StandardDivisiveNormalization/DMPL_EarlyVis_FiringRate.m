function  retval_resp_divNorm = DMPL_EarlyVis_FiringRate(M,stimuli)

%DMPL_EarlyVis_FiringRate - Compute response (firing rate) of the DNM
%
% See also DMPL_EarlyVis_PoolRectifiedLinChannels, DMPL_EarlyVis_MeanChannelResp

% (c) Laboratory for Cognitive Modeling and Computational Cognitive
% Neuroscience at the Ohio State University, http://cogmod.osu.edu
%
% 1.0.0 2017-07-25 TS -- Wrote it for simplicity


%% Calculate outputs of linear filters 
[resp_CentComp,...
 resp_CentSimp,...
 resp_SuppChan] = DMPL_EarlyVis_PoolRectifiedLinChannels(M, stimuli);


%% Operate divisive normalization process
[resp_DivNorm,...
 resp_SuppDrv,...
 resp_StimDrv]  = DMPL_EarlyVis_MeanChannelResp(M,resp_CentComp,...
                                                  resp_CentSimp,...
                                                  resp_SuppChan); % [neOri,neSpf,neRfl,neTyp,img]
                                              
retval_resp_divNorm = resp_DivNorm; % [neOri,neSpf,neRfl,neTyp,img]

%% Return retval_resp_divNorm
end  %%% of file
