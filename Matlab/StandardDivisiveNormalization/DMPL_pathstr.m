function  pathstr = DMPL_pathstr

%DMPL_pathstr - Pathstring to the root of the current Dimple implementation
%
% pathstr = DMPL_pathstr
%
% Example:  DMPL_pathstr --> '/Users/apetrov/work/models/Dimple0'
%
% See also cd, fileparts, fullfile, pwd, work_pathstr.

% (c) Laboratory for Cognitive Modeling and Computational Cognitive
% Neuroscience at the Ohio State University, http://cogmod.osu.edu
%
% 1.0.0 2013-06-15 AAP -- Wrote it based on wort_pathstr.m


%% Call standard function 'which' to find the Dimple version that is 
% currently on the Matlab path.
%pathstr = fileparts(which('DMPL_pathstr')) ;  % '~/work/models/Dimple0/Dimple0sa'
pathstr = fileparts(which(mfilename())) ;      % '~/work/models/Dimple0/Dimple0sa'


%% For the time being, we have minor versions (e.g., Dimple0sa) in subdirectories.
%  This necessitates going one level up the directory tree.
idx = find(pathstr==filesep) ;  % indices of all separator characters
cutoff = idx(end) ;             % 
pathstr(cutoff:end) = [] ;                    % '~/work/models/Dimple0'


%% Return pathstr
end  %%% of file
