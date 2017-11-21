%% Installing software of Sawada & Petrov (2017, doi:10.1152/jn.00821.2016) for Matlab

% (c) Laboratory for Cognitive Modeling and Computational Cognitive
% Neuroscience at the Ohio State University, http://cogmod.osu.edu
%
% Please read the LICENSE and NO WARRANTY statement in:
% SawadaPetrov_License.txt

mfilepath=fileparts(mfilename('fullpath'));
cd(mfilepath);

addpath(mfilepath);
addpath(fullfile(mfilepath,'Utility'));
addpath(fullfile(mfilepath,'Sim_EachExperiment'));
addpath(fullfile(mfilepath,'Sim_Figures'));
addpath(fullfile(mfilepath,'StandardDivisiveNormalization'));

movefile('Matlab\Sim_*')
