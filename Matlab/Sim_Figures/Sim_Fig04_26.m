%% Simulation Experiments using the Hyperbolic-ratio model

% (c) Laboratory for Cognitive Modeling and Computational Cognitive
% Neuroscience at the Ohio State University, http://cogmod.osu.edu
%
% Please read the LICENSE and NO WARRANTY statement in:
% SawadaPetrov_License.txt


%% Clean up for debugging
cd(fileparts(mfilename('fullpath'))); % Change directory of this script
clear;
close all;


%%
log10Contrast = 0:0.01:+2; % Figure 04
contrast = (10.^log10Contrast)/100; % Figure 04
% contrast = 0:0.001:1; % Figure 26

numContrast = size(contrast,2);

hbr = @(x, hbrM, hbrA, hbrN) ( hbrM*(x.^hbrN)./(hbrA^hbrN + x.^hbrN) );

tmpN = 2;
tmpA = 0.1;
tmpM = 1;

varA = [0.05,0.1,0.2]; %[0.01,0.02,0.05,0.1,0.2,0.5];
varN = [1,2,3]; %[1,1.5,2,2.5,3];


dataA = NaN(size(varA,2),numContrast);
dataN = NaN(size(varN,2),numContrast);

for i=1:size(varN,2)
    dataN(i,:) = hbr(contrast, tmpM,tmpA,varN(i));
end

for i=1:size(varA,2)
    dataA(i,:) = hbr(contrast, tmpM,varA(i),tmpN);
end


figure('Name','dataA');
%semilogx( contrast, dataA );
plot( contrast, dataA );


figure('Name','dataN');
%semilogx( contrast, dataN );
plot( contrast, dataN );