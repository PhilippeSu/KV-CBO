
addpath(pwd)

S = dir(fullfile(pwd,'*'));
N = setdiff({S([S.isdir]).name},{'.','..'}); 

for i = 1:numel(N)    
    addpath(genpath(N{i}))
end
