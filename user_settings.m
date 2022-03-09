% If 'pwd' is used instead of manually adding the rootpath of the experiment folder, user_settings
% has to be called from that rootpath folder (the one with your
% experiment name) everytime you call it from a script.

% Option 1 (manual)
path.rootpath = 'C:\Users\User\Documents\MatLab\MOT';

% Option 2 (pwd) - to be execute in the experiment folder

%  path.rootpath = pwd; 
 temp.idcs   = strfind(path.rootpath,'\');
 temp.newdir = path.rootpath(1:temp.idcs(end));

path.data = fullfile(path.rootpath, 'data');
path.grouplist = fullfile(path.rootpath);
path.list =fullfile(path.rootpath, 'data', 'BVlists');


addpath(genpath(fullfile(temp.newdir, 'VSDI_ourToolbox', 'functions')));
addpath(genpath(path.rootpath));


clear temp
        