function [output] = MOT(action, object, object_feature)
% PERFORMS BASIC LOAD/SAVE FUNCTIONS REFERENCE TO MOT
% action = 'save' or 'load' or 'savemovie' or 'loadmovie'

% Input and output will depend on the case ('action')
% Use according to 'action':
%   [~]= MOT('save', VSDI) or  MOT('save', VSDI)-uses internal VSDI.ref to save in appropiate
%	.mat
%   [VSDI] = MOT('load', nfish)
%  [~]= MOT('savemovie', VSDmov, movierefernce) - uses moviereference
%  (~char) to name the matfile

%  [VSDmov]= MOT('loadmovie', nfish, moviereference) - uses moviereference
%  [VSDroiTS]= MOT('loadmovie', nfish, moviereference) - uses moviereference
%  [spike]= MOT('loadspike', nfish, moviereference) - uses moviereference

% nsubject = MOT('nsubject',    ref)


datapath = 'C:\Users\User\Documents\MatLab\MOT\data';

VSDIpath = fullfile(datapath,'dataVSDI');
moviepath = fullfile(datapath,'datamovies');
wavespath = fullfile(datapath,'datawaves');
spikepath = fullfile(datapath,'dataspike');

expref = 'MOT';
nchar = length(expref);



% Input control
switch action
    case 'save'
        if  ~isstruct(object)
            disp('the input is not what expected'); end
        if  ~strcmpi(object.expref, expref)
            warning('It cannot be saved: the structure"s experiment does not match the function')
            return
        end
    case 'load'
        %         assert(mod(object, 1) == 0 && , 'input to load must be a single number');
        
        try
            load(fullfile(datapath, 'grouplist.mat'))
        catch
            warning('fish cannot be load because "grouplist.mat" does not exist')
        end
        
        
    case 'savemovie'
        %             if ~exist('object_feature')
        %                 error('input a proper reference name for the movie (as 3rd argument)'); end
        object_feature = object.movieref;
        if  ~strcmpi(object.expref, expref)
            warning('It cannot be saved: the structure"s experiment does not match the function')
            return
        end
        
end % input control

%% FUNCTION CODE:

switch action
    case 'save'
        VSDI = object;
        %saveVSDI saves current VSDI structure respect to the current rootpath
        pathname = fullfile(VSDIpath,[expref '_' num2str(object.ref) '.mat']);
        save(pathname, 'VSDI')
        
    case 'load'
        load(fullfile(datapath, 'grouplist')) %load structure list to take the fish reference
        if object > length(grouplist)
            warning('there are not so many fish in grouplist')
            return
        end
        
        load(fullfile(VSDIpath,[grouplist{object},'.mat'])) %load the fish VSDI
        if object > length(grouplist)
            warning('there are not so many fish in grouplist')
            return
        end
        
        disp(strcat (grouplist{object}, '_loaded'));
        output= VSDI;
        
    case 'savemovie'
        VSDmov= object;
        %saveVSDI saves current VSDI structure respect to the current rootpath
        pathname = fullfile(moviepath,['MOTMov_',num2str(VSDmov.ref),object_feature,'.mat']);
        save(pathname,'VSDmov','-v7.3')
        
    case 'loadmovie'
        load(fullfile(datapath, 'grouplist'))
        if object > length(grouplist)
            warning('there are not so many fish in grouplist')
            return
        end
        
        fishref = grouplist{object}(nchar+2:end);
        if object > length(grouplist)
            warning('there are not so many fish in grouplist')
            return
        end
        
        %saveVSDI saves current VSDI structure respect to the current rootpath
        movieref = [expref,'Mov_',fishref,object_feature,'.mat'];
        load(fullfile(moviepath,movieref))
        output= VSDmov;
        disp([movieref, '_loaded']);
        
    case 'savewave'
        VSDroiTS = object;
        %saveVSDI saves current VSDI structure respect to the current rootpath
        pathname = fullfile(wavespath,[expref 'RoiTS_' num2str(object.ref) '.mat']);
        save(pathname, 'VSDroiTS')
        
    case 'loadwave'
        load(fullfile(datapath, 'grouplist')) %load structure list to take the fish reference
        if object > length(grouplist)
            warning('there are not so many fish in grouplist')
            return
        end
        
        fishref = grouplist{object}(nchar+2:end);
        load(fullfile(wavespath,[expref 'RoiTS_',fishref,'.mat'])) %load the fish VSDI
        disp(strcat ('ROIs timeseries for fish',grouplist{object}, '_loaded'));
        output= VSDroiTS;
        
    case 'savespike'
        spike = object;
        %saveVSDI saves current VSDI structure respect to the current rootpath
        pathname = fullfile(spikepath,[expref '_spike' num2str(object.ref) '.mat']);
        save(pathname, 'spike')
        
    case 'loadspike'
        load(fullfile(datapath, 'grouplist')) %load structure list to take the fish reference
        fishref = grouplist{object}(nchar+2:end);
        load(fullfile(spikepath,[expref '_spike',fishref,'.mat'])) %load the fish VSDI
        disp(strcat ('Spike structure for fish',grouplist{object}, '_loaded'));
        output= spike;
        
    case 'who'
        load(fullfile(datapath, 'grouplist')) %load structure list to take the fish reference
        nref = object; if  isnumeric(nref); nref = num2str(nref); end
        name = [expref, '_', num2str(nref)];
        output = find(strcmpi([grouplist(:)],name)); %number of fish
        
end %switch
end

% function T = isIntegerValue(X)
% T = (mod(X, 1) == 0);
% end

%% Created: 31/01/2021
% Updated: 08/02/21
% 21/02/2022 - add safety functionality: save only if the structure.expref
% matches the function's
