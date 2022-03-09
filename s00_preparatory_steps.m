%% PREPARATORY STEPS

% ROOT FOLDER content: 
% pipeline files, 'savecode'
% 'functions' folder
% 'data_structures' folder - for VSDI structures, struct_list
% 'data_bigstructures' folder -for movies and big mat files
% 'plots' folder

%load prespecify settings (change path.rootpath for each computer)
clear
working = pwd;
cd C:\Users\User\Documents\MatLab\MOT
user_settings
cd(working)
%% 1. CREATE LIST OF FISH
% It will be needed for looping across fish
% It has to match the name of the .mat that cointain the 

% load(fullfile(path.data,'grouplist.mat'),'grouplist');
 
grouplist{1} = 'MOT_210318' ;
grouplist{2} = 'MOT_210322' ;
grouplist{3} = 'MOT_210326' ;
grouplist{4} = 'MOT_210405' ;
grouplist{5} = 'MOT_210526' ;
grouplist{6} = 'MOT_210604' ;
save(fullfile(path.data,'grouplist.mat'),'grouplist');

% load(fullfile(path.data,'grouplist.mat'),'grouplist');
%% 2. CREATE VSDI structure (FOR EACH FISH) 
% VSDI.info.Sonset = 590; % ms (from start of recording) %@ SET
VSDI.ref = 210604 ; %@ SET
VSDI.info.stime = 5; %ms (sampling time) @ SET
VSDI.info.Sonset = 300; % ms (from start of recording) %@ SET
VSDI.info.Sdurat = 150; % ms (from start of recording) %@ SET
VSDI.expref = 'MOT'; 

% IMPORT LIST. Have to be saved from Brainvision. See notes for details 
listpath =  path.list; %@ SET in user_settings
listpath = fullfile(listpath,strcat('filelist',num2str(VSDI.ref),'.csv'));

list_table=  readtable(listpath);
VSDI.list = table2struct(list_table);
for triali = 1:length(VSDI.list)
VSDI.trialref(triali,1) = str2num(VSDI.list(triali).Name(1:end-5)); %save references for the trials
end

%Get the preceeding ITI for each trial
[VSDI.iti] = get_iti(VSDI); 
MOT('save',VSDI);

%% @@@@@@@@ SET
% ADD MANUALLY THE CONDITION FOR EACH TRIAL (in new fields -name them to be able to copy them)...
% Non-included trials: NaN
disp('MANUALLY IDENTIFY CONDITIONS')
return

% Set all to NaN and later add the conditions
for ii = 1:length(VSDI.list)
    VSDI.list(ii).code= NaN;
end

% BLANK = 0; BOCA = 1; LOMO = 2; COLA = 3;
for ii = [36:44 50:57 84:91 113:119 ]
    VSDI.list(ii).code= 3;
end


MOT('save', VSDI);

%AND THEN COPY IT INTO NEW FIELDS (as many as condition - columns you have)

for triali = 1:length(VSDI.list)
VSDI.condition(triali,1) = VSDI.list(triali).code; %@ SET the name of the field so it c an be c opied
end

VSDI.conditionlabels {1,1} = 0; VSDI.conditionlabels {1,2} = 'blank'; 
VSDI.conditionlabels {2,1} = 1; VSDI.conditionlabels {2,2} = 'boca'; 
VSDI.conditionlabels {3,1} = 2; VSDI.conditionlabels {3,2} = 'lomo'; 
VSDI.conditionlabels {4,1} = 3; VSDI.conditionlabels {4,2} = 'cola'; 

VSDI.nonanidx= find(~isnan(VSDI.condition(:,1))) ;

% Get BV times from 'Date' info
for ii = 1:length(VSDI.list)
    date = VSDI.list(ii).Date; %get from VSDI structure
    date.Format = 'HH:mm:ss'; % capital 'H' so it's in 24h format
    date= cellstr(date); %turn into string to keep only hours
    hour = duration(date, 'Format','hh:mm:ss'); % turn again into duration
   
    %store in variable
    VSDI.trialtime(ii).trialref = VSDI.trialref(ii);
    VSDI.trialtime(ii).hour = date;
    VSDI.trialtime(ii).trialref = VSDI.trialref(ii);
end

MOT('save',VSDI)


% % VSDI.info.Sdur = 200; % ms(duration of stimulus) %@ SET
for triali = 1:length(VSDI.list)
VSDI.info.Sdur(triali,1) = 0.15; %@ SET in 's' the name of the field so it c an be c opied
end

% Save changes
MOT('save',VSDI);

% test saving
clear
[VSDI] = MOT('load',2);
