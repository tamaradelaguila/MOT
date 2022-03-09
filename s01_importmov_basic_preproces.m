%% s01 IMPORT MOVIES and BASIC PREPROCESSING: FOR EACH FISH
% 00 - IMPORT RAW MOVIE (from python-extracted matfiles) MATFILES ('_00raw' )
% 01 - BRAIN ALINEATION THROUGH TRIALS ('_01registered')
% 02 - DIFFERENTIAL VALUES ('_02diff')
% 03 - MASKED MOVIES ('_03crop' )

clear
working = pwd;
cd C:\Users\User\Documents\MatLab\MOT
user_settings
cd(working)

% nfish = 1;%@ SET
% [VSDI] = MOT('load',nfish);

%% PYTHON dml extraction
% adjust paths in '00automatic_extraction.py) and execute

%% 00 - IMPORT RAW MOVIE (from python-extracted matfiles) MATFILES
% folder with all movies extracted with python
% VSDmov.data (x*y*t*trial)

for nfish = [2:6]%@ SET
[VSDI] = MOT('load',nfish);
clearvars -except VSDI nfish path
py_output = ['C:\Users\User\Documents\MatLab\MOT\data\temp_spyder\MOT', num2str(VSDI.ref)]; %@SET where matfiles have been output from python %@ SET

[movies4D, times] = assemble4Draw(VSDI,py_output,VSDI.info.stime);
movieref = '_00raw';

%Save movie structure
VSDmov.ref = VSDI.ref;
VSDmov.expref = 'MOT'; 
VSDmov.movieref= movieref;
VSDmov.data = movies4D;
VSDmov.times = times;
VSDmov.hist{1} = 'raw: py-imported + assemble4Dmovies';
MOT('savemovie', VSDmov, movieref);

% [Save times  in VSDI]
VSDI.timeabs = VSDmov.times;
VSDI.timebase = VSDI.timeabs - VSDI.info.Sonset;
MOT('save',VSDI);

% Test
VSDmov = MOT('loadmovie',nfish,movieref);
disp(['fish ' num2str(VSDI.ref) ' finished'])
end
blob()
%% 01 - BRAIN ALINEATION THROUGH TRIALS
nfish = 6;%@ SET
[VSDI] = MOT('load',nfish);

clearvars -except VSDI nfish path
% 1. REFERENCES for input/output movies (see 'z_notes.txt', point 5 for
% complete list)
inputRef =  '_00raw';
outputRef = '_01registered';

%load input movie (non-aligned raw)
[inputStruct] = MOT('loadmovie',nfish,inputRef);
rawmov=inputStruct.data;

% 2. PERFORM COMPUTATIONS:  REGISTER (monomodal)
% ref_frame = rawmov(:,:,1,VSDI.nonanidx(1)) ; %background from 1st included trial
% [registermov] = register_wrap(rawmov,ref_frame, VSDI.nonanidx);

% (ALTERNATIVELY) PERFORM COMPUTATIONS:  REGISTER ()
ref_frame = rawmov(:,:,1,VSDI.nonanidx(1)) ; %background from 1st included trial
[registermov] = register_wrap2multimodal(rawmov,ref_frame, VSDI.nonanidx);

blob()

% @MANUAL VISUAL INSPECTION OF THE RESULT

disp('MANUALLY CHECK REGISTERED BACKGROUNDS')
return

for triali= makeRow(VSDI.nonanidx) %150:VSDI.nonanidx(end)
    imagesc(squeeze(registermov(:,:,1,triali))); colormap('bone')
    title(['trial=' num2str(triali)])
    pause
end

% [Save reference frame in VSDI]
VSDI.info.register_ref = ref_frame;
MOT('save',VSDI,nfish);

% 3.SAVE NEW MOVIE STRUCTURE: Save new registered movie, copying some references from the movie
% structure used to apply new changes in
VSDmov.ref = inputStruct.ref;
VSDmov.expref = 'MOT'; 
VSDmov.movieref= outputRef;
VSDmov.data = registermov;
VSDmov.times = inputStruct.times;
VSDmov.hist = inputStruct.hist;
VSDmov.hist{length(VSDmov.hist)+1,1} = 'register'; %append a new cell with new info
MOT('savemovie', VSDmov, VSDmov.movieref);
clear inputStruct


% Save background value in VSDI
for triali = makeRow(VSDI.nonanidx) %import only included trials
    VSDI.backgr(:,:,triali) = VSDmov.data(:,:,1,triali); % store background
end

MOT('save',VSDI,nfish);


%% 02 - DIFFERENTIAL VALUES F0= 10 preS
clearvars -except path
user_settings

for nfish = [2:6]
clearvars -except nfish

[VSDI] = MOT('load',nfish);

clearvars -except VSDI nfish


% 1. REFERENCES for input/output movies
inputRef =  '_01registered';
outputRef = '_02diff_f0pre';

inputStruct = MOT('loadmovie', nfish, inputRef);
% [VSDmov] = MOT('loadmovie',nfish,outputRef); % BORRAR

% 2. PERFORM COMPUTATIONS: %DIFFERENTIAL VALUES

inputdata = inputStruct.data;

n_preSframes = 10;  % @SET! nº of frames pre-Stimulus to use as F0 in differential formula

% Turn into string to save later in History:
baseltext = strcat(num2str(n_preSframes),'frames_preS');


% Preallocate in NaN
inputdim = size(inputdata);
diffmovies = NaN(inputdim(1),inputdim(2),inputdim(3)+1,inputdim(4));
F0 = NaN(inputdim(1),inputdim(2),inputdim(4)); 

for triali = makeRow(VSDI.nonanidx) %import only included trials
    inputmovie = squeeze(inputdata(:,:,:,triali));
    
    idx_preS = dsearchn(VSDI.timeabs,VSDI.info.Sonset); %
    ...if there are multiple Sonset, change 'VSDI.info.Sonset(:,triali)'
        
% GET IDX OF FRAMES PREVIOUS TO THE
b(1) = idx_preS-n_preSframes;
b(2) = idx_preS-1;
baseframe = b(1):b(2);

[diffmovies(:,:,:,triali), F0(:,:,triali)] = raw2diff(inputmovie, baseframe);

% VSDI.backgr(:,:,triali) = diffmovies(:,:,end,triali); % store background
end

% 3.SAVE NEW MOVIE STRUCTURE:  copying some references from the movie
% structure used to apply new changes in
VSDmov.ref = inputStruct.ref;
VSDmov.expref = 'MOT'; 
VSDmov.movieref= outputRef;
VSDmov.data = diffmovies;
VSDmov.times = inputStruct.times;
VSDmov.hist = inputStruct.hist;
VSDmov.hist{length(VSDmov.hist)+1,1} = ['raw2diff_' baseltext]; %append a new cell with new info
VSDmov.F0 = F0;

MOT('savemovie', VSDmov, VSDmov.movieref);

% SUGGESTION: if different F0 are ,keep the basic reference + info about the F0, e.g. outputRef = '_02diffbase10';
end 
blob()
 %% OLD 02 - DIFFERENTIAL VALUES !!!!!!!!!!!!!!!!!!!!!!!!! DO NOT USE
% clearvars -except VSDI nfish 
% 
% % 1. REFERENCES for input/output movies
% inputRef =  '_01registered';
% % outputRef = '_02diff_f0pre';
% 
% inputStruct = MOT('loadmovie', nfish, inputRef);
% 
% % 2. PERFORM COMPUTATIONS: %DIFFERENTIAL VALUES
% 
% inputdata = inputStruct.data;
% 
% baseframe = 1; % @SET! idx of frames to use as F0 in differential formula
% % Turn into string to save later in History:
% baseltext = strcat(num2str(baseframe(1)),'to',num2str(baseframe(end)));
% 
% % Preallocate in NaN
% inputdim = size(inputdata);
% diffmovies = NaN(inputdim(1),inputdim(2),inputdim(3)+1,inputdim(4));
% 
% for triali = makeRow(VSDI.nonanidx) %import only included trials
%     inputmovie = squeeze(inputdata(:,:,:,triali));
%     diffmovies(:,:,:,triali) = raw2diff(inputmovie, baseframe);
%     
%     VSDI.backgr(:,:,triali) = diffmovies(:,:,end,triali); % store background
% end
% 
% % 3.SAVE NEW MOVIE STRUCTURE:  copying some references from the movie
% % structure used to apply new changes in
% VSDmov.ref = inputStruct.ref;
% VSDmov.expref = 'MOT'; 
% VSDmov.movieref= outputRef;
% VSDmov.data = diffmovies;
% VSDmov.times = inputStruct.times;
% VSDmov.hist = inputStruct.hist;
% VSDmov.hist{length(VSDmov.hist)+1,1} = 'raw2diff'; %append a new cell with new info
% MOT('savemovie', VSDmov, VSDmov.movieref);
% 
% MOT('save', VSDI);
% blob()
% % SUGGESTION: if different F0 are ,keep the basic reference + info about the F0, e.g. outputRef = '_02diffbase10';


%% 03 - PERCENT DIFFERENTIAL VALUES
clearvars -except VSDI nfish

% 1. REFERENCES for input/output movies
inputRef =  '_01registered';
outputRef = '_03diff_perc';

inputStruct = MOT('loadmovie', nfish, inputRef);

% 2. PERFORM COMPUTATIONS: %DIFFERENTIAL VALUES

inputdata = inputStruct.data;

baseframe = 1:10; % @SET! idx of frames to use as F0 in differential formula
% Turn into string to save later in History:
baseltext = strcat(num2str(baseframe(1)),'to',num2str(baseframe(end)));

% Preallocate in NaN
inputdim = size(inputdata);
diffmovies = NaN(inputdim(1),inputdim(2),inputdim(3)+1,inputdim(4));

for triali = makeRow(VSDI.nonanidx) %import only included trials
    inputmovie = squeeze(inputdata(:,:,:,triali));
    diffmovies(:,:,:,triali) = raw2diffperc2(inputmovie, baseframe);
    
    VSDI.backgr(:,:,triali) = diffmovies(:,:,end,triali); % store background
    disp(triali)
end

VSDI.F0 = mean(inputdata(:,:,baseframe,:),3); %store F0


% 3.SAVE NEW MOVIE STRUCTURE:  copying some references from the movie
% structure used to apply new changes in
VSDmov.ref = inputStruct.ref;
VSDmov.expref = 'MOT'; 
VSDmov.movieref= outputRef;
VSDmov.data = diffmovies;
VSDmov.times = inputStruct.times;
VSDmov.hist = inputStruct.hist;
VSDmov.hist{length(VSDmov.hist)+1,1} = [outputRef baseltext]; %append a new cell with new info
MOT('savemovie', VSDmov, VSDmov.movieref);

MOT('save', VSDI);
blob()
% SUGGESTION: if different F0 are ,keep the basic reference + info about the F0, e.g. outputRef = '_02diffbase10';

%% 04 - [1] MAKE CROPMASK, [2] CROPPED-MOVIES (optional)
clearvars -except path

nfish = 4; 
VSDI = MOT('load', nfish);
% MAKE MASK FROM REFERENCE FRAME AND SAVE IN VSDI
ref_frame = VSDI.backgr(:,:,VSDI.nonanidx(1)); %the background from the first included trial

%%  Before cropping: check all backgrounds to take into account if there is
%  much movements (and, for instance, leave out of the mask the margins)
for triali = makeRow(VSDI.nonanidx)
    imagesc(VSDI.backgr(:,:,triali)); colormap('bone'); axis image
    title(strcat('trial=',num2str(triali)))
    pause %to advance to the next frame, press any key; to skip to the end, press 'Ctrl+C'
end
close 

% alternative to look at them all at the same time in tiled figures
inspect_allbackgrounds(VSDI.backgr)

%% DRAW & SAVE CROPMASK:
[crop_poly, crop_mask] = roi_draw(ref_frame);
[crop_poly, crop_mask] = roi_draw(VSDI.backgr(:,:,50));

%% View the result on a all trial
for trialsel= makeRow(VSDI.nonanidx)% makeRow(VSDI.nonanidx) %@ SET (if you want to check the mask onto any specific frame)
    roi_preview(VSDI.backgr(:,:,trialsel), crop_poly{1});
    title(['trial' num2str(trialsel)])
    pause
    close
end

    roi_preview(VSDI.backgr(:,:,10), crop_poly{1}); % en un solo trial
    imagesc(crop_mask)
    

%% IF WE ARE HAPPY WITH THE MASK: SAVE in structure
VSDI.crop.mask = crop_mask;
VSDI.crop.poly = crop_poly{1}; %stored in rows

% Save one cropped image to help in drawing the different ROI
cropframe =roi_crop(ref_frame, VSDI.crop.mask);
% imagesc(cropframe); colormap bone

VSDI.crop.preview = cropframe; % save the crop movie of the first included frame
MOT('save', VSDI);

% [2] CROPPED-MOVIES (mute if you don't want to extract them)
%
% % 1. REFERENCES for input/output movies
% inputRef =  '_02diff';
% % inputRef =  '_03diff_perc';
% outputRef = '_04crop';
%
% inputStruct = MOT('loadmovie', nfish, inputRef);
% inputdata = inputStruct.data;
%
% % 2. PERFORM COMPUTATIONS: APPLY MASK
% cropmovies = NaN (size(inputdata));
%
% %initialize
% for triali = makeRow(VSDI.nonanidx)
% inputmovie = squeeze(inputdata(:,:,:,triali));
% cropmovies(:,:,:,triali)= roi_crop(inputmovie, VSDI.crop.mask);
% end
%
% % 3.SAVE NEW MOVIE STRUCTURE:  copying some references from the movie
% % structure used to apply new changes in
% VSDmov.ref = inputStruct.ref;
% VSDmov.expref = 'MOT'; 
% VSDmov.movieref= outputRef;
% VSDmov.data = cropmovies;
% VSDmov.times = inputStruct.times;
% VSDmov.hist = inputStruct.hist;
% VSDmov.hist{length(VSDmov.hist)+1,1} = 'crop_backgr'; %append a new cell with new info
%
% MOT('savemovie', VSDmov, VSDmov.movieref);

%% LAST UPDATE: 
% 09/02/22 - in MOT experiment
% 21/11/21 - add F0 to the function raw2diff and update the code to save it
% in the movies VSDmov
% 09/02/22 - _16diff_f0pre is now  _02diff_f0pre

