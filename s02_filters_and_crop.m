%% s02 FILTERING (further preprocessing) and crop
% Crop masks should have been drawn in the s01_importmov_basic_preproces

clear
working = pwd;
cd C:\Users\User\Documents\MatLab\MOT
user_settings
cd(working)


%% _04 FILT1:
tic
for nfish = [1 2]

clearvars -except nfish path
 
VSDI = MOT('load',nfish);

%@ SET FILTERS PARAMETERS :
tcnst = 10;% Time-constant for TEMPORAL SMOOTHING
% mean_filter = 3; %  gaussian smoothing kernel for SPATIAL SMOOTHERING
meanpix = 9;
medianpix = 3; % pixels size for MEDIAN SPATIAL FILTER


% 1. REFERENCES for input/output movies (see 'z_notes.txt', point 5 for
% complete list)
inputRef =  '_02diff_f0pre';
outputRef = '_04filt1'; %@ SET

% Load input movie
[inputStruct] = MOT('loadmovie',nfish,inputRef);
inputdata=inputStruct.data;


% 2. PERFORM COMPUTATIONS: %DIFFERENTIAL VALUES

% Preallocate in NaN
filtmov1 = NaN(size(inputdata));
filtmov2 = NaN(size(inputdata));
filtmov3 = NaN(size(inputdata));


% 2.1. Tcnst
for triali = makeRow(VSDI.nonanidx)
    tempmov = inputdata(:,:,:,triali);
    filtmov1(:,:,:,triali) = filter_Tcnst(tempmov,tcnst);
    clear tempmov
    disp(triali)
end

% 2.2. Spatial Filter (mean)
for triali = makeRow(VSDI.nonanidx)
    tempmov = filtmov1(:,:,:,triali);
    filtmov2(:,:,:,triali) = filter_spatial2(tempmov, meanpix);
    clear tempmov
end

% 2.3. Median spatial filter
for triali = makeRow(VSDI.nonanidx)
    tempmov = filtmov2(:,:,:,triali);
    filtmov3(:,:,:,triali) = filter_median(tempmov, medianpix);
    clear tempmov
end

% % 2.5. Crop background
% for triali = makeRow(VSDI.nonanidx)
%     tempmov = filt4(:,:,:,triali);
%     filtmov(:,:,:,triali)= roi_crop(tempmov, VSDI.crop.mask);
%     clear tempmov 
% end

% SET definitive filtered movie that will be stored (do not forget to set the filters
% accordingly
filtmov_def = filtmov3;

% 3.SAVE NEW MOVIE STRUCTURE:  copying some references from the movie
% structure used to apply new changes in
VSDmov.ref = inputStruct.ref;
VSDmov.expref = 'MOT'; 
VSDmov.movieref= outputRef;
VSDmov.data = filtmov_def;
VSDmov.times = inputStruct.times;
%@ SET !!! according to the filters applie (append as many as needeed)
VSDmov.hist = inputStruct.hist;
VSDmov.hist{length(VSDmov.hist)+1,1} = ['tcnst = ' num2str(tcnst)]; %@ SET
VSDmov.hist{length(VSDmov.hist)+1,1} = ['spatialmean =' num2str(meanpix)]; %@ SET
VSDmov.hist{length(VSDmov.hist)+1,1} = ['median ='  num2str(medianpix)]; %@ SET
VSDmov.F0 = inputStruct.F0; 
% VSDmov.hist{length(VSDmov.hist)+1,1} = 'crop-background'; %@ SET
MOT('savemovie', VSDmov, VSDmov.movieref);

disp(['fish ' num2str(VSDI.ref) ' finished'])

end
t = toc
blob()

%% _05 FILT1:
tic
for nfish = 1:2

clearvars -except nfish path
 
VSDI = MOT('load',nfish);

%@ SET FILTERS PARAMETERS :
tcnst = 10;% Time-constant for TEMPORAL SMOOTHING
% mean_filter = 3; %  gaussian smoothing kernel for SPATIAL SMOOTHERING
meanpix = 9;
medianpix = 3; % pixels size for MEDIAN SPATIAL FILTER
cutoff = 0.1 ; %Hz for the butterworth high-pass BLEACH-LIKE FILTER 


% 1. REFERENCES for input/output movies (see 'z_notes.txt', point 5 for
% complete list)
inputRef =  '_02diff_f0pre';
outputRef = '_05filt2'; %@ SET

% Load input movie
[inputStruct] = MOT('loadmovie',nfish,inputRef);
inputdata=inputStruct.data;


% 2. PERFORM COMPUTATIONS: %DIFFERENTIAL VALUES

% Preallocate in NaN
filtmov1 = NaN(size(inputdata));
filtmov2 = NaN(size(inputdata));
filtmov3 = NaN(size(inputdata));
filtmov4 = NaN(size(inputdata));


% 2.1. Tcnst
for triali = makeRow(VSDI.nonanidx)
    tempmov = inputdata(:,:,:,triali);
    filtmov1(:,:,:,triali) = filter_Tcnst(tempmov,tcnst);
    clear tempmov
    disp(triali)
end

% 2.2. Spatial Filter (mean)
for triali = makeRow(VSDI.nonanidx)
    tempmov = filtmov1(:,:,:,triali);
    filtmov2(:,:,:,triali) = filter_spatial2(tempmov, meanpix);
    clear tempmov
end

% 2.3. Median spatial filter
for triali = makeRow(VSDI.nonanidx)
    tempmov = filtmov2(:,:,:,triali);
    filtmov3(:,:,:,triali) = filter_median(tempmov, medianpix);
    clear tempmov
end

% 2.4. Slow drift removal
for triali = makeRow(VSDI.nonanidx)
    tempmov = filtmov3(:,:,:,triali);
    filtmov4(:,:,:,triali) = filter_bleach_butterhigh(tempmov, cutoff, VSDI.info.stime);
    clear tempmov
end

% SET definitive filtered movie that will be stored (do not forget to set the filters
% accordingly
filtmov_def = filtmov4;

% 3.SAVE NEW MOVIE STRUCTURE:  copying some references from the movie
% structure used to apply new changes in
VSDmov.ref = inputStruct.ref;
VSDmov.expref = 'MOT'; 
VSDmov.movieref= outputRef;
VSDmov.data = filtmov_def;
VSDmov.times = inputStruct.times;
%@ SET !!! according to the filters applie (append as many as needeed)
VSDmov.hist = inputStruct.hist;
VSDmov.hist{length(VSDmov.hist)+1,1} = ['tcnst = ' num2str(tcnst)]; %@ SET
VSDmov.hist{length(VSDmov.hist)+1,1} = ['spatialmean =' num2str(meanpix)]; %@ SET
VSDmov.hist{length(VSDmov.hist)+1,1} = ['median ='  num2str(medianpix)]; %@ SET
VSDmov.F0 = inputStruct.F0; 
VSDmov.hist{length(VSDmov.hist)+1,1} = ['butterhigh=' num2str(cutoff) 'Hz']; %@ SET
MOT('savemovie', VSDmov, VSDmov.movieref);

disp(['fish ' num2str(VSDI.ref) ' already filtered'])

end
t = toc
blob()

%% _18filt6 (from diff_f0preS-10frames) + bleach removal
% clear
% user_settings;
% 
% for nfish =  [2]
%     user_settings;
%     
%     VSDI = MOT('load',nfish);
%     
%     %@ SET FILTERS PARAMETERS :
%     tcnst = 10;% Time-constant for TEMPORAL SMOOTHING
%     % mean_filter = 3; %  gaussian smoothing kernel for SPATIAL SMOOTHERING
%     meanpix = 9;
%     medianpix = 3; % pixels size for MEDIAN SPATIAL FILTER
%     cutoff = 0.1 ; %Hz for the butterworth high-pass BLEACH-LIKE FILTER 
%     
%     % 1. REFERENCES for input/output movies (see 'z_notes.txt', point 5 for
%     % complete list)
%     inputRef =  '_16diff_f0pre';
%     outputRef = '_18filt6'; %@ SET
%     
%     % Load input movie
%     [inputStruct] = MOT('loadmovie',nfish,inputRef);
%     inputdata=inputStruct.data;
%     
%     
%     % 2. PERFORM COMPUTATIONS: %DIFFERENTIAL VALUES
%     
%     % Preallocate in NaN
%     filtmov1 = NaN(size(inputdata));
%     filtmov2 = NaN(size(inputdata));
%     filtmov3 = NaN(size(inputdata));
%     filtmov4 = NaN(size(inputdata));
%     
%     
%     % 2.1. Tcnst
%     for triali = makeRow(VSDI.nonanidx)
%         tempmov = inputdata(:,:,:,triali);
%         filtmov1(:,:,:,triali) = filter_Tcnst(tempmov,tcnst);
%         clear tempmov
%         disp(triali)
%     end
%     
%     % 2.2. Spatial Filter (mean)
%     for triali = makeRow(VSDI.nonanidx)
%         tempmov = filtmov1(:,:,:,triali);
%         filtmov2(:,:,:,triali) = filter_spatial2(tempmov, meanpix);
%         clear tempmov
%     end
%     
%     % 2.3. Median spatial filter
%     for triali = makeRow(VSDI.nonanidx)
%         tempmov = filtmov2(:,:,:,triali);
%         filtmov3(:,:,:,triali) = filter_median(tempmov, medianpix);
%         clear tempmov
%     end
%     
%     % 2.4. Slow drift removal
%     for triali = makeRow(VSDI.nonanidx)
%         tempmov = filtmov3(:,:,:,triali);
%         filtmov4(:,:,:,triali) = filter_bleach_butterhigh(tempmov, cutoff, VSDI.info.stime);
%         clear tempmov
%     end
%     % % 2.5. Crop background
%     % for triali = makeRow(VSDI.nonanidx)
%     %     tempmov = filt4(:,:,:,triali);
%     %     filtmov(:,:,:,triali)= roi_crop(tempmov, VSDI.crop.mask);
%     %     clear tempmov
%     % end
%     
%     % SET definitive filtered movie that will be stored (do not forget to set the filters
%     % accordingly
%     filtmov_def = filtmov4;
%     
%     % 3.SAVE NEW MOVIE STRUCTURE:  copying some references from the movie
%     % structure used to apply new changes in
%     VSDmov.ref = inputStruct.ref;
%     VSDmov.expref = 'MOT'; 
%     VSDmov.movieref= outputRef;
%     VSDmov.data = filtmov_def;
%     VSDmov.times = inputStruct.times;
%     %@ SET !!! according to the filters applie (append as many as needeed)
%     VSDmov.hist = inputStruct.hist;
%     VSDmov.hist{length(VSDmov.hist)+1,1} = ['tcnst =' num2str(tcnst)]; %@ SET
%     VSDmov.hist{length(VSDmov.hist)+1,1} = ['spatialmean = ' num2str(meanpix)]; %@ SET
%     VSDmov.hist{length(VSDmov.hist)+1,1} = ['median =' num2str(medianpix)]; %@ SET
%     VSDmov.hist{length(VSDmov.hist)+1,1} = ['butterhigh=' num2str(cutoff) 'Hz']; %@ SET
%     VSDmov.F0 = inputStruct.F0; %@ SET
%     
%     % VSDmov.hist{length(VSDmov.hist)+1,1} = 'crop-background'; %@ SET
%     MOT('savemovie', VSDmov, VSDmov.movieref);
%     
%     blob()
%     clear
% end
% blob(); pause(0.1); blob()
% 
% %% _19filt7 (from diff_f0preS-10frames) only bleach removal
% clear
% user_settings;
% 
% for nfish =  [2 8 9 10 11 12]
%     user_settings;
%     
%     VSDI = MOT('load',nfish);
%     
%     %@ SET FILTERS PARAMETERS :
%     cutoff = 0.1 ; %Hz for the butterworth high-pass BLEACH-LIKE FILTER 
%     
%     % 1. REFERENCES for input/output movies (see 'z_notes.txt', point 5 for
%     % complete list)
%     inputRef =  '_16diff_f0pre';
%     outputRef = '_19filt7'; %@ SET
%     
%     % Load input movie
%     [inputStruct] = MOT('loadmovie',nfish,inputRef);
%     inputdata=inputStruct.data;
%     
%     
%     % 2. PERFORM COMPUTATIONS:
%     
%     % Preallocate in NaN
%     filtmov1 = NaN(size(inputdata));
%     
%     % 2.4. Slow drift removal
%     for triali = makeRow(VSDI.nonanidx)
%         tempmov = inputdata(:,:,:,triali);
%         filtmov1(:,:,:,triali) = filter_bleach_butterhigh(tempmov, cutoff, VSDI.info.stime);
%         clear tempmov
%     end
%     
%     % SET definitive filtered movie that will be stored (do not forget to set the filters
%     % accordingly
%     filtmov_def = filtmov1;
%     
%     % 3.SAVE NEW MOVIE STRUCTURE:  copying some references from the movie
%     % structure used to apply new changes in
%     VSDmov.ref = inputStruct.ref;
%     VSDmov.expref = 'MOT'; 
%     VSDmov.movieref= outputRef;
%     VSDmov.data = filtmov_def;
%     VSDmov.times = inputStruct.times;
%     %@ SET !!! according to the filters applie (append as many as needeed)
%     VSDmov.hist = inputStruct.hist;
%     VSDmov.hist{length(VSDmov.hist)+1,1} = ['butterhigh=' num2str(cutoff) 'Hz']; %@ SET
%     VSDmov.F0 = inputStruct.F0; %@ SET
%     
%     % VSDmov.hist{length(VSDmov.hist)+1,1} = 'crop-background'; %@ SET
%     MOT('savemovie', VSDmov, VSDmov.movieref);
%     
%     blob()
%     clear
% end
blob(); pause(0.1); blob()