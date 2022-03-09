%% s03 DEFINE ROIS

clear
working = pwd;
cd C:\Users\User\Documents\MatLab\MOT
user_settings
cd(working)


nfish = 4;
VSDI = MOT('load',nfish);

% Define roi (for each case, R or L somatic stim)
if strcmpi(VSDI.info.Sside, 'R')
    
VSDI.roi.labels = {...
   'dm4m_R', 'dm4m_L',...
   'dm4l_R', 'dm4l_L',...
   'dm3_R', 'dm3_L',...
   'dm2_R', 'dm2_L',...
   'dldm_L',...
   'dm1_R', 'dm1_L'}';

elseif strcmpi(VSDI.info.Sside, 'L')
    VSDI.roi.labels = {...
   'dm4m_R', 'dm4m_L',...
   'dm4l_R', 'dm4l_L',...
   'dm3_R', 'dm3_L',...
   'dm2_R', 'dm2_L',...
   'dldm_R',...
   'dm1_R', 'dm1_L'}';
end

for ii = 1:size(VSDI.roi.labels,1)
    roiname  =VSDI.roi.labels{ii,1};
    VSDI.roi.labels{ii,2} = roiname_ipsicontra(roiname, VSDI.info.Sside); 
end

MOT('save', VSDI)

%% ANATOMICAL ROI
% Draw ROIs and store in structure
close all
% DEFINE with the help of the function:
labels = VSDI.roi.labels(:,1);
manual_poly= [];
manual_mask= [];
[manual_poly, manual_mask] = roi_draw(VSDI.crop.preview,labels); %in each column, one set of rois

% View the result
%  roi_preview_multiple(VSDI.crop.preview, manual_poly);
  roi_preview_multiple(VSDI.backgr(:,:,VSDI.nonanidx(1)), manual_poly); %see in a selected frame
%  roi_preview_multiple(VSDI.backgr(:,:,,VSDI.nonanidx(end), manual_poly); %see in a selected frame

% % TO REDO A SINGLE ROI, DRAW THEM ALL AND COPY MANUALLY THE COORDINATES
% ... AND MASK
%         [redo_poly, redo_mask] = roi_draw(VSDI.crop.preview,VSDI.roi.labels); 
%          roi_preview_multiple(VSDI.backgr(:,:,VSDI.nonanidx(1)), redo_poly); %see in a selected frame
%         manual_poly{15,1} = redo_poly{1,1};
%         manual_mask(:,:,15:16) = redo_mask(:,:,1:2);
%         roi_preview_multiple(VSDI.backgr(:,:,VSDI.nonanidx(end)), manual_poly); %see in a selected frame

VSDI.roi.manual_poly  = manual_poly;
VSDI.roi.manual_mask = manual_mask; 
MOT('save',VSDI);

roi_preview_multiple(VSDI.backgr(:,:,VSDI.nonanidx(1)),  VSDI.roi.manual_poly); %see in a selected frame

%% STORE IN WAVES-STRUCTURE that we create here (for extraction in '_extract_ROItimeseries):

% VSDroiTS.ref = VSDI.ref; 
% VSDroiTS.roi = VSDI.roi; 
% MOT('savewave', VSDroiTS);

%% CIRCULAR ROIS
% R = 10;
% [manual_poly, manual_mask] = roicir_draw(VSDI.crop.preview,VSDI.roi.labels, R); %in each column, one set of rois
% 
% figure
% imagesc(VSDI.crop.preview); colormap('bone');
% drawnroi = images.roi.Circle('InteractionsAllowed', 'translate','Radius',R, 'LineWidth',1.5);
% draw(drawnroi);
% 
% 
% figure
% imagesc(VSDI.crop.preview); colormap('bone');
% drawncircle = drawcircle('Radius',6,'InteractionsAllowed', 'translate', 'LineWidth',1.5);


% drawnroi = drawcircle([],'Radius',40, 'InteractionsAllowed', 'translate', 'LineWidth',1.5);


%% DRAW CIRCULAR ROI WITH FUNCTION
clear
user_settings
nfish = 4;
VSDI = MOT('load',nfish);


% close all
r = 4;
labels = VSDI.roi.labels(:,1);
figure
[coord, mask] = roicir_draw(VSDI.crop.preview,VSDI.roi.labels,r); 

circle.center = coord;
circle.R = r;
circle.mask = mask;

figure
%  roicirc_preview_multiple(VSDI.crop.preview, circle.center, circle.R); 
 roicirc_preview_multiple(VSDI.backgr(:,:,VSDI.nonanidx(end)), circle.center, circle.R); 

imagesc(circle.mask(:,:,1))
axis image

VSDI.roi.circle.center = circle.center;
VSDI.roi.circle.R = circle.R;
VSDI.roi.circle.mask = circle.mask;
MOT('save',VSDI)

% 
% %view one by one
% for roi = 1:length(VSDI.roi.labels)
%      roicirc_preview_multiple(VSDI.crop.preview, VSDI.roi.circle.center(roi,:), VSDI.roi.circle.R); 
% title([num2str(VSDI.ref) ':' VSDI.roi.labels{roi}])
% pause
% end

%% TO CHECK ALL ROI INDIVIDUALLY
% TO CHECK ALL ROI
nfish = 6;
VSDI = MOT('load', nfish);

disp(['side of somatic side:' VSDI.info.Sside])
for roi = 1:length(VSDI.roi.labels)
     roicirc_preview_multiple(VSDI.crop.preview, VSDI.roi.circle.center(roi,:), VSDI.roi.circle.R); 
title([num2str(VSDI.ref) ':' VSDI.roi.labels{roi}])
pause
close
end

disp(['side of somatic side:' VSDI.info.Sside])
for roi = 1:length(VSDI.roi.labels)
roi_preview(VSDI.backgr(:,:,VSDI.nonanidx(1)),  VSDI.roi.manual_poly{roi}); %see in a selected frame
title([num2str(VSDI.ref) ':' VSDI.roi.labels{roi}])
pause
close
end

%% ADD NEW CIRCULAR ROI TO THE ALREADY COLLECTION - !!! ONE BY ONE
clear
user_settings
nfish = 12;
VSDI = MOT('load',nfish);

newroi_label = {'dldr_R2'};

% newroi_label = {'dld_rR'};
newroi_row = numel(VSDI.roi.labels)+1; %make sure it's an empty row

close all
r = VSDI.roi.circle.R;
[coord, mask] = roicir_draw(VSDI.crop.preview,newroi_label,r); 

circle.center = coord;
circle.R = r;
circle.mask = mask;

 roicirc_preview_multiple(VSDI.crop.preview, circle.center, circle.R); 
 roicirc_preview_multiple(VSDI.backgr(:,:,VSDI.nonanidx(end)), circle.center, circle.R); 

VSDI.roi.labels{newroi_row} =  newroi_label{1};
VSDI.roi.circle.center(newroi_row,:) = circle.center;
VSDI.roi.circle.mask(:,:,newroi_row) = circle.mask;
% MOT('save',VSDI)


