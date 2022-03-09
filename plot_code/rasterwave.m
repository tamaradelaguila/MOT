clear
working = pwd;
cd C:\Users\User\Documents\MatLab\MOT
user_settings
cd(working)

%----------------------------------------------------------------
... @SET: fish + conditions + parameters
    %----------------------------------------------------------------
nfish = 2
conditions = [1:4];

movie_ref = '_05filt2'; % '_04filt1' '_05filt2'

saveraster = 1;
rasterclim = [0 0.003] ; %when [], the script is useful for outliers identification
savein = 'C:\Users\User\Documents\MatLab\MOT\plot\rasterwaves';

roi2plot = 'dm4m_c';

roikind = 'circle'; %

trange = [-300 1300]; %ms Range of analysis

%----------------------------------------------------------------
... @SET: MEASURE (OR LOOP THROUGH ALL MEASURES)
    %----------------------------------------------------------------
reject_on= 1;

setting.manual_visual = 1; %
% setting.manual_reject = 1; %
% setting.GSmethod_reject = 1;  %
% setting.GSabsthres_reject = 1; %
% setting.force_include = 0; %

%% LOAD / COMPUTE SETTINGS
%----------------------------------------------------------------
... LOAD DATA
    %----------------------------------------------------------------
VSDI = MOT('load', nfish);
VSDmov = MOT('loadmovie',nfish,movie_ref);

%----------------------------------------------------------------
... GET INDEXES OF TIMERANGE AND ADJUSTED TIMEBASE
    %----------------------------------------------------------------
idxrange = dsearchn(makeCol(VSDI.timebase), makeCol(trange));
idxrange = idxrange(1) : idxrange(end); % robust code in case we input both range or two-values

idx0 = dsearchn(makeCol(VSDI.timebase), 0);
timebase_adj = VSDI.timebase(idxrange);

%----------------------------------------------------------------
... COMPUTE REJECTION IDX FROM REJECT-OPTIONS
    %----------------------------------------------------------------
rej = 'reject' ;
if reject_on > 1
    rej = [rej num2str(reject_on)];
end

rejectidx = [];

if setting.manual_visual
    rejectidx = [rejectidx  makeRow(VSDI.(rej).visual)];
end
%
% if setting.manual_reject
%     rejectidx = [rejectidx  makeRow(VSDI.reject.manual)];
% end
%
% if setting.GSabsthres_reject
%     rejectidx = [rejectidx  makeRow(VSDI.(rej).GSabs025)];
%
% end
%
% if setting.GSmethod_reject
%     rejectidx = [rejectidx makeRow(VSDI.(rej).GSdeviat2sd)];
%
% end
%
% if setting.force_include
%     rejectidx = setdiff(rejectidx, VSDI.reject.forcein);
%
% end

rejectidx = sort(unique(rejectidx));


%----------------------------------------------------------------
... SELECT ROI
%----------------------------------------------------------------
roi_idx = find(strcmpi(VSDI.roi.labels(:,2), roi2plot));

roilabels = VSDI.roi.labels(:,2);

switch roikind
    case 'circle'
        masks =  VSDI.roi.circle.mask;
        
    case 'anat'
        masks = VSDI.roi.manual_mask;
end

%% ----------------------------------------------------------------
... BUILD RASTER
    %----------------------------------------------------------------
for condition = makeRow(conditions)
    condi = VSDI.conditionlabels{condition,1};
    [sel_trials] = find(VSDI.condition(:,1)==condi);
    
    if reject_on
        sel_trials = setdiff(sel_trials, rejectidx);
        disp('trials rejected')
    end
    
    n = numel(timebase_adj)-1;
    roiraster = NaN(numel(sel_trials), n,  numel(roi_idx));
        
    
    roimask = masks(:,:,roi_idx);
    tri = 0; 
    for triali = makeRow( sel_trials)
        tri = tri+1;
        %to plot single trial
        movie2plot = squeeze(VSDmov.data(:,:,idxrange,triali));
        meanF0 = squeeze(VSDmov.F0(:,:,triali));
        
        roiraster(tri,:) =  roi_TSave_percF_roiwise(movie2plot,roimask, meanF0);

    end % for triali
    
    %% ----------------------------------------------------------------
    ... PLOT RASTER
        %------------------------------------------------------------------
    
    nroi = numel(roi_idx);
    
    figure
    h = subplot(2,1,1); 
    imagesc(roiraster);

    if isempty(rasterclim)
        set(h(ploti),'xtick',[],'ytick',[])
        displ('ATT: raster color limits adjusted to its own maximum')
    else
        set(h,'xtick',[],'ytick',[], 'clim', rasterclim)
    end
    
    xline(idx0, 'color', 'w');
    ylabel(roi2plot)
    colormap(jet)
    
    subplot(2,1,2)
    plot(timebase_adj(1:end-1), mean(roiraster), 'linewidth', 1.8, 'color', 'k')
    % colorbar
    
    if condi > 0
        savename= ['Rasterwave' num2str(VSDI.ref) movie_ref '_' VSDI.conditionlabels{condition,2} VSDI.info.Sside  '_rej'  num2str(reject_on) '.jpg'];
        titulo = [num2str(VSDI.ref) '.All trials from:' VSDI.conditionlabels{condition,2} VSDI.info.Sside];
        
    elseif condi == 0
        savename= ['Rasterwave' num2str(VSDI.ref) movie_ref '_' VSDI.conditionlabels{condition,2}  '_rej'  num2str(reject_on) '.jpg'];
        titulo = [num2str(VSDI.ref) '.All trials from:' VSDI.conditionlabels{condition,2}];
    end
    
    sgtitle(titulo)
    
    if saveraster
        %           saveas(gcf, fullfile(savein, savename), 'jpg')
        set(gcf,'units','normalized','outerposition',[0 0 1 1]) %set in total screen
        print(fullfile(savein, savename),'-r600','-djpeg') % prints it as you see them
        close
    end
    
end % for condi