clear
working = pwd;
cd C:\Users\User\Documents\MatLab\MOT
user_settings
cd(working)

%----------------------------------------------------------------
... @SET: fish + conditions + parameters
%----------------------------------------------------------------
nfish = 6;
cond_rows = 3; 

movie_ref = '_05filt2'; %

saveraster = 0;
fact_clim = 0.5; % percet respect to the maximum to set the upper colorlimit
% the script is useful for outliers identification
savein = 'C:\Users\User\Documents\MatLab\MOT\plot\rasters';

% side = 'ipsi';
side = 'contra';

if strcmpi(side, 'ipsi')
%     selroinames_ini = {'dm4m_i', 'dm4l_i', 'dm3_i', 'dm2_i', 'dm1_i', 'dld_i'};
    selroinames_ini = {'dm4m_i'}; 

elseif strcmpi(side, 'contra')
%     selroinames_ini = {'dm4m_c', 'dm4l_c', 'dm3_c', 'dm2_c', 'dm1_c', 'dld_c'};
    selroinames_ini = {'dm4m_c'};
end

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
        % ADJUST SELECTED ROI  (according to whether the roi exists in the
        % fish
        %----------------------------------------------------------------
        ii = 1;
        selroinames = [];
        for roii = 1:length(selroinames_ini)
            idx =  find(contains( VSDI.roi.labels(:,2),selroinames_ini{roii}));
            if sum(idx)>0
                selroinames{ii} = VSDI.roi.labels{idx,2};
                ii = ii+1;
            end
            clear idx
        end

%----------------------------------------------------------------
... SELECT ROI
%----------------------------------------------------------------
selroi =name2idx(selroinames, VSDI.roi.labels(:,2));
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
for condition = makeRow(cond_rows)
    condi = VSDI.conditionlabels{condition,1};
[sel_trials] = find(VSDI.condition(:,1)==condi);

if reject_on
    sel_trials = setdiff(sel_trials, rejectidx);
    disp('trials rejected')
end

n = numel(timebase_adj)-1;
roiraster = NaN(numel(sel_trials), n,  numel(selroi));

ii = 0; %counter for all rasters in 2D

for roii = makeRow(selroi)
    
    roimask = masks(:,:,roii);
    
    tri = 0;
    for triali = makeRow( sel_trials)
        tri = tri +1;
        ii = ii+1;
        %to plot single trial
        movie2plot = squeeze(VSDmov.data(:,:,idxrange,triali));
        meanF0 = squeeze(VSDmov.F0(:,:,triali));
        
        roiraster(tri,:,roii) =  roi_TSave_percF_roiwise(movie2plot,roimask, meanF0);
        rasters2D (ii,1:n) =  roi_TSave_percF_roiwise(movie2plot,roimask, meanF0);
    end % for triali
    ii = ii +1; 
    rasters2D (ii:ii+1,:, condi+1) =  NaN(2,n);
    ii = ii +1; 
end % for roii

%% ----------------------------------------------------------------
... PLOT RASTER
%------------------------------------------------------------------

nroi = numel(selroi);
condmax = max(roiraster(:)); %maximum of that condition
% to calculate as % of maximum value
upperlim = condmax*fact_clim;
% to override the previous one
rasterclim = [0 upperlim] ; %when [], clim is automatically set and 
% rasterclim = [0 0.003] % uncomment to override 

figure
ploti = 0;
for roii = makeRow(selroi)
    ploti = ploti+1;
    h(ploti) = subplot(nroi,1,ploti);
    imagesc(roiraster(:,:,roii)); 
    axis tight
    if isempty(rasterclim)
          set(h(ploti),'xtick',[],'ytick',[])
          displ('ATT: raster color limits adjusted to its own maximum')
    else
    set(h(ploti),'xtick',[],'ytick',[], 'clim', rasterclim)
    end

    xline(idx0, 'color', 'w');
    ylabel(roilabels{roii})
    colormap(jet)
    colorbar
end
% colorbar
        
        if condi > 0
            savename= ['RoiRasters_' num2str(VSDI.ref) movie_ref '_' VSDI.conditionlabels{condition,2} VSDI.info.Sside  '_rej'  num2str(reject_on) '_H' side '.jpg'];
            titulo = [num2str(VSDI.ref) '.All trials from:' VSDI.conditionlabels{condition,2} VSDI.info.Sside '. H' side];

        elseif condi == 0
            savename= ['RoiRasters_' num2str(VSDI.ref) movie_ref '_' VSDI.conditionlabels{condition,2}  '_rej'  num2str(reject_on) '_H' side '.jpg'];
            titulo = [num2str(VSDI.ref) '.All trials from:' VSDI.conditionlabels{condition,2} '. H' side];
        end
        
                sgtitle(titulo)

        if saveraster
%           saveas(gcf, fullfile(savein, savename), 'jpg')
            set(gcf,'units','normalized','outerposition',[0 0 1 1]) %set in total screen
            print(fullfile(savein, [savename]),'-r600','-djpeg') % prints it as you see them
            close
        end

end % for condi