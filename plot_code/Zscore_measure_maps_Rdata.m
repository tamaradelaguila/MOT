%%  MEASURE MAPS + R data

% THIS CODE IS TRICKY IN THAT NOT ALL FISH HAVE THE SAME TIMEBASE, SO IT
% HAS TO BE ADJUSTED. IT IS IMPORTANT THEREFORE TO TAKE THIS INTO ACCOUNT
% WHEN USING THE BASELINE (VSDI.baseline) THAT IT HAS TO BE ALSO ADPATED

% BLANK SUBSTRACTION OF THE WAVES

% average movie
% extract wave for each pixel and make the measure maps (for all
% conditions): for condition for xi for yi

% z-spatial all conditions for each fish
% extract roi from z-spatially maps >>> definite measures

clear
working = pwd;
cd C:\Users\User\Documents\MatLab\MOT
user_settings
cd(working)

% ///////////////////////////////////////////////////////////
% SETTINGS

sourcedata = 'normal';
selroinames_ini = {'dm4m_c', 'dm4l_c','dm3_c','dm2_c','dm1_c','dldm_c','dm4m_i', 'dm4l_i','dm3_i','dm2_i','dm1_i'}; %4roi
% 
roikind = 'circle'; %
% roikind = 'anat';

ref_movie= '_05filt2';% '

activ_unit = 'diffF'; % @ MANUALLY SET (just for display purposes)
analysisref = 'group2_'; % MANUALLY SET!!! extra info for the name. Set group according to the rows of 'groupplot' selected in;   for suji  =  [1 3 4 9]
reject_on = 1; %@ SET


% @SET: WHICH OUTPUT TO PLOT
plotmaps = 0;
savemaps = 0; % if plotmaps =1
getR = 1; %RO print Rdata in long format
getTmat = 1;
saveTmat = 1;

savein = 'C:\Users\User\Documents\MatLab\MOT\plot\Zscore_measure_maps' ;%@ SET

% Time range (to fit all the waves so they span the same timewindow) !!!
% all baselines should be also but to this range
trange = [-300 1300]; %ms Range of analysis

% FUNCTION SETTINGS
feedf.window.min = [-100 100]; %'feed-function' structure
feedf.window.max = [0 600]; % where to find the max peak
feedf.window.movsum = 50; %ms
feedf.window.basel = [-100 0]; %cambiar a  -100
feedf.window.slope=50;
feedf.window.wmean=[0 350]; %ms

% feedf.noise.fr_abovenoise = 30;
% feedf.noise.SDfactor = 2;
% feedf.noise.basel = [-200 0]; %it won't be used by the function, but will be used to manually get the baseline
feedf.method = 'movsum';

slope.window = [0, 200]; %ms

% COPY FUNCTIONS SETTINGS IN CELL STRUCTURE TO OUTPUT IN EXCEL
params{1,1} = 'window.max (ms)';
params{1,2} = feedf.window.max;

params{2,1} = 'basel (ms)';
params{2,2} = feedf.window.basel;

params{3,1} = 'movsum (ms)';
params{3,2} = feedf.window.movsum;

params{4,1} = 'wmean (ms)';
params{4,2} = feedf.window.wmean;

params{4,1} = 'slope window (ms)';
params{4,2} = [slope.window(1) slope.window(2)];

    %----------------------------------------------------------------
    % @SET: REJECT SETTINGS
    %----------------------------------------------------------------
    
    % Subsettings:
    setting.visual_reject = 1;
    setting.manual_reject = 0; %@ SET CHA-CHA-CHA-CHANGEEEEEEEEE
    setting.GSmethod_reject = 0;  %@ SET
    setting.GSabsthres_reject = 0; %@ SET+
    setting.force_include = 0; %@ SET
    
    % END
% ///////////////////////////////////////////////////////////////

    

    
    %% FIRST 'suji' LOOP :   GET WAVES (for each subject and roi)
    i = 2; % counter for long format rows: the first will be the labels
    si = 0; % counter for subjects list
    fi = 0; %counter for n fish
    
    for nfish  =  1:5 %1:size(groupplot,1) SELECT included fish+condition

        si = si+1;
        fi = fi+1;

        cond = 0:3;
        
        VSDI = MOT('load', nfish);
        VSDmov = MOT('loadmovie',nfish,ref_movie);
        movies = VSDmov.data ;
        F0 = VSDmov.F0;
        
        %----------------------------------------------------------------
        % ADJUST SELECTED ROI  (according to whrether the roi exists in the
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
        
        % clear newnames
        %----------------------------------------------------------------
        % SELECT ROI
        %----------------------------------------------------------------
       selroi = [];
       selroi =name2idx(selroinames, VSDI.roi.labels(:,2));
       
       roilabels = [];
       roilabels = VSDI.roi.labels; 

        switch roikind
            case 'circle'
                masks =  VSDI.roi.circle.mask;
                
            case 'anat'
                masks = VSDI.roi.manual_mask;
        end
        
        %----------------------------------------------------------------
        % CONTROL ROI PICTURE
        %----------------------------------------------------------------
        %     roicirc_preview_multiple(VSDI.crop.preview, VSDI.roi.circle.center(selroi,:), VSDI.roi.circle.R);
        %     title([num2str(VSDI.ref) 'roi preview:' selroinames{:}])
        %     saveas(gcf, [num2str(VSDI.ref)'roipreview'], 'jpg')
        
        
        %----------------------------------------------------------------
        % GET INDEXES OF TIMERANGE
        %----------------------------------------------------------------
        idxrange = dsearchn(makeCol(VSDI.timebase), makeCol(trange));
        idxrange = idxrange(1) : idxrange(end); % robust code in case we input both range or two-values
        
        timebase_adj = VSDI.timebase(idxrange);
        
        % AND INDEXES FOR SLOPE RELATIVE TO THE ADJUSTED TIMEBASE
        
        slope.windowidx = dsearchn(timebase_adj, makeCol(slope.window));
        slope.windowidx = [slope.windowidx(1) slope.windowidx(end)];
        
        
        %----------------------------------------------------------------
        % COMPUTE REJECTION IDX FROM REJECT-OPTIONS
        %----------------------------------------------------------------
        rej = 'reject' ;
        if reject_on > 1
            rej = [rej num2str(reject_on)];
        end
        
        rejectidx = [];
        
        if setting.visual_reject
                rejectidx = [rejectidx  makeRow(VSDI.(rej).visual)];
        end

        
%         if setting.manual_reject
%             try
%                 rejectidx = [rejectidx  makeRow(VSDI.(rej).manual)];
%             catch
%                 rejectidx = [rejectidx  makeRow(VSDI.reject.manual)];
%                 disp(['rejec.manual was used for fish' num2str(VSDI.ref) 'because there is no reject' num2str(reject_on) '.manual'])
%             end
%         end
%         
%         if setting.GSabsthres_reject
%             rejectidx = [rejectidx  makeRow(VSDI.(rej).GSabs025)];
%             
%         end
%         
%         if setting.GSmethod_reject
%             rejectidx = [rejectidx makeRow(VSDI.(rej).GSdeviat2sd)];
%         end
%         
%         if setting.force_include
%             rejectidx = setdiff(rejectidx, VSDI.reject.forcein);
%             
%         end
        
        rejectidx = sort(unique(rejectidx));
        
        % -------------------------------------------
        % FIRST LOOP THROUGH CONDITIONS AND PIXELS TO GET PIXEL-MEASURES
        % -------------------------------------------
        
        ci = 0;
        for condition =  makeRow(cond)
            ci = ci+1;
            
            cond_blank = force0ending(condition);
            % -------------------------------------------
            % SELECT CASES  AND AVERAGE MOVIE
            % -------------------------------------------
            sel_trials = find(VSDI.condition(:,1)==condition);
            sel_blank = find(VSDI.condition(:,1)==cond_blank);
            
            if reject_on
                sel_trials = setdiff(sel_trials, rejectidx);
                sel_blank = setdiff(sel_blank, rejectidx);
            end
            
            % --------------------------------------------------------------------------
            % GET AVERAGE MOVIE. Use timerange set
            % --------------------------------------------------------------------------
            
            switch sourcedata
                
                case 'normal'
                    refcase = '';
                    
                    movieave = mean(movies(:,:,idxrange,sel_trials),4);
                    movieave(:,:,end+1) = NaN; % add a last NaN frame because the 'roi_TSave' function will delete the last point from each wave
                    ...(that normally corresponds to the background)
                        
                %                     case  '%F'
                %                         refcase = '%F';
                %
                %                         % Step-1: calculate %F trial-wise
                %                         for triali =1:size(movies,4)
                %                             movie = movies(:,:,idxrange,triali);
                %
                %                             movie(:,:,end+1) = NaN; % add a last NaN frame because the 'roi_TSave' function will delete the last point from each wave
                % %                             trialwave(:,triali) =   roi_TSave_percF_roiwise(movie,roimask, F0(:,:,triali))*1000; % scaling factor
                %                         end
                %                         % Step-2: average waves from selected trials
                %
                %                         % Step-3 : get non %F wave to get the color limit
                %                         % for the tile
                %                         movieave = mean(movie, 4);
                %                         localwave =  roi_TSave(movieave,roimask);
                %                         maxval = max(localwave(:));
                
                case  'blank-s'
                    refcase = 'blankS';
                    
                    movieave = mean(movies(:,:,idxrange,sel_trials),4);
                    movieave(:,:,end+1) = NaN; % add a last NaN frame because the 'roi_TSave' function will delete the last point from each wave
                    
                    movieblank = mean(movies(:,:,idxrange,sel_blank),4);
                    movieblank(:,:,end+1) = NaN; % add a last NaN frame because the 'roi_TSave' function will delete the last point from each wave
                    %                         blankwave = roi_TSave(movieblank,roimask);
                    %                         roiwave =  roi_TSave(movieave,roimask)-blankwave;
                    %             case '%deltaF blank-s'
                    %                 refcase = '%deltaF blanks';
                    %                         maxval = max(roiwave(:));
                    
            end %case
            
            
            % -------------------------------------------
            % CALCULATE MEASURES for each pixel and condition (and
            % store in 'maps')
            % -------------------------------------------
            switch sourcedata
                
                case 'normal'
                    
                    for xi = 1:size(movieave,1)
                        for yi = 1:size(movieave,2)
                            pixelwave = movieave(xi, yi,:);
                            
                            temp = devo_peak2peak(pixelwave, timebase_adj, feedf.window, [], feedf.method, 0, 0);
                            
                            maps.peak(xi,yi,ci) = temp.peakminusbasel;
                            maps.wmean(xi,yi,ci) = temp.wmean;
                            
                            idx0 = slope.windowidx(1);
                            idxend = slope.windowidx(end);
                            waveW = pixelwave(idx0:idxend);
                            slopemean = mean(diff(waveW));
                            maps.slope(xi,yi,ci) = slopemean;
                            
                            clear pixelwave waveW slopemean 

                        end % for yi
                    end % for xi
                    
                case  'blank-s'
                    
            end % switch sourcedata
            
        end %for condition
        
        %% SPATIAL Z-SCORE AMONG CONDITIONS (for each measure)
        dim = size(maps.peak);
        temp_peak = reshape(maps.peak, [dim(1)*dim(2) dim(3)]);
        temp_wmean = reshape(maps.wmean, [dim(1)*dim(2) dim(3)]);
        temp_slope = reshape(maps.slope, [dim(1)*dim(2) dim(3)]);
        
        Zpeak = zscore(maps.peak, 0, 'all');
        Zwmean = zscore(maps.wmean, 0, 'all');
        Zslope = zscore(maps.slope, 0, 'all');
        
        mapsZ.peak =reshape(Zpeak, [dim(1) dim(2) dim(3)]);
        
        mapsZ.wmean =reshape(Zwmean, [dim(1) dim(2) dim(3)]);
        
        mapsZ.slope = reshape(Zslope, [dim(1) dim(2) dim(3)]);


%         % -------------------------------------------
%         % GET MEASURES OF ROI WAVES FOR EACH TRIAL - for
%         % checking
%         % -------------------------------------------
%         !!!!! AT THE MOMENT IT CANNOT BE RUN BECAUSE ROILABELS MIGHT CHANGE FOR
%        EACH FISH 
%         % -------------------------------------------
%         if get_excel_indiv
%             ci = 0;
%             ti = 1; %trial counter
%             for condition =  makeRow(cond)
%                 ci = ci +1;
%                 
%                 sel_trials = find(VSDI.condition(:,1)==condition);
%                 if reject_on
%                     sel_trials = setdiff(sel_trials, rejectidx);
%                 end
%                 
%                 for triali = makeRow(sel_trials) %each trial will be saved in a row
%                     ti = ti+1; % trial marker independent from triali (because will loop through trials across conditions)
%                     movie2plot = movies(:,:,idxrange,triali);
%                     ri= 3; % row counter
%                     for roii = makeRow(selroi) %each roi in a column
%                         
%                         ri= ri+1; 
%                         roimask = masks(:,:, roii);
%                         %                             meanF0 = squeeze(mean(VSDmov.F0(:,:,triali),3));
%                         roiwave = roi_TSave(movie2plot,roimask);
%                         temp = devo_peak2peak(roiwave, VSDI.timebase(idxrange), feedf.window, [], feedf.method, 0, 0);
%                         
%                         idx0 = slope.windowidx(1);
%                         idxend = slope.windowidx(end);
%                         waveW = roiwave(idx0:idxend);
%                         slopemean = mean(diff(waveW));
%                         
%                         temp.slope = slopemean;
%                         clear waveW slopemean
%                         
%                         % trial identification / kind
%                         trialmeasure.peak{ti,1} = triali;
%                         trialmeasure.peak{ti,2} = VSDI.trialref(triali);
%                         trialmeasure.peak{ti,3} = condition;
%                         % measures: PEAK
%                         trialmeasure.peak{ti,ri} = round(temp.peakminusbasel,2);
%                         
%                         % trial identification / kind
%                         trialmeasure.wmean{ti,1} = triali;
%                         trialmeasure.wmean{ti,2} = VSDI.trialref(triali);
%                         trialmeasure.wmean{ti,3} = condition;
%                         % measures: WMEAN
%                         trialmeasure.wmean{ti,ri} = round(temp.wmean,2);
%                         
%                         % trial identification / kind
%                         trialmeasure.slope{ti,1} = triali;
%                         trialmeasure.slope{ti,2} = VSDI.trialref(triali);
%                         trialmeasure.slope{ti,3} = condition;
%                         % measures: PEAK
%                         trialmeasure.slope{ti,ri} = round(temp.slope,2);
%                         
%                         clear temp
%                         
%                         clear roimask roiwave meanF0
%                         labels{1,1} = 'idx';
%                         labels{1,2} = 'trial';
%                         labels{1,3} = 'condition';
%                         labels{1,ri} = roilabels{roii};
%                         
% %                         % TRIALS TABLE TO LATER REJECT BASED ON STD
% %                         % trial identification / kind
% %                         trialmat.peak(ti,1) = triali;
% %                         trialmat.peak(ti,2) = VSDI.trialref(triali);
% %                         trialmat.peak(ti,3) = condition;
% %                         % measures: PEAK
% %                         trialmat.peak(ti,ri) = round(temp.peakminusbasel,2);
% %                         
% %                         % trial identification / kind
% %                         trialmat.wmean(ti,1) = triali;
% %                         trialmat.wmean(ti,2) = VSDI.trialref(triali);
% %                         trialmat.wmean(ti,3) = condition;
% %                         % measures: WMEAN
% %                         trialmat.wmean(ti,ri) = round(temp.wmean,2);
%                     end % forroii
%                 end % for triali
%             end %for condition
%             % add labels
%             
%             trialmeasure.peak(1,1:length(labels)) = labels(1:end);
%             trialmeasure.wmean(1,1:length(labels)) = labels(1:end);
%             trialmeasure.slope(1,1:length(labels)) = labels(1:end);
% 
%             excelname = [analysisref 'Zscore_individual_trials' ref_movie '_' refcase 'rej' num2str(reject_on) '_' num2str(numel(selroi)) 'roi.xls'];
%             labels = {'idx' 'trial' 'cond'};
%             
%             % write output (new sheet for each fish
%             writecell (trialmeasure.peak, fullfile(savein , excelname), 'sheet',  [num2str(VSDI.ref) 'peak'])
%             writecell (trialmeasure.wmean, fullfile(savein,  excelname), 'sheet',  [num2str(VSDI.ref) 'wmean'])
%             writecell (trialmeasure.slope, fullfile(savein,  excelname), 'sheet',  [num2str(VSDI.ref) 'slope'])
% 
%             clear trialmeasure
%             
%         end %if get_excel_indiv
        
        
        %% EXTRACT ROI MEASURES AND STORE (for R)
        
        
        ci = 0;
        for condition =  makeRow(cond)
            ci = ci +1;
            
            if getR || getTmat
                % -------------------------------------------------------
                % CALCULATE MEASURE FOR EACH ROI AND STORE IN LONG FORMAT
                % -------------------------------------------------------
                for roi_i = makeRow(selroi)
                    roimask= masks(:,:,roi_i);
                    roipeak = sum(maps.peak(:,:,ci).*roimask) / sum(roimask) ;
                    roiwmean = sum(maps.wmean(:,:,ci).*roimask) / sum(roimask) ;
                    roislope = sum(maps.slope(:,:,ci).*roimask) / sum(roimask) ;
                    
                    roipeakZ = sum(mapsZ.peak(:,:,ci).*roimask) / sum(roimask) ;
                    roiwmeanZ = sum(mapsZ.wmean(:,:,ci).*roimask) / sum(roimask) ;
                    roislopeZ = sum(mapsZ.wmean(:,:,ci).*roimask) / sum(roimask) ;

                    % GET ROI IDX RELATIVE TO THE 'selroinames' (instead of
                    % absolute roi idx) so the contra/ipsi lateral roi
                    % share the same idx for the stats
                    roi_number = name2idx(VSDI.roi.labels{roi_i,2}, selroinames);
                
                    % FACTORS
                    longF{i,1} = MOT('who', VSDI.ref) ; %subject id
                    longF{i,2} = roi_number; % roi %% !!! NAME RESPECT TO 'selroinames' (instead of absolute roi nº)
                    longF{i,3} = roilabels{roi_i,2}; % roi
                    longF{i,4} = ci; % condition (1 = low; 2 = hi; 3 = blank)
                    % MEASURES
                    longF{i,5} = round(roipeak,2); % outputP(roi_i, condi)
                    longF{i,6} = round(roiwmean,2);%
                    longF{i,7} = round(roislope,2);%
                    
                    longF{i,8} = round(roipeakZ,2); % outputP(roi_i, condi)
                    longF{i,9} = round(roiwmeanZ,2);%
                    longF{i,10} = round(roislopeZ,2);%
                    
                    i = i+1;
                    
                    % T-TEST MATRIX
                    tmat.peakZ(roi_i, fi, ci) =  round(roipeakZ,2);
                    tmat.wmeanZ(roi_i, fi, ci) =  round(roiwmeanZ,2);
                    tmat.slopeZ(roi_i, fi, ci) =  round(roislopeZ,2);
                    
                end %for roi_i
            end % for if getR
        end % for condition
        
        %% PLOT MAPS
        if plotmaps
            
            savename = [num2str(VSDI.ref) '_' ref_movie '_' roikind '_' refcase '_rej' num2str(reject_on) '_'];
            %
            ncond = length(cond);
            
            %             % ------------------------------------------------------------------
            %             % PLOT MAPS OF AVERAGE MEASURES CONDITION-WISE
            %             % ------------------------------------------------------------------
            %             figure
            %             ci = 0;
            %             for condition =  makeRow(cond)
            %                 % peak in the first row
            %                 ci = ci+1;
            %                 ax(ci) = subplot(2,ncond,ci);
            %                 imagesc(maps.peak(:,:,ci))
            %                 axis image
            %                 colorbar
            %                 %
            %                 set(gca, 'clim', [0 max(maps.peak(:))*0.8])
            %                 condidx = find(VSDI.condition(:,1) ==condition); % get idx from condition
            %                 tempmA = VSDI.condition(condidx(1),4); %get mA from first trial that meet the condition
            %                 title(['c',num2str(condition), '(', num2str(tempmA),'mA)'])
            %
            %                 % wmean in the second row
            %                 ax(ci+ncond) = subplot(2,ncond,ci+ncond);
            %                 imagesc(maps.wmean(:,:,ci))
            %                 axis image
            %                 colorbar
            %                 set(gca, 'clim', [0 max(maps.wmean(:))*0.8])
            %
            %
            %             end
            %
            %             sgtitle([num2str(VSDI.ref), '(',ref_movie,  ')', 'up: peak; down: onsetA'])
            %             localname = ['MAPS'] ;
            %
            %             if savemaps
            %                 saveas(gcf, fullfile(savein, [analysisref savename localname '.jpg']), 'jpg')
            %                 close
            %             end
            
            % ------------------------------------------------------------------
            % PLOT MAPS OF SPATIAL Z-SCORE OF AVERAGE MEASURES CONDITION-WISE
            % ------------------------------------------------------------------
            
            % colormap like jet but changing initial values 
            ccmap = jet;
            % remove darker colors: 
            darkblue = ccmap(4,:);
            ccmap = removerows(ccmap,'ind',[1:8, 1:3]);

            % change first color and interpolate
            ccmap(1,:) = darkblue;
            flag = 3;
            R = linspace(ccmap(1,1), ccmap(flag,1), flag);
            G = linspace(ccmap(1,2), ccmap(flag,2), flag);
            B = linspace(ccmap(1,3), ccmap(flag,3), flag);
            ccmap(1:flag,:) = [R; G; B]'; 
            
            
            cclim = [0 4];
            %         ccmap(1,:) = [0 0 0];
            
            % PLOT: PEAK + WMEAN
            % ------------------------------------------------------------------
            figure
            ci = 0; % counter in case that condition is not 0:3
            for condition =  makeRow(cond)
                % peak in the first row
                ci = ci+1;
                ax(ci) = subplot(2,ncond,ci);
                
                im1 = mapsZ.peak(:,:,ci);
                %                 im1 = interp2(im1, 5, 'nearest');
                %                 im1(~VSDI.crop.preview) = 0;
                
                alphamask = ones(size(im1))*0.6;
                
                %                 alphamask= im1  > 1 ; %VISUALIZATION THRESHOLD
                %                 alphamask =  interp2(alphamask, 5, 'nearest');
                
                back = VSDI.backgr(:,:,VSDI.nonanidx(1));
                back = interp2(back,5);
                %                 imagesc(im1)
                %                 axis image
                %                 set(gca, 'clim',cclim)
                %                 colorbar; colormap(ccmap)
                %
                plot_framesoverlaid2(im1, back, alphamask, 0, ax(ci), cclim, 1, 0, ccmap)
                ax(ci).Visible = 'off';
                % STOPS THE CODE FOR CHECKING ROI CENTERS
                %             if ci ==2
                %                 return
                %             end
                
                condidx = find(VSDI.condition(:,1) ==condition); % get idx from condition
                title(VSDI.conditionlabels{condition+1,2})
                
                % wmean in the second row
                ax(ci+ncond) = subplot(2,ncond,ci+ncond);
                
                im2 = mapsZ.wmean(:,:,ci);
                %                 im2(~VSDI.crop.preview) = 0;
                
                alphamask = ones(size(im2))*0.6;
                %                 alphamask= im2  >0;
                
                back = VSDI.backgr(:,:,VSDI.nonanidx(1));
                back = interp2(back,5);
                plot_framesoverlaid2(im2, back, alphamask, 0, ax(ci+ncond), cclim, 1 , 0, ccmap)
                ax(ci+ncond).Visible = 'off';
                
            end
            
            sgtitle([num2str(VSDI.ref), 'rej' , num2str(reject_on), '(',ref_movie(2:end),  ')', 'Zmap. up: peak; down: onsetA'])
            localname = ['MAPS_Zscore'] ;
            
            if savemaps
                %                 saveas(gcf, fullfile(savein, [analysisref savename localname '.jpg']), 'jpg')
                set(gcf,'units','normalized','outerposition',[0 0 1 1]) %set in total screen
                print(fullfile(savein, [analysisref savename localname '.jpg']),'-r600','-djpeg') % prints it as you see them
                
                close
            end
            
            % PLOT: SLOPE
            % ------------------------------------------------------------------
            figure
            ci = 0; % counter in case that condition is not 0:3
            for condition =  makeRow(cond)
                % peak in the first row
                ci = ci+1;
                ax(ci) = subplot(1,ncond,ci);
                
                im1 = mapsZ.slope(:,:,ci);
                %                 im1 = interp2(im1, 5, 'nearest');
                %                 im1(~VSDI.crop.preview) = 0;
                
                alphamask = ones(size(im1))*0.6;
                
                %                 alphamask= im1  > 1 ; %VISUALIZATION THRESHOLD
                %                 alphamask =  interp2(alphamask, 5, 'nearest');
                
                back = VSDI.backgr(:,:,VSDI.nonanidx(1));
                back = interp2(back,5);
                %                 imagesc(im1)
                %                 axis image
                %                 set(gca, 'clim',cclim)
                %                 colorbar; colormap(ccmap)
                %
                plot_framesoverlaid2(im1, back, alphamask, 0, ax(ci), cclim, 1, 0, ccmap)
                ax(ci).Visible = 'off';
                
            end
            
            sgtitle(['MOT' num2str(VSDI.ref), 'rej' , num2str(reject_on), '(',ref_movie(2:end),  ')', 'Zmap. slope in ms:' num2str(slope.window)])
            localname = ['MAPS_Zscore_slope'] ;
            
            if savemaps
                %                 saveas(gcf, fullfile(savein, [analysisref savename localname '.jpg']), 'jpg')
                set(gcf,'units','normalized','outerposition',[0 0 1 1]) %set in total screen
                print(fullfile(savein, [analysisref savename localname '.jpg']),'-r600','-djpeg') % prints it as you see them
                
                close
            end
            
%             % ------------------------------------------------------------------
%             % PLOT ROIS USED
%             % ------------------------------------------------------------------
%             ax1 = subplot(1,3,1);
%             
%             switch roikind
%                 case 'circle'
%                     centers = VSDI.roi.circle.center(selroi, :) ;
%                     roicirc_preview_multiple(VSDI.crop.preview, centers, VSDI.roi.circle.R, ax1);
%                     
%                 case 'anat'
%                     roi_preview_multiple(VSDI.crop.preview, VSDI.roi.manual_poly(selroi,:), ax1);
%             end
%             
%             localname = [num2str(VSDI.ref), ref_movie,'roipreview',num2str(numel(selroi)) ];
%             
%             if savemaps
%                 saveas(gcf, fullfile(savein, [analysisref localname '.jpg']), 'jpg')
%                 close
%             end
        end
%         
%         
%         %         % PLOT PEAK VALUES
%         %         ax2 = subplot(1,3,2);
%         %         temp = barplot.peak(suji,:,selroi)';
%         %         legend(selroinames)
%         %         title(num2str(suji))
%         %
%         %         localref = [num2str(VSDI.ref) '- cond=' num2str(condition) '-' ref_movie '-' refcase  '-reject' num2str(reject_on) ];
%         %         sgtitle(localref)
%         %
%         %         savename = [ localref   '-' num2str(numel(selroi)) roikind 'ROI - PART2 WAVES'];
%         % %         localname = 'tilewaves'];
%         %         saveas(gcf, fullfile(savein, [savename localname '.jpg']), 'jpg')
%         %         close
%         
%         %
        clear maps mapsZ
        
    end % for nfish
    
    %% -------------------------------------------
    % EXPORT FOR R
    % -------------------------------------------
    if getR
        excelname = fullfile(savein, [analysisref 'Zscore_long_format_forR_group' ref_movie '_' roikind '_' refcase 'rej' num2str(reject_on) '_' num2str(numel(selroi)) 'roi.xls']);
        labels = {'id' 'roi' 'roi n' 'cond' 'peak' 'wmean' 'slope' 'peak_z' 'wmean_z' 'slope_z'};
        
        for col = 1:numel(labels)
            longF{1,col}= labels{col};
        end
        
        % write output (new sheet for each fish
        writecell (longF, excelname, 'sheet',  [roikind ref_movie 'rej' num2str(reject_on)])
        
        writecell (labels, excelname, 'sheet', 'labels')
%         writecell (groupplot_print, excelname, 'sheet','group')
        writecell (params, excelname, 'sheet','param')
        
        clear longF
    end % if getR
    
    blob()
    
    
    %% -------------------------------------------------------
    % T-TEST MATRIX
    % -------------------------------------------------------

    if getTmat
        
        ci = 0;
        for condition =  makeRow(cond)
            ci = ci +1;
            
            for roi1 = makeRow(selroi)
                for roi2 = makeRow(selroi)
                    
                    % Peak
                    peak1=  tmat.peakZ(roi1, :, ci);
                    peak2=  tmat.peakZ(roi2, :, ci);
                    
                    tresults.peak(roi1,roi2, ci) = ttest2(peak1, peak2);
                    
                    % Slope
                    slope1=  tmat.slopeZ(roi1, :, ci);
                    slope2=  tmat.slopeZ(roi2, :, ci);
                    
                    tresults.slope(roi1,roi2, ci) = ttest2(slope1, slope2);                    % Wmean
                   
                    % Wmean
                    wmean1=  tmat.wmeanZ(roi1, :, ci);
                    wmean2=  tmat.wmeanZ(roi2, :, ci);
                    
                    tresults.wmean(roi1,roi2, ci) = ttest2(wmean1, wmean2);                    % Wmean                    
                    
                end
            end
        end
        
    %% plot ttest matrix
    % -------------------------------------------------------
        ci = 0;
        nroi = numel(selroi);
        
        for condition =  makeRow(cond)
            ci = ci +1; %corresponds to condition row (in this experiment)
            
            figure
            
            % PEAK
            subplot(1,3,1)
            imagesc(tresults.peak(:,:,ci))
            title('peak')
            
            axis image
            xticks(1:nroi)
            xticklabels(selroinames); xtickangle(90)
            yticks(1:nroi)
            yticklabels(selroinames); 
            set(gca, 'clim', [0 1])

            % WMEAN
            subplot(1,3,2)
            imagesc(tresults.wmean(:,:,ci))
            title('wmean')
            axis image
            
            xticks(1:nroi)
            xticklabels(selroinames); xtickangle(90)
            yticks(1:nroi)
            yticklabels(selroinames);
            set(gca, 'clim', [0 1])

            % SLOPE
            subplot(1,3,3)
            imagesc(tresults.slope(:,:,ci))
            title('slope')
            axis image
            
            xticks(1:nroi)
            xticklabels(selroinames); xtickangle(90)
            yticks(1:nroi)
            yticklabels(selroinames); ytickangle(90)
            set(gca, 'clim', [0 1])

            
            sgtitle(['Tmap.'  VSDI.conditionlabels{ci,2} '-rej' , num2str(reject_on), '(',ref_movie(2:end),  '). yellow= *; blue = n.s.'])
            
            if saveTmat
                % saveas(gcf, fullfile(savein, [analysisref savename localname '.jpg']), 'jpg')
                set(gcf,'units','normalized','outerposition',[0 0 1 1]) %set in total screen
                print(fullfile(savein, [analysisref 'Tmaps' '.rej' num2str(rej) '_' VSDI.conditionlabels{ci,2}   '.jpg']),'-r600','-djpeg') % prints it as you see them
                
                close
                
            end
        end
    end % if get getTmat
    
                % Created: 09/02/22
% Based on TORus script: /home/tamara/Documents/MATLAB/VSDI/TORus/plot_code/informes_code/03_figure_sketch/groupplot_measures/group_plot_brain_reprogram_4_Zscore_working.m

% Last Update:
% 25/02/22: add Tmaps
% 10/02/22: update with slope code, and flexible number of roi (some fish
% lack dm1 for instance)
