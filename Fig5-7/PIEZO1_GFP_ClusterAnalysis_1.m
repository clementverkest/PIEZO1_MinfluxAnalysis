%% Script for the analysis of MINFLUX signals from DNA-Paint labelled GFP-tagged PIEZOs
    % Dependencies: DBSCAN.m, SetGFPanalysisParameters.m, PlotRawData.m; 
    % Plot raw and filtered data
    
    clear all
    close all

%%  %%%%%%%%%%%%%%%%%%%%%%%% – CHOOSE OPTIONS – %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    opt.plotRaw             = false; % set to 'true' to plot raw data overview 
    opt.plotFiltered        = true; % set to true to plot filtered data overview
    aggregateTraces         = true; % Perform aggregation of X number of loc (group of 3 here) for each trace



%%  %%%%%%%%%%%%%%%%%%%%%%%% – Load data – %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Chose which data to load: 0 = CTL, all data /// 1 = OSMO, all data /// 2
% = CTL Image Fig 1C /// 3 = OSMO image Fig2a /// 4 = CTL image 2a

    DataSource = 1; 

%
if DataSource == 0;
  load("PIEZO1mGL_CTL_RawData.mat");
  myStruct = PIEZO1mGL_CTL_RawData;
  myRange = [1:length(myStruct)];
  elseif DataSource ==1;
  load("PIEZO1mGL_OSMO_RawData.mat");
  myStruct = PIEZO1mGL_OSMO_RawData;
  myRange = [1:length(myStruct)];
  elseif DataSource ==2;
  load("PIEZO1mGL_CTL_RawData.mat");
  myStruct = PIEZO1mGL_CTL_RawData;
  myRange = 13;
  elseif DataSource ==3;
  load("PIEZO1mGL_OSMO_RawData.mat");
  myStruct = PIEZO1mGL_OSMO_RawData;
  myRange = 3;
   elseif DataSource ==4;
  load("PIEZO1mGL_CTL_RawData.mat");
  myStruct = PIEZO1mGL_CTL_RawData;
  myRange = 17;
end


%% create empty result arrays
 LocPerTraceAll = []; numTracesperHour = [];%% create empty result arrays
 

    %% Start loop 
     for k = myRange                  % loop through indicated files

     traces_RAW = myStruct(k).data;
     myfile = myStruct(k).name;

     traces_RAW(:,1:3)=1e9*traces_RAW(:,1:3); % Convert to nm
   

%%  %%%%%%%%%%%%%%%%%%%%%%%% – set analysis parameters – %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [p]=SetGFPanalysisParameters_V2(); % set analysis paramters

%%  %%%%%%%%%%%%%%%%%%%%%%%% – plot raw data – %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if opt.plotRaw
        status = '-raw-';
        PlotRawData(traces_RAW, myfile, status);
    end

%%  %%%%%%%%%%%%%%%%%%%%%%%% – filtering – %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    traces_FILT = traces_RAW;
    traces_FILT(:,3)=0.7*traces_FILT(:,3); % Adjust z coordinates for refractive index mismatch


   % aggregate traces
    if aggregateTraces
    p.loc_per_trace_threshold = 0; traces_Aggregate = []; j=1;
    [w, ~, ~] = unique(traces_FILT(:,4));
    agg = 3;
    for i=1:size(w,1)
    temp = traces_FILT(traces_FILT(:,4)==w(i,1),1:8);
        for k=1:agg:size(temp,1)
            if size(temp,1)<(agg+1)
                continue
            elseif k+(agg-1)>size(temp,1)
                traces_Aggregate(j,1:8) = mean(temp(k:size(temp,1),:),1);
            elseif k+(agg-1)<=size(temp,1)
                traces_Aggregate(j,1:8) = mean(temp(k:k+(agg-1),:));
            end
            j=j+1;
        end
    end
    traces_FILT = traces_Aggregate;
    clear w agg temp;
    end


    % efo, cfr filtering
        traces_FILT = traces_FILT(traces_FILT(:,6) <= p.cfr_threshold, :);
        traces_FILT = traces_FILT(traces_FILT(:,5) <= p.efo_threshold, :); 
        traces_FILT = traces_FILT(traces_FILT(:,5) >= 20000, :);


    % time filtering 
        traces_FILT = traces_FILT(traces_FILT(:,8) <= p.time_threshold, :);  % filter by recording time
    
    % filter by standard deviation
        [uv_tid, ~, id_tid] = unique(traces_FILT(:,4));
        StDevALL=[accumarray(id_tid,traces_FILT(:,1),[],@std) accumarray(id_tid,traces_FILT(:,2),[],@std) accumarray(id_tid,traces_FILT(:,3),[],@std)];
        % create arrays for graph coloring purposes only 
        STDforeachiteration = zeros(size(id_tid,1),3);
        for Tidx = 1:size(id_tid,1)
            STDforeachiteration(Tidx,:)=StDevALL(id_tid(Tidx,1),:);
        end
        StdMAX = max(STDforeachiteration,[],2); % can be used for loc coloring
        traces_FILT = traces_FILT(StdMAX(:,1)<p.stdev_trace_threshold,:);
        StdMAX = StdMAX(StdMAX(:,1)<p.stdev_trace_threshold,:);
   
        
    % filter by localisations per trace 
        [uv_tid, ~, id_tid] = unique(traces_FILT(:,4));     
        n_tid = histcounts(id_tid,size(uv_tid,1)); % How many of each?
        
        % create arrays for graph coloring purposes only 
        LocPerTrace=n_tid';
        Numtraceforeachiteration = id_tid;
        for Tidx = 1:size(id_tid,1)
            Numtraceforeachiteration(Tidx)=LocPerTrace(id_tid(Tidx,1));
        end

        % filtering code
        traces_FILT = traces_FILT(ismember(traces_FILT(:,4), uv_tid(n_tid > p.loc_per_trace_threshold)),:); % Keep ones with more than threshold.

        % filtering code
        % traces_filt = traces_filt(Numtraceforeachiteration(:,1)>p.loc_per_trace_threshold,:);
        Numtraceforeachiteration = Numtraceforeachiteration(Numtraceforeachiteration(:,1)>p.loc_per_trace_threshold,:);
        numTracesperHour(k,1)=size(LocPerTrace,1);
        LocHist = histcounts(LocPerTrace,"BinWidth",1)';
        LocPerTraceAll = cat(1,LocPerTraceAll,LocPerTrace);
        n_tid=n_tid(n_tid(:)>p.loc_per_trace_threshold)';



%%  %%%%% – calculate center of mass for each trace and merge signals from same fluorophore– %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [uv_tid, ~, id_tid] = unique(traces_FILT(:,4)); 
    traces_AVG = [accumarray(id_tid,traces_FILT(:,1),[],@mean) accumarray(id_tid,traces_FILT(:,2),[],@mean) accumarray(id_tid,traces_FILT(:,3),[],@mean)];

    % run DBSCAN to find signals from the same fluorophore 
    ClusterIDs = dbscan(traces_AVG(:,1:3),p.MaxDistProtamers,p.MinNumProtamers);
    
    % add column with cluster ID to main coordinates table
    clust_XYZ = [traces_AVG ClusterIDs n_tid];
    
    % split clust_XYZ into 2 submatrices - one with single signals and one with multiple signals per channel 
    singleloc = clust_XYZ(clust_XYZ(:,4) == -1,:);
    multiloc = clust_XYZ(clust_XYZ(:,4) > -1,:);
    
    % calculate averages of clusters (i.e. signals from one channel)
    [uv, ~, id] = unique(multiloc(:,4));
    MyMeanX = [accumarray(id,multiloc(:,1),[],@mean)];
    MyMeanY = [accumarray(id,multiloc(:,2),[],@mean)];
    MyMeanZ = [accumarray(id,multiloc(:,3),[],@mean)];
    ClID = [accumarray(id,multiloc(:,4),[],@mean)];
    numLocs = [accumarray(id,multiloc(:,5),[],@sum)];
    multiloc = [MyMeanX, MyMeanY, MyMeanZ, ClID, numLocs];
    
    % merge submatrices into single matrix
    traces_AVG = cat(1,singleloc(:,1:5),multiloc(:,1:5));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  %%%%%%%%%%%%%%%%%%%%%%%% – plot filtered data – %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if opt.plotFiltered  
status = '-filtered-';
PlotRawData(traces_FILT,myfile, status)
end



 end % end of main loop
  clearvars -except PIEZO1mGL_CTL_RawData PIEZO1mGL_OSMO_RawData traces_RAW traces_FILT traces_Aggregate traces_AVG StDevALL myfile LocPerTraceAll



