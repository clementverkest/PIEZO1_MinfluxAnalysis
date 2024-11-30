%% Script for the analysis of MINFLUX signals from DNA-Paint labelled ALFA-tagged PIEZO1
    % Dependencies: DBSCAN.m; 
    % SetALFAanalysisParameters.m, PlotRawData.m; PIEZO1Superparticle;
    % PlotTrimerAnalysisResult.m; CalcTrimerAngle1.m
    
    clear all;
    close all;
    Exclude = [];
   

    %%  %%%%%%%%%%%%%%%%%%%%%%%% – CHOOSE OPTIONS – %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    plotRaw             = false; % set to 'true' for raw data summary plot
    filter              = true; % set to 'true' if data should be filtered 
    plotFiltered        = false; % plot filtered data summary plot
    plotTraceAVGs       = false;
    aggregateTraces     = false; % aggregate group of loc for every trace
    TrimTrace           = true; % trim first 2 loc of every trace
    PlotResult          = true; % plot final graphs

    % create empty result arrays    
    result.Trimers = [];
    result.Angles = [];
    result.InterBlades =[];
    result.Counts = [];
    result.NumTraces = [];
    result.InterBladesWithNN = [];
    result.TrimerCount=0;
    result.TrimerCoord = [];
    result.LocPerTrace = [];
    result.AllStdTrimer = [];
    numTracesperHour = [];
    AllStdTrimer = [];


%% %%%%%%%%%%%%%%%%%%%% – Load data – %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% choose which data to analyse: 0 = soma, 1 = neurite, 2 = cytochalasin-D
   
    DataSource = 0;

    if DataSource == 0
        load("PIEZO1_ALFAmGL_SOMA_RawData.mat");
         myStruct = PIEZO1_ALFAmGL_SOMA_RawData;
    result.color = [0 0.5 0];
    elseif DataSource == 1
        load("PIEZO1_ALFAmGL_NEURITE_RawData.mat");
         myStruct = PIEZO1_ALFAmGL_NEURITE_RawData;
        result.color = [0.6 0.4 0.8];
    Exclude = [14 16]; %% trimers excluded because the signals are in the soma and not in neurite
    else 
    load("PIEZO1_ALFAmGL_CYTOD_RawData.mat");
             myStruct = PIEZO1_ALFAmGL_CYTOD_RawData;
             result.color = [1 0 0];
        Exclude = [2 3 16 28 36]'; %% trimers excluded because the signals are not likely in the plasma membrane
    end



%% %%%%%%%%%%%%%%%%%%%% – Start loop – %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for k = 1:length(myStruct)      % loop through all files
  % for k = 5                  % loop through indicated files
     traces_RAW = myStruct(k).data;
     myfile = myStruct(k).name;
    
    
    traces_RAW(:,3) = 0.7*traces_RAW(:,3); % adjust for refractive index mismatch
    traces_RAW(:,1:3)=1e9*traces_RAW(:,1:3); % Convert to nm


 %%  %%%%%%%%%%%%%%%%%%%%%%%% – set analysis parameters – %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [p]=SetALFAanalysisParameters(); % set analysis paramters
    p.z_Upper_threshold = 1e9*myStruct(k).UpperZ;  % Z-filt threshold
    p.z_Lower_threshold = 1e9*myStruct(k).LowerZ;  % Z-filt threshold     
      
 %%  %%%%%%%%%%%%%%%%%%%%%%%% – create empty results arrays for iteration – %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Trimers = [];
    TrimWithNN = [];
    InterBlade = [];
    Angles =[];
    TrimerCoord = [];
     tempALL = [];     
      
   %%  %%%%%%%%%%%%%%%%%%%%%%%% – plot raw data – %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if plotRaw
        status = 'raw'
        PlotRawData(traces_RAW,myfile, status);
    end   
  
  %%  %%%%%%%%%%%%%%%%%%%%%%%% – filtering – %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
  traces_FILT = traces_RAW;

   % aggregate traces
    if aggregateTraces
    p.loc_per_trace_threshold = 0; traces_Aggregate = []; j=1;
    [a, ~, ~] = unique(traces_FILT(:,4));
    agg = 3;
    for i=1:size(a,1)
    temp = traces_FILT(traces_FILT(:,4)==a(i,1),1:8);
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
    clear a agg temp;
    end

% trim traces
    if TrimTrace
        j=1;
        [a, ~, ~] = unique(traces_FILT(:,4));
        for i=1:size(a,1)
        temp = traces_FILT(traces_FILT(:,4)==a(i,1),1:8);
        if size(temp,1)>1
        temp(1:2,:)=[];
        tempALL = cat(1,tempALL,temp);
        end
        end
        traces_FILT = tempALL;
        clear j temp tempALL;
    end

 
if filter
        % Z-filtering
        traces_FILT(traces_FILT(:, 3) < p.z_Lower_threshold, :)= []; % filter by Z upper limit 160,:) = []; 
        traces_FILT(traces_FILT(:, 3) > p.z_Upper_threshold, :)= []; % filter by Z lower limit
        
 
        

        % efo, cfr filtering

        traces_FILT = traces_FILT(traces_FILT(:,5) <= p.efo_threshold, :);                          % filter by efo
        traces_FILT = traces_FILT(traces_FILT(:,6) <= p.cfr_threshold, :);                          % filter by cfr
            

        % time filtering 
        traces_FILT = traces_FILT(traces_FILT(:,8) <= p.time_threshold, :);                         % filter by recording time
        

        % filter by standard deviation
        [~, ~, id_tid] = unique(traces_FILT(:,4));
        StDevALL=[accumarray(id_tid,traces_FILT(:,1),[],@std) accumarray(id_tid,traces_FILT(:,2),[],@std) accumarray(id_tid,traces_FILT(:,3),[],@std)];
        % create arrays for graph coloring purposes only 
        STDforeachiteration = zeros(size(id_tid,1),3);
        for Tidx = 1:size(id_tid,1)
            STDforeachiteration(Tidx,:)=StDevALL(id_tid(Tidx,1),:);
        end

        StdMIN = mean(STDforeachiteration,2);
        traces_FILT = traces_FILT(StdMIN(:,1)<p.stdev_trace_threshold,:);




        % filter by localisations per trace 
        [uv_tid, ~, id_tid] = unique(traces_FILT(:,4));         % Unique elements and locations in third column  
        n_tid = histcounts(id_tid,size(uv_tid,1)); % How many of each?
        % create arrays for graph coloring purposes only 
        LocPerTrace=n_tid';
        Numtraceforeachiteration = id_tid;
        for Tidx = 1:size(id_tid,1)
            Numtraceforeachiteration(Tidx)=LocPerTrace(id_tid(Tidx,1));
        end
        

        traces_FILT = traces_FILT(ismember(traces_FILT(:,4), uv_tid(n_tid > p.loc_per_trace_threshold)),:); % Keep ones with more than threshold.
        STDforeachiteration = STDforeachiteration(Numtraceforeachiteration(:,1)>p.loc_per_trace_threshold,:);
        Numtraceforeachiteration = Numtraceforeachiteration(Numtraceforeachiteration(:,1)>p.loc_per_trace_threshold,:);
        
        numTracesperHour(k,1)=size(LocPerTrace,1);

        result.LocPerTrace = cat(1,result.LocPerTrace,LocPerTrace);
        n_tid=n_tid(n_tid(:)>p.loc_per_trace_threshold)';
                 


    else % the following operations must be performed when no filtering is applied in order for the script to work
        p.loc_per_trace_threshold=0;
        [uv_tid, ~, id_tid] = unique(traces_FILT(:,4));
        StDevALL=[accumarray(id_tid,traces_FILT(:,1),[],@std) accumarray(id_tid,traces_FILT(:,2),[],@std) accumarray(id_tid,traces_FILT(:,3),[],@std)];
        n_tid = histcounts(id_tid,"BinWidth",1);                                                            % How many of each?
        traces_FILT = traces_FILT(ismember(traces_FILT(:,4), uv_tid(n_tid > p.loc_per_trace_threshold)),:); % Keep ones with more than threshold.
        n_tid=n_tid(n_tid(:)>p.loc_per_trace_threshold)';
        myNumLoc = size(n_tid);
        result.NumTraces = cat(1,result.NumTraces, myNumLoc);
        result.Counts = cat(1,result.Counts, n_tid);
    end


    % clean up
    [~, ~, Std_tid] = unique(traces_FILT(:,4));
    StDevALL=[accumarray(Std_tid,traces_FILT(:,1),[],@std) accumarray(Std_tid,traces_FILT(:,2),[],@std) accumarray(Std_tid,traces_FILT(:,3),[],@std)];
    clear StdMIN LocPerTrace


     %% plot filered data - optional
    if plotFiltered
        status = 'filtered'
        PlotRawData(traces_FILT,myfile, status);
    end


    %%  %%%%%%%%%%%%%%%%%%%%%%%% – calculate center of mass for each trace – %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [~, ~, id_tid] = unique(traces_FILT(:,4)); 
    CenterX = [accumarray(id_tid,traces_FILT(:,1),[],@mean)];
    CenterY = [accumarray(id_tid,traces_FILT(:,2),[],@mean)];
    CenterZ = [accumarray(id_tid,traces_FILT(:,3),[],@mean)];
    OrigTID = [accumarray(id_tid,traces_FILT(:,4),[],@mean)];
    traces_AVG = [CenterX CenterY CenterZ];

    % clean up
    clear CenterX CenterY CenterZ;


%%  %%%%%%%%%%%%%%%%%%%%%%%% – run DBSCAN to find signals from the same fluorophore – %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ClusterIDs = dbscan(traces_AVG(:,1:3),p.MaxDistProtamers,p.MinNumProtamers);
    
    % add column with cluster ID to main coordinates table
    traces_AVG_noMerge = [traces_AVG ClusterIDs n_tid OrigTID];
    
    % split clust_XYZ into 2 submatrices - one with single signals and one with multiple signals per channel 
    singleloc = traces_AVG_noMerge(traces_AVG_noMerge(:,4) == -1,:);
    multiloc = traces_AVG_noMerge(traces_AVG_noMerge(:,4) > -1,:);
    
    % calculate averages of clusters (i.e. signals from one channel)
    [~, ~, id] = unique(multiloc(:,4));
    MyMeanX = [accumarray(id,multiloc(:,1),[],@mean)];
    MyMeanY = [accumarray(id,multiloc(:,2),[],@mean)];
    MyMeanZ = [accumarray(id,multiloc(:,3),[],@mean)];
    OrigTIDmulti = [accumarray(id,multiloc(:,6),[],@min)];
    ClID = [accumarray(id,multiloc(:,4),[],@mean)];
    numLocs = [accumarray(id,multiloc(:,5),[],@sum)];
    multiloc = [MyMeanX, MyMeanY, MyMeanZ, ClID, numLocs,OrigTIDmulti];
    
    % merge submatrices into single matrix
    traces_AVG = cat(1,singleloc(:,1:6),multiloc(:,1:6));

    % clean up
    traces_AVG_noMerge = cat(2,traces_AVG_noMerge,StDevALL);
    clear MyMeanX MyMeanY MyMeanZ singleloc multiloc ClusterIDs ClID StDevALL numLocs


    %%  %%%%%%%%%%%%%%%%%%%%%%%% – find trimers – %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % create distance matrix
    DistMatrix=squareform(pdist(traces_AVG(:,1:3)));
    [H, Hidx] = sort(DistMatrix,1);
    NearestNeighbor = H(2,:)';

    % find trimers
    for i=1:length(H)
        % check if third nearest neighbors are further away than p.MinTrimerDist
        if(H(4,i)>p.MinTrimerDist && H(4,Hidx(2,i))>p.MinTrimerDist && H(4,Hidx(3,i))>p.MinTrimerDist)
            % check if all three triangle sides are shorter than p.MaxInterBladeDist
            if(H(2,i)<p.MaxInterBladeDist && H(3,i)<p.MaxInterBladeDist && DistMatrix(Hidx(2,i),Hidx(3,i))<p.MaxInterBladeDist)    
                 Trimers = cat(1,Trimers, sort(Hidx(1:3,i))');
                 TrimWithNN = cat(1, TrimWithNN, [H(2,i) H(3,i) DistMatrix(Hidx(2,i),Hidx(3,i)) mean(H(4:6,i))]);
                 InterBlade = cat(1, InterBlade, [H(2,i) H(3,i) DistMatrix(Hidx(2,i),Hidx(3,i))]);
                 TrimerCoord = cat(1,TrimerCoord, [traces_AVG(Hidx(1,i),1:3) traces_AVG(Hidx(2,i),1:3) traces_AVG(Hidx(3,i),1:3)]);
            end
        end
    end

    % evaluate trimers
    if size(Trimers,1)>0
        % find and remove duplicates
        [~,ia,~] = unique(Trimers(:,1:3),'rows'); % find duplicates
        Trimers = Trimers(ia,:); % remove duplicates
        TrimWithNN = TrimWithNN(ia,:);
        InterBlade = InterBlade(ia,:); % remove duplicates
        TrimerCoord = TrimerCoord(ia,:);

        % calculate angles
        for j = 1:size(Trimers,1)
            PointA = traces_AVG(Trimers(j,1),1:3);
            PointB = traces_AVG(Trimers(j,2),1:3);
            PointC = traces_AVG(Trimers(j,3),1:3);
            [TrimerAngles] = CalcTrimerAngle1(PointA, PointB, PointC);           
            Angles = cat(1,Angles, TrimerAngles);
        end
        MaxAngleTrimer = max(Angles,[],2);

        % filter by angle
            maxAngle = 120;
            Trimers=Trimers(MaxAngleTrimer(:,1)<maxAngle,:);        
            InterBlade=InterBlade(MaxAngleTrimer(:,1)<maxAngle,:);
            TrimerCoord=TrimerCoord(MaxAngleTrimer(:,1)<maxAngle,:);
            Angles=Angles(MaxAngleTrimer(:,1)<maxAngle,:);
            TrimWithNN=TrimWithNN(MaxAngleTrimer(:,1)<maxAngle,:); 
        
            % get STD of trimers
            for n = 1:size(Trimers,1)
                idA = traces_AVG(Trimers(n,1),6);
                idB = traces_AVG(Trimers(n,2),6);
                idC = traces_AVG(Trimers(n,3),6);
                SizeA = size(traces_FILT(traces_FILT(:,4)==idA,1),1);
                SizeB = size(traces_FILT(traces_FILT(:,4)==idB,1),1);
                SizeC = size(traces_FILT(traces_FILT(:,4)==idC,1),1);
                StdA = [std(traces_FILT(traces_FILT(:,4)==idA,1)) std(traces_FILT(traces_FILT(:,4)==idA,2)) std(traces_FILT(traces_FILT(:,4)==idA,3)) SizeA];
                StdB = [std(traces_FILT(traces_FILT(:,4)==idB,1)) std(traces_FILT(traces_FILT(:,4)==idB,2)) std(traces_FILT(traces_FILT(:,4)==idB,3)) SizeB];
                StdC = [std(traces_FILT(traces_FILT(:,4)==idC,1)) std(traces_FILT(traces_FILT(:,4)==idC,2)) std(traces_FILT(traces_FILT(:,4)==idC,3)) SizeC];
                StdTrimer = [StdA; StdB; StdC];
                AllStdTrimer = cat(1,AllStdTrimer, StdTrimer);
                % extract raw coordinates
                CoordA = traces_FILT(traces_FILT(:,4)==idA,1:3);
                CoordB = traces_FILT(traces_FILT(:,4)==idB,1:3);
                CoordC = traces_FILT(traces_FILT(:,4)==idC,1:3);
                IndivTrimersRAW.(strcat('Trimer',string(k),'sub',string(n))).A = CoordA;
                IndivTrimersRAW.(strcat('Trimer',string(k),'sub',string(n))).B = CoordB;
                IndivTrimersRAW.(strcat('Trimer',string(k),'sub',string(n))).C = CoordC;
            end

            % clean up
            clear idA idB idC StdA StdB StdC StdTrimer

        % write results
            result.Trimers = cat(1,result.Trimers, Trimers);
            result.InterBlades = cat(1, result.InterBlades, InterBlade);
            result.TrimerCoord = cat(1,result.TrimerCoord,TrimerCoord);
            result.Angles = cat(1, result.Angles, Angles);
            result.InterBladesWithNN = cat(1, result.InterBladesWithNN,TrimWithNN);
    end

    % clean up
    clear TrimWithNN PointA PointB PointC InterBlade TrimerCoord Angles MaxAngleTrimer maxAngle TrimerAngles

%%  %%%%%%%%%%%%%%%%%%%%%%%% – plot processed data – %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if plotTraceAVGs
        figure;
        scatter3(traces_AVG(:,1), traces_AVG(:,2), traces_AVG(:,3),20, [0.8 0.8 0.8], 'filled');
        colormap parula; axis equal; title(myfile); hold on;
    
        % add trimers to plot
        for k=1:size(Trimers,1)
            result.TrimerCount=result.TrimerCount+1;
            scatter3(traces_AVG(Trimers(k,1),1), traces_AVG(Trimers(k,1),2), traces_AVG(Trimers(k,1),3),50, result.color, 'filled')
            scatter3(traces_AVG(Trimers(k,2),1), traces_AVG(Trimers(k,2),2), traces_AVG(Trimers(k,2),3),50, result.color, 'filled')
            scatter3(traces_AVG(Trimers(k,3),1), traces_AVG(Trimers(k,3),2), traces_AVG(Trimers(k,3),3),50, result.color, 'filled')
            text(traces_AVG(Trimers(k,1),1), traces_AVG(Trimers(k,1),2), traces_AVG(Trimers(k,1),3),sprintfc(' %d',result.TrimerCount))
        end
    end


  % main loop end here
  end
  



%% - post processing - %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% combine filtered data into single array
traces_FILT = cat(2,traces_FILT,STDforeachiteration, Numtraceforeachiteration);

% create array with distances of nearest neighbors in trimers
result.NearestNeighborInTrimer = sort(result.InterBlades,2);
result.NearestNeighborInTrimer = cat(1,result.NearestNeighborInTrimer(:,1),result.NearestNeighborInTrimer(:,2),result.NearestNeighborInTrimer(:,1));

if PlotResult;

    
% plot trimer analysis results
result.TrimerCoord(Exclude,:) = [];
result.InterBlades(Exclude,:) = [];
result.Angles(Exclude,:) = [];
result.InterBladesWithNN(Exclude,:) = [];
result.AllStdTrimer = [AllStdTrimer];
[result] = PlotTrimerAnalysisResult(result, DataSource);

end




%% – final clean up – %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear uv_tid uv n_tid
clear plotRaw  filter...
    plotFiltered FindClusters...
    FitSurface plotTraceAVGs MakeConfocalOverlay
clear i j k K ia id id_tid H Hidx FileIdx
clear myfile UpperZ LowerZ include result.TrimerCount Tidx DistMatrix Trimers data
clear Std_tid STDforeachiteration Numtraceforeachiteration NearestNeighbor aggregateTraces TrimTrace SizeA SizeB SizeC 
