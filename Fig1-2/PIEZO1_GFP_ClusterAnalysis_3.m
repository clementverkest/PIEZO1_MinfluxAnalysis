%% Script for the analysis of MINFLUX signals from DNA-Paint labelled GFP-tagged PIEZOs
    % Dependencies: DBSCAN.m  
    % Analyze and classify identified PIEZO GFP clusters

clear all
close all


%% %%%%%%%%%%%%%%%%%%%% - Load Data - %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Structures containing selected PIEZOmGL clusters - CTL and OSMO

load 'GFP_CTL_all_selected_clusters.mat'; 
load 'GFP_OSMO_all_selected_clusters.mat';

%% %%%%%%%%% -choose data source - %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for DataSource = 1:2; %% 1 = CTL all clusters, 2 = OSMO all clusters

    if DataSource == 1
    SelectedPointsTID = SelectedPointsTID_allGFP;
    fns = fieldnames(SelectedPointsTID); 
    myRange = [1:size(fns,1)];
    elseif DataSource == 2
    SelectedPointsTID = SelectedPointsTID_OSMO_v2;
    fns = fieldnames(SelectedPointsTID); 
    myRange = [1:size(fns,1)];
    Exclude = [61]; % Cluster removed from depth calc because of outlier traces detected 
    end

%% %%%%%%%%%%%%%%% - set options - %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    OpenCutOff = 0.25; % cutoff value for cluster to be open (normalized distance from cluster center)
    SegmentSize = 5; % segment of the cluster use for classification: here first 5th of cluster total height 
    normalize = false; % use absolute or normalized coordinates for clusers  
    PlotRawData1 = false; % plot 3D view of individual cluster ! Might slow down ! Or adjust loop length !
    PlotDonut = false; % plot overlay of all clusters
    PlotSummaryResults = true % plot summary graphs

    % create empty results sheets and variables
    AllOpen = []; AllClosed = []; AllOpenBOTTOM = [];AllClosedBOTTOM = [];
    AllOpenTOP = []; AllClosedTOP = []; AllMeanDistance = []; AllHeight = []; AllAbsoluteDepth = [];
    AllNumChannels = []; AllNearestNeighbor = []; AllClosedDepth = []; 

    
    % create structure with trace means of clusters
    [IndivClustersRAW_means] = CalculateTraceMean(SelectedPointsTID);

    % create empty array for open/closed cluster assignment
    ClusterCategory = zeros(size(fns,1),1);
  

%%  %%%%%%%%%%%%%%%%%% - loop through data - %%%%%%%%%%%%%%%%%%%%%%%%%%%


for k = myRange % adjust range if desired / myRange = full dataset

    clear ClusXYZRAWinvert;

    ClusXYZRAW = SelectedPointsTID.(fns{k});
    ClusXYZtraceMeans = IndivClustersRAW_means.(fns{k});

    %Calculate cluster depth
    ClusZAbs = sort(ClusXYZtraceMeans(:,3),'descend');
    DepthAbs = mean(ClusZAbs(1:3,1));
    AllAbsoluteDepth = cat(1,AllAbsoluteDepth,DepthAbs);

    % find cluster borders and calculate center coordinates
    minX = min(ClusXYZtraceMeans(:,1),[],1);
    maxX = max(ClusXYZtraceMeans(:,1),[],1);
    minY = min(ClusXYZtraceMeans(:,2),[],1);
    maxY = max(ClusXYZtraceMeans(:,2),[],1);
    CenterA = [minX+(maxX-minX)/2 minY+(maxY-minY)/2 min(ClusXYZtraceMeans(:,3),[],1)];
    ClusXYZtraceMeans=ClusXYZtraceMeans-CenterA;
    OuterClusterRadius = sqrt((maxX-minX)^2 + (maxY-minY)^2)/2;

   % find RAW cluster borders and calculate center coordinates
    minXraw = min(ClusXYZRAW(:,1),[],1);
    maxXraw = max(ClusXYZRAW(:,1),[],1);
    minYraw = min(ClusXYZRAW(:,2),[],1);
    maxYraw = max(ClusXYZRAW(:,2),[],1);
    CenterAraw = [minX+(maxX-minX)/2 minY+(maxY-minY)/2 min(ClusXYZRAW(:,3),[],1)];

    % shift cluster to origin
    ClusXYZRAW(:,1:3) = ClusXYZRAW(:,1:3)-CenterAraw;
    
    % normalize if chosen
    if normalize
        ClusXYZtraceMeans = ClusXYZtraceMeans/OuterClusterRadius; %% new
        CutoffOpenCluster = OpenCutOff;
    else
        CutoffOpenCluster = OuterClusterRadius*OpenCutOff;
    end

        maxZ = max(ClusXYZtraceMeans(:,3),[],1);
        ClusZ = sort(ClusXYZtraceMeans(:,3),'descend');
        ClusterDepth = mean(ClusZ(1:3,1));
        ClusNumSignals = size(ClusZ,1);
        ClusterRim = mean(ClusZ((ClusNumSignals-2):ClusNumSignals,1));
        ClusterDepth = mean(ClusZ(1:3,1))-ClusterRim;


        % split cluster into segments
        Z_range4evaluation = maxZ/SegmentSize;
        ClusXYZ_TOP = ClusXYZtraceMeans(ClusXYZtraceMeans(:,3)<Z_range4evaluation,:);
        ClusXYZ_BOTTOM = ClusXYZtraceMeans(ClusXYZtraceMeans(:,3)>maxZ-Z_range4evaluation,:);

        % add index column for later identification and debugging
        ClusIDX_TOP = zeros(size(ClusXYZ_TOP,1),1)+k;
        ClusXYZ_TOP = [ClusXYZ_TOP ClusIDX_TOP];
    
        ClusIDX_ALL = zeros(size(ClusXYZtraceMeans,1),1)+k;
        ClusXYZtraceMeans = [ClusXYZtraceMeans ClusIDX_ALL];
        
        % calculate distance matrix
        distance = sqrt(ClusXYZ_TOP(:,1).^2+ClusXYZ_TOP(:,2).^2);
        distance2 = mean(distance,1);
        distAll = sqrt(ClusXYZtraceMeans(:,1).^2+ClusXYZtraceMeans(:,2).^2);
        distVSz = [distAll ClusXYZtraceMeans(:,3)];
 
   %% classify clusters %%
        % pit-shaped cluster
        if min(distance)>CutoffOpenCluster && size(ClusXYZ_TOP,1)>0
            AllOpen = cat(1, AllOpen, ClusXYZtraceMeans);
            AllOpenTOP = cat(1, AllOpenTOP, ClusXYZ_TOP);
            AllOpenBOTTOM = cat(1, AllOpenBOTTOM, ClusXYZ_BOTTOM);
            AllMeanDistance = cat(1,AllMeanDistance,distance2);
            AllHeight= cat(1,AllHeight, ClusterDepth); 
               
            ClusterCategory(k,1)=1;
          
            % analyse cluster
            NumChannels = size(ClusXYZtraceMeans,1);
            DistMatrix=squareform(pdist(ClusXYZtraceMeans(:,1:3)));
            [H, Hidx] = sort(DistMatrix,1);
            NearestNeighbor = H(2,:)';
            if normalize
                NearestNeighbor = NearestNeighbor*OuterClusterRadius; %% new
            end
                    
            % append results from this iteration to result array
            AllNumChannels = cat(1,AllNumChannels,NumChannels);
            AllNearestNeighbor = cat(1,AllNearestNeighbor,NearestNeighbor);

            % plot cluster as scatter3 - optional    
            if PlotRawData1
                PlotCluster3D(ClusXYZRAW, k)
            end

        % spherical cluster    
        else
            AllClosed = cat(1, AllClosed, ClusXYZtraceMeans);
            AllClosedTOP = cat(1, AllClosedTOP, ClusXYZ_TOP);
            AllClosedBOTTOM = cat(1, AllClosedBOTTOM, ClusXYZ_BOTTOM);
            AllClosedDepth = cat(1,AllClosedDepth, ClusterDepth);
            % PlotCluster3D(ClusXYZRAW, k)
        end


end % end of loop through cluster list




%% %%%%%%%%%%%%%%%%%% - count clusters per category - %%%%%%%%%%%%%%%%%%%%%%

        if size(AllOpen,2)>0
        NumClusterOpen = size(unique(AllOpen(:,4)),1)
        end
        
        if size(AllClosed,2)>0
        NumClusterClosed = size(unique(AllClosed(:,4)),1)
        end
        ClusterProportions = [100*NumClusterOpen/(NumClusterClosed+NumClusterOpen), 100*NumClusterClosed/(NumClusterClosed+NumClusterOpen)];
 




%% %%%%%%%%%%%%%%%%% - plot data (donut plots) - %%%%%%%%%%%%%%%%%%%%%%%%%%

if PlotDonut
SymbolSize = 40;
SymbolColor = [0.3 0.3 0.3];
    if DataSource == 1  ;
    SymbolColorPit = [0 0.35 0.75];
    SymbolTransparency = 0.4;
    condition = 'CTL';
    elseif DataSource == 2  ;
    SymbolColorPit = [0 0.25 0.60];
    SymbolTransparency = 1;
    condition = 'OSMO';
    end
ColorMax = 3.5;
figure('Position', [0 300 1000 170]);
    tiledlayout(1,5);
    set(gcf,'renderer','Painters');
annotation('textbox', [0.02 0.5 0.2 0.1], 'String', condition, ...
    'FitBoxToText', 'on', ...
    'EdgeColor', 'none', ...
    'FontSize', 12);
nexttile
    scatter3(AllOpenTOP(:,1), AllOpenTOP(:,2),AllOpenTOP(:,3),SymbolSize,SymbolColorPit,'filled','MarkerFaceAlpha',SymbolTransparency);
    colormap parula; axis equal; axis tight
    clim([0 ColorMax]);
    view(0,90)
    if normalize
        xlim([-1 1]);ylim([-1 1]);
    end
    title('pit-shaped cluster (top)');
nexttile
    scatter3(AllOpenBOTTOM(:,1), AllOpenBOTTOM(:,2),AllOpenBOTTOM(:,3),SymbolSize,SymbolColorPit,'filled', 'MarkerFaceAlpha',SymbolTransparency);
    colormap parula; axis equal; axis tight
    clim([0 ColorMax]);
    view(0,90)
    if normalize
        xlim([-1 1]);ylim([-1 1]);
    end
    title('pit-shaped cluster (bottom)');
    if size(AllClosed,2)>0
nexttile
    scatter3(AllClosedTOP(:,1), AllClosedTOP(:,2),AllClosedTOP(:,3),SymbolSize,SymbolColor,'filled', 'MarkerFaceAlpha',SymbolTransparency);
    colormap parula; axis equal; axis tight
    clim([0 ColorMax]);
    view(0,90)
    if normalize
        xlim([-1 1]);ylim([-1 1]);
    end
    title('spherical cluster (top)');

nexttile
    scatter3(AllClosedBOTTOM(:,1), AllClosedBOTTOM(:,2),AllClosedBOTTOM(:,3),SymbolSize, SymbolColor,'filled', 'MarkerFaceAlpha',SymbolTransparency);
    colormap parula; axis equal; axis tight
    clim([0 ColorMax]);
    view(0,90)
    if normalize
        xlim([-1 1]);ylim([-1 1]);
    end
    title('spherical cluster (bottom)');
    end
nexttile
% barX = 1;
ba = bar(1,ClusterProportions,'stacked');
set(ba, 'FaceColor', 'Flat')
set(ba, 'FaceAlpha', 0.4)
title('Proportion (%)');
ba(1).CData = [0 0.35 0.75];  
ba(2).CData = [0.3 0.3 0.3];  
end



clear SymbolColor SymbolSize SymbolTransparency Z_range4evaluation SegmentSize PlotDonut PlotRawData1 OuterClusterRadius OpenCutOff
clear normalize minYraw minY minXraw minX maxYraw maxY maxX maxXraw k i Hidx H DistMatrix distance  CutoffOpenCluster
clear ColorMax ClusZ ClusterRim ClusterDepth ClusNumSignals ClusIDX_ALL ClusIDX_TOP CenterAraw CenterA 



%% %%%%%%%%%%%%%%%%%% - Append and plot results - %%%%%%%%%%%%%%%%%%%%%%%%%%

if PlotSummaryResults

    % result.ClusterDepth = [AllHeight];
    % result.ChannelsPerCluster = [AllNumChannels];
    % result.Proportion =[ClusterProportions];
    % result.Radius = [AllMeanDistance];
    
    if DataSource == 1  ;
        result.ClusterDepth_CTL = [AllHeight];
    result.ChannelsPerCluster_CTL = [AllNumChannels];
    result.Proportion_CTL =[ClusterProportions];
    result.Radius_CTL = [AllMeanDistance];
        result.color_CTL = [0 0.35 0.75];
        result.alpha_CTL = 0.4;
    elseif DataSource == 2  ;
                if DataSource == 2;
                AllHeight(Exclude) = [];
                end
         result.ClusterDepth_OSMO = [AllHeight];
    result.ChannelsPerCluster_OSMO = [AllNumChannels];
    result.Proportion_OSMO =[ClusterProportions];
    result.Radius_OSMO = [AllMeanDistance];
        result.color_OSMO = [0 0.25 0.60];
        result.alpha_OSMO = 1;
    end
end



end

[result] = PlotClusterAnalysisResult(result);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% – END – %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% – subroutines – %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% calculate trace means
function [IndivClustersRAW_means] = CalculateTraceMean(IndivClustersRAW)
    fns = fieldnames(IndivClustersRAW);

    for k=1:size(fns,1)
    % read cluster data from structure
    ClusXYZRAW = IndivClustersRAW.(fns{k});
    ClusXYZ = ClusXYZRAW(:,1:4);

    % calculate trace center
    [uv_tid, ~, id_tid] = unique(ClusXYZ(:,4)); 
    ClusterTraces_AVG = [accumarray(id_tid,ClusXYZ(:,1),[],@mean) accumarray(id_tid,ClusXYZ(:,2),[],@mean) accumarray(id_tid,ClusXYZ(:,3),[],@mean)];

    % run DBSCAN to find signals from the same fluorophore – %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ClusterIDs = dbscan(ClusterTraces_AVG(:,1:3),25,2);
    
    % add column with cluster ID to main coordinates table
    ClusXYZ = [ClusterTraces_AVG ClusterIDs];
    
    % split ClusXYZ into 2 submatrices - one with single signals and one with multiple signals per channel 
    singleloc = ClusXYZ(ClusXYZ(:,4) == -1,:);
    multiloc = ClusXYZ(ClusXYZ(:,4) > -1,:);
    
    % calculate averages of multiloc clusters (i.e. signals from one channel)
    [uv, ~, id] = unique(multiloc(:,4));
    MyMeanX = [accumarray(id,multiloc(:,1),[],@mean)];
    MyMeanY = [accumarray(id,multiloc(:,2),[],@mean)];
    MyMeanZ = [accumarray(id,multiloc(:,3),[],@mean)];
    ClID = [accumarray(id,multiloc(:,4),[],@mean)];
    multiloc = [MyMeanX, MyMeanY, MyMeanZ, ClID];
    
    % merge submatrices into single matrix
    ClusterTraces_AVG = cat(1,singleloc(:,1:4),multiloc(:,1:4));
    IndivClustersRAW_means.(fns{k}) = ClusterTraces_AVG(:,1:3);

    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot cluster as 3D scatter
function PlotCluster3D(ClusXYZRAW, k)
            ClusXYZRAWinvert(:,3)=-1*ClusXYZRAW(:,3);
            figure;
            set(gcf,'renderer','Painters');
            tiledlayout(1,2);
            nexttile
                scatter3(ClusXYZRAW(:,1),ClusXYZRAW(:,2),ClusXYZRAWinvert(:,3),10,ClusXYZRAWinvert(:,3),'filled', 'MarkerFaceAlpha', 0.6);
                title(string(k)); colormap parula; caxis([-250 0]); axis equal;
                % xlim([-200 +200]); ylim([-200 +200]);
                zlim ([-320 10]);
                xlabel('x (nm)'); ylabel('y (nm)'); zlabel('z (nm)');
                ax = gca;
                ax.XAxis.FontSize = 14; ax.YAxis.FontSize = 14; ax.ZAxis.FontSize = 14;
                ax.GridColor = 'white';
                ax.GridAlpha = 0.4;
                ax.GridLineWidth = 0.5;
                ax.Color = 'black';
                hold on;
            nexttile
                scatter3(ClusXYZRAW(:,1),ClusXYZRAW(:,2),ClusXYZRAWinvert(:,3),10,ClusXYZRAWinvert(:,3),'filled', 'MarkerFaceAlpha', 0.6);
                title(string(k)); colormap(flipud(parula)); colorbar;           
                axis equal; xlabel('x (nm)'); ylabel('y (nm)'); zlabel('z (nm)');
                ax = gca;
                ax.XAxis.FontSize = 14; ax.YAxis.FontSize = 14; ax.ZAxis.FontSize = 14;
                ax.GridColor = 'White';
                ax.GridAlpha = 0.4;
                ax.GridLineWidth = 0.5;
                ax.Color = 'black';
                view(0, 90); 
                c = colorbar; c.Label.String = 'z (nm)'; c.Label.FontSize = 14; caxis([-250 0]);
                xlim([-200 +200]); ylim([-200 +200]);
end
