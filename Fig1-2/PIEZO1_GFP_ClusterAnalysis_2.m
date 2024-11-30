%% Script for the analysis of MINFLUX signals from DNA-Paint labelled GFP-tagged PIEZOs
    % Dependencies: DBSCAN.m  
    % Analyze and classify identified PIEZO GFP clusters

clear all
close all


%% %%%%%%%%%%%%%%%%%%%% - Load Data - %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Structures containing selected PIEZOmGL clusters - CTL and OSMO

load 'GFP_CTL_all_selected_clusters.mat'; 
load 'GFP_OSMO_all_selected_clusters.mat';


%% %%%%%%%%%%%%%%% - set options - %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    OpenCutOff = 0.25; % cutoff value for cluster to be open (normalized distance from cluster center)
    SegmentSize = 5; % segment of the cluster use for classification: here first 5th of cluster total height 
    normalize = true; % use absolute or normalized coordinates for clusers  
    PlotRawData1 = false; % plot 3D view of individual cluster ! Might slow down ! USE DataSource == 3 or 4 !
    PlotDonut = true; % plot overlay of all clusters

%% %%%%%%%%% -choose data source - %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DataSource = 2; %% 1 = CTL all clusters, 2 = OSMO all clusters, 3 = CTL exmpales to plot, 4 = OSMO examples to plot

if DataSource == 1
    SelectedPointsTID = SelectedPointsTID_allGFP;
    fns = fieldnames(SelectedPointsTID); 
    myRange = [1:size(fns,1)];
elseif DataSource == 2
    SelectedPointsTID = SelectedPointsTID_OSMO_v2;
    fns = fieldnames(SelectedPointsTID); 
    myRange = [1:size(fns,1)];
elseif DataSource == 3
    SelectedPointsTID = SelectedPointsTID_allGFP;
    fns = fieldnames(SelectedPointsTID); 
    myRange = [14 32 67 75 93 106 113 117 122 124 131 133 141 143 144 153 154 156 161 174 182 187 188 206 222 224 228 242 245 246 249];
elseif DataSource == 4
    SelectedPointsTID = SelectedPointsTID_OSMO_v2;
    fns = fieldnames(SelectedPointsTID); 
    myRange = [4 9 17 24 27 35 36 66 68 71 87 106 111 113];
end




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
            AllHeight = cat(1,AllHeight, ClusterDepth);
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
if DataSource == 1 | DataSource == 2;
        if size(AllOpen,2)>0
        NumClusterOpen = size(unique(AllOpen(:,4)),1)
        end
        
        if size(AllClosed,2)>0
        NumClusterClosed = size(unique(AllClosed(:,4)),1)
        end
        ClusterProportions = [100*NumClusterOpen/(NumClusterClosed+NumClusterOpen), 100*NumClusterClosed/(NumClusterClosed+NumClusterOpen)];
else
end



%% %%%%%%%%%%%%%%%%% - plot data (donut plots) - %%%%%%%%%%%%%%%%%%%%%%%%%%
if DataSource == 1 | DataSource == 2;
if PlotDonut
SymbolSize = 40;
SymbolColor = [0.3 0.3 0.3];
SymbolColorPit = [0 0.35 0.75];
SymbolTransparency = 0.2;
ColorMax = 3.5;
figure('Position', [0 300 1000 170]);
    tiledlayout(1,5);
    set(gcf,'renderer','Painters');
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
else
end


clear SymbolColor SymbolSize SymbolTransparency Z_range4evaluation SegmentSize PlotDonut PlotRawData1 OuterClusterRadius OpenCutOff
clear normalize minYraw minY minXraw minX maxYraw maxY maxX maxXraw k i Hidx H DistMatrix distance DataSource CutoffOpenCluster
clear ColorMax ClusZ ClusterRim ClusterDepth ClusNumSignals ClusIDX_ALL ClusIDX_TOP CenterAraw CenterA 

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
                zlim ([-200 0]);
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
