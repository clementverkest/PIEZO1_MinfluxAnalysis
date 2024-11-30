function [result] = PlotClusterAnalysisResult(result)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot analysis results

fig = figure('Position',[0, 50 1100 300]);
tiledlayout(1,5);

%% plot mean Cluster depth
nexttile
 
    bar(1, mean(result.ClusterDepth_CTL(:,1)),'FaceColor',result.color_CTL, 'FaceAlpha', result.alpha_CTL);
    hold on;
    swarmchart(ones(size(result.ClusterDepth_CTL,1)),result.ClusterDepth_CTL(:,1),'MarkerFaceColor',[0.7 0.7 0.7], 'MarkerEdgeColor','k');
    ylabel('CTL Cluster depth (nm)');
    title('Pit-shape Cluster depth');

    bar(2, mean(result.ClusterDepth_OSMO(:,1)),'FaceColor',result.color_OSMO, 'FaceAlpha', result.alpha_OSMO);
    swarmchart(2*ones(size(result.ClusterDepth_OSMO,1)),result.ClusterDepth_OSMO(:,1),'MarkerFaceColor',[0.7 0.7 0.7], 'MarkerEdgeColor','k');
    
xticks([1 2]);
xticklabels({'CTL', 'OSMO'});


    hold off
%% plot Number of PIEZO channels per cluster
    nexttile
     bar(1, mean(result.ChannelsPerCluster_CTL(:,1)),'FaceColor',result.color_CTL, 'FaceAlpha', result.alpha_CTL);
     title('# Channels / pit-shape clust.');
    hold on;
    swarmchart(ones(size(result.ChannelsPerCluster_CTL,1)),result.ChannelsPerCluster_CTL(:,1),'MarkerFaceColor',[0.7 0.7 0.7], 'MarkerEdgeColor','k');
    ylabel('Channels/cluster');

bar(2, mean(result.ChannelsPerCluster_OSMO(:,1)),'FaceColor',result.color_OSMO, 'FaceAlpha', result.alpha_OSMO);
swarmchart(2*ones(size(result.ChannelsPerCluster_OSMO,1)),result.ChannelsPerCluster_OSMO(:,1),'MarkerFaceColor',[0.7 0.7 0.7], 'MarkerEdgeColor','k');
    
xticks([1 2]);
xticklabels({'CTL', 'OSMO'});

hold off
%% plot proportion
nexttile
ba_ctl = bar(1,result.Proportion_CTL,'stacked');
hold on
set(ba_ctl, 'FaceColor', 'Flat')
set(ba_ctl, 'FaceAlpha', result.alpha_CTL)
title('Proportion (%)');
ba_ctl(1).CData = result.color_CTL;  
ba_ctl(2).CData = [0.3 0.3 0.3];  

ba_osmo = bar(2,result.Proportion_OSMO,'stacked');
set(ba_osmo, 'FaceColor', 'Flat')
set(ba_osmo, 'FaceAlpha', result.alpha_OSMO)
title('Proportion (%)');
ba_osmo(1).CData = result.color_OSMO;  
ba_osmo(2).CData = [0.3 0.3 0.3]; 
xticks([1 2]);
xticklabels({'CTL', 'OSMO'});

 hold off
%% plot Cluster radius histogram
nexttile
    medianRadius = median(result.Radius_CTL(:,1));
    histogram(result.Radius_CTL, [0:10:250],'FaceColor',result.color_CTL, 'FaceAlpha', result.alpha_CTL);
    line([medianRadius, medianRadius], [0, 22], 'Color', 'red', 'LineWidth', 2, 'LineStyle', '--');
    title('CTL Pit-shape cluster radius'); xlabel('cluster radius (nm)'); ylabel('counts');


    nexttile
    medianRadius = median(result.Radius_OSMO(:,1));
    histogram(result.Radius_OSMO, [0:10:250],'FaceColor',result.color_OSMO, 'FaceAlpha', result.alpha_OSMO);
    line([medianRadius, medianRadius], [0, 22], 'Color', 'red', 'LineWidth', 2, 'LineStyle', '--');
    title('OSMO Pit-shape cluster radius'); xlabel('cluster radius (nm)'); ylabel('counts');

end