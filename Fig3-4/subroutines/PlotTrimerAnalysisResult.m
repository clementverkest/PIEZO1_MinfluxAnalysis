
function [result] = PlotTrimerAnalysisResult(result, DataSource)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot analysis results
if size(result.Trimers,1) > 0
screenSize = get(0, 'ScreenSize');
figure('Position',[0 50 1000 300]);
tiledlayout(1,7);
    if DataSource == 0
        condition = 'soma';
    elseif DataSource == 1
        condition = 'neurite'
    else condition = 'cytoD'
    end   
annotation('textbox', [0.02 0.5 0.2 0.1], 'String', condition, ...
    'FitBoxToText', 'on', ...
    'EdgeColor', 'none', ...
    'FontSize', 12);
%% plot mean interblade
nexttile
    result.InterBlades = cat(2, result.InterBlades, mean(result.InterBlades(:,1:3),2));
    bar(mean(result.InterBlades(:,4)),'FaceColor',result.color);
    hold on;
    swarmchart(ones(size(result.InterBlades,1)),result.InterBlades(:,4),'k','filled');
    ylabel('interblade distance(nm)');
    title('interblade dist.');

%% plot interblade angle histogram
    result.AnglesSingleColumn = cat(1, result.Angles(:,1), result.Angles(:,2),result.Angles(:,3));
    nexttile
    histogram(result.AnglesSingleColumn, [0:10:120],'FaceColor',result.color);
    title('angle distrib.'); xlabel('interblade angle (°)'); ylabel('counts');

%% plot superparticle
[result.Superparticle] = PIEZO1Superparticle(result); % second parameter: choose 1 align or no align
nexttile(3,[1,2])  
hold on;
    scatter(result.Superparticle(:,1), result.Superparticle(:,2),200, "o", "filled", 'MarkerFaceColor', result.color, 'MarkerFaceAlpha',[0.3]);
    scatter(result.Superparticle(:,3), result.Superparticle(:,4),200, "o", "filled", 'MarkerFaceColor', result.color, 'MarkerFaceAlpha',[0.3]);
    scatter(result.Superparticle(:,5), result.Superparticle(:,6),200, "o", "filled", 'MarkerFaceColor', result.color, 'MarkerFaceAlpha',[0.3]);
    grid on; axis equal; xlim([-30 30]); ylim([-30 30]);
    xlabel('(nm)');
    ylabel('(nm)');
 title('all trimers overlay')
    %% STD 
    nexttile
 histogram(result.AllStdTrimer(:,1),'BinWidth', 0.5,'FaceColor',result.color);
       xlim([0 10])
       title("STD X");
        ax=gca;
       ax.TickDir="out";
  xlabel('\sigma X(nm)')
  ylabel('counts')
  nexttile
histogram(result.AllStdTrimer(:,2),'BinWidth', 0.5,'FaceColor',result.color);
       xlim([0 10]) ;   
       title("STD Y");
        ax=gca;
       ax.TickDir="out";
  xlabel('\sigma Y(nm)')
  nexttile
histogram(result.AllStdTrimer(:,3),'BinWidth', 0.5,'FaceColor',result.color);
       xlim([0 10]) ;
       title("STD Z");
        ax=gca;
       ax.TickDir="out";
 xlabel('\sigma Z(nm)')


end
