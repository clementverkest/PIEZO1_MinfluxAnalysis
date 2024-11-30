function PlotRawData(traces, myfile, status)

% traces = traces(traces(:,8)<3600,:);

fig1 = figure('Position', [0 300 1200 400]);
tiledlayout(2,5);
% set(gcf,'renderer','Painters');
nexttile(1,[2 2]);
        scatter3(traces(:,1), traces(:,2), traces(:,3),5,traces(:,3), 'filled');  % color by z
        colormap parula
        axis equal
        title(strcat(myfile, status, 'data'));
        ax=gca;
       ax.TickDir="out";
nexttile
       histogram(traces(:,6)); %CFR
       title(strcat('CFR histogram'));
        ax=gca;
       ax.TickDir="out";
       ylabel('count')
nexttile
      histogram(traces(:,5),200); %EFO
       title(strcat('EFO histogram'));
        ax=gca;
       ax.TickDir="out";
xlabel('Freq (Hz)')
ylabel('count')
nexttile
    [A1, ~, A2] = unique(traces(:,4)); 
    [A3,~] = histcounts(A2,size(A1,1)); 
    LocDistribution = A3';
    histogram(LocDistribution, 'BinWidth', 1); %number of localisations per trace
    title(strcat('loc/trace histogram'));
    [uv_tid, ~, id_tid] = unique(traces(:,4));
    StDevALL=[accumarray(id_tid,traces(:,1),[],@std) accumarray(id_tid,traces(:,2),[],@std) accumarray(id_tid,traces(:,3),[],@std)];
    StDevHist = StDevALL(StDevALL(:,1)>0,:);  
        ax=gca;
       ax.TickDir="out";
    
nexttile
       histogram(StDevHist(:,1),'BinWidth', 0.5);
       xlim([0 10])
       title("STD-X");
        ax=gca;
       ax.TickDir="out";
       xlabel('std (nm)')
nexttile
        histogram(StDevHist(:,2),'BinWidth', 0.5);
        xlim([0 10])
        title("STD-Y");
         ax=gca;
       ax.TickDir="out";
       xlabel('std (nm)')
nexttile
        histogram(StDevHist(:,3),'BinWidth', 0.5);
        xlim([0 10])
        title("STD-Z");
        ax=gca;
       ax.TickDir="out";
       xlabel('std (nm)')
%% optional code for plotting raw data in 2D as fig 10 for measring neurite diameters with getcorsor.m script
    %     fig1=10;
    % figure(fig1);
    %     scatter(traces(:,1), traces(:,2),10,'red', 'filled');  % color options: ClusterIDs, myXYZ(:,5)
    %     axis equal




% Dateiname = replace(myfile,".mat","_SOMA_RAW.fig");
% savefig(Dateiname)
end