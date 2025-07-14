%% plot trimer raw projection
StdTrimALL = [];
fns = fieldnames(IndivTrimersRAW);

%% options
colorscheme = 3; %1 = green, 2 = purple, 3 = grey
PlotTrimer3D = false;
PlotContours = false;

for k=1:size(fns,1) 

    % extract protamer coordinates    
    TrimA = IndivTrimersRAW.(fns{k}).A;
    TrimB = IndivTrimersRAW.(fns{k}).B;
    TrimC = IndivTrimersRAW.(fns{k}).C;
    
    % calculate standard deviation of each protamer signal
    StdA = [std(TrimA(:,1)) std(TrimA(:,2)) std(TrimA(:,3))];
    StdB = [std(TrimB(:,1)) std(TrimB(:,2)) std(TrimB(:,3))];
    StdC = [std(TrimC(:,1)) std(TrimC(:,2)) std(TrimC(:,3))];
    StdTrim = [StdA; StdB; StdC];
    StdTrimALL = cat(1,StdTrimALL, StdTrim);
    
    
    % shift coordinates to center of first protamer
    CenterA = [mean(TrimA(:,1)) mean(TrimA(:,2)) mean(TrimA(:,3))];
    TrimAtrans = TrimA-CenterA;
    TrimBtrans = TrimB-CenterA;
    TrimCtrans = TrimC-CenterA;
    MyTrimer = [TrimAtrans; TrimBtrans;TrimCtrans];
    
    % create copy of trimer and shit center of trimer to origin for plotting
    CenterTrim = [mean(MyTrimer(:,1),1) mean(MyTrimer(:,2),1) mean(MyTrimer(:,3),1)];
    TrimAtrans2 = TrimAtrans-CenterTrim;
    TrimBtrans2 = TrimBtrans-CenterTrim;
    TrimCtrans2 = TrimCtrans-CenterTrim;
    
    % calculate center of each protamer
    ProtA = [mean(TrimAtrans(:,1)) mean(TrimAtrans(:,2)) mean(TrimAtrans(:,3))];
    ProtB = [mean(TrimBtrans(:,1)) mean(TrimBtrans(:,2)) mean(TrimBtrans(:,3))];
    ProtC = [mean(TrimCtrans(:,1)) mean(TrimCtrans(:,2)) mean(TrimCtrans(:,3))];
    TrimerOrig = [ProtA; ProtB; ProtC];
    
    %% rotate data to create projection to plane formed by ABC
    % rotate around Z
    a = ProtB(:,1);
    b = ProtB(:,2);
    c = sqrt(a^2+b^2);
    angle = asind(b/c);
    if a>0 && b>0
        angle = angle;
    elseif a>0 && b<0
        angle = 360+angle;
    elseif a<0 && b>0
        angle = 180-abs(angle);
    elseif a<0 && b<0
        angle = 180 + abs(angle);
    end
    
    RotMatrixZ = rotz(angle);
    ProtA = ProtA*RotMatrixZ;
    ProtB = ProtB*RotMatrixZ;
    ProtC = ProtC*RotMatrixZ;
    
    % % rotate around y
    a = ProtB(:,1);
    b= ProtB(:,3);
    c = sqrt(a^2+b^2);
    angle = -asind(b/c);
    if a>0 && b>0
        angle = angle;
    elseif a>0 && b<0
        angle = 360+angle;
    elseif a<0 && b>0
        angle = 180-abs(angle);
    elseif a<0 && b<0
        angle = 180 + abs(angle);
    end
    RotMatrixY = roty(angle);
    ProtA = ProtA*RotMatrixY;
    ProtB = ProtB*RotMatrixY;
    ProtC = ProtC*RotMatrixY;
    % 
    % rotate around X
    a = ProtC(:,2);
    b = ProtC(:,3);
    c = sqrt(a^2+b^2);
    angle = asind(b/c);
    if a>0 && b>0
        angle = angle;
    elseif a>0 && b<0
        angle = 360-abs(angle);
    elseif a<0 && b>0
        angle = 180 - abs(angle);
    elseif a<0 && b<0
        angle = 180 + abs(angle);
    end
    RotMatrixX = rotx(angle);
    ProtA = ProtA*RotMatrixX;
    ProtB = ProtB*RotMatrixX;
    ProtC = ProtC*RotMatrixX;
    
    % create array with rotated protamer means
    TrimerRot = [ProtA; ProtB; ProtC];
    
    % apply rotation matrices to indicidual localisaitons
    TrimAtrans = TrimAtrans*RotMatrixZ;
    TrimAtrans = TrimAtrans*RotMatrixY;
    TrimAtrans = TrimAtrans*RotMatrixX;
    
    TrimBtrans = TrimBtrans*RotMatrixZ;
    TrimBtrans = TrimBtrans*RotMatrixY;
    TrimBtrans = TrimBtrans*RotMatrixX;
    
    TrimCtrans = TrimCtrans*RotMatrixZ;
    TrimCtrans = TrimCtrans*RotMatrixY;
    TrimCtrans = TrimCtrans*RotMatrixX;
    
    % calculate center positions of final rotated trimer (projection)
    CenterTrimTransA = [mean(TrimAtrans(:,1)) mean(TrimAtrans(:,2)) mean(TrimAtrans(:,3))];
    CenterTrimTransB = [mean(TrimBtrans(:,1)) mean(TrimBtrans(:,2)) mean(TrimBtrans(:,3))];
    CenterTrimTransC = [mean(TrimCtrans(:,1)) mean(TrimCtrans(:,2)) mean(TrimCtrans(:,3))];
    AllCenters =[CenterTrimTransA; CenterTrimTransB; CenterTrimTransC];
    
    % additional rotation for figure
    if k==1
    RotMatrixZ = rotz(135);
    TrimAtrans = TrimAtrans*RotMatrixZ;
    TrimBtrans = TrimBtrans*RotMatrixZ;
    TrimCtrans = TrimCtrans*RotMatrixZ;
    AllCenters = AllCenters*RotMatrixZ;
    end
    
    % shift final coordinates such that the trimer center is in the origin
    FinalCenter = mean(AllCenters,1);
    TrimAtrans = TrimAtrans-FinalCenter;
    TrimBtrans = TrimBtrans-FinalCenter;
    TrimCtrans = TrimCtrans-FinalCenter;
    AllCenters = AllCenters-FinalCenter;
    
    %% fit gaussian
    % Define the range for the plot
    MyGridSpacing = 61;
    x = linspace(-30, +30, MyGridSpacing);
    y = linspace(-30, +30, MyGridSpacing);

    x_a = linspace(min(TrimAtrans(:,1)), max(TrimAtrans(:,1)), MyGridSpacing);
    y_a = linspace(min(TrimAtrans(:,2)), max(TrimAtrans(:,2)), MyGridSpacing);
    x_b = linspace(min(TrimBtrans(:,1)), max(TrimBtrans(:,1)), MyGridSpacing);
    y_b = linspace(min(TrimBtrans(:,2)), max(TrimBtrans(:,2)), MyGridSpacing);
    x_c = linspace(min(TrimCtrans(:,1)), max(TrimCtrans(:,1)), MyGridSpacing);
    y_c = linspace(min(TrimCtrans(:,2)), max(TrimCtrans(:,2)), MyGridSpacing);

    [X, Y] = meshgrid(x, y);
    TrimAtrans2D = TrimAtrans(:,1:2);
    if size(TrimAtrans2D,1)>2
        GausA = fitgmdist(TrimAtrans2D,1);
        % Create a grid
        [Xa, Ya] = meshgrid(x, y);
        % Combine the grid points for evaluation
        XY = [Xa(:) Ya(:)];
        % Evaluate the GMM PDF at the grid points
        Za = reshape(pdf(GausA, XY), size(Xa));
    end
    
    TrimBtrans2D = TrimBtrans(:,1:2);
    if size(TrimBtrans2D,1)>2
        GausB = fitgmdist(TrimBtrans2D,1);
        % Create a grid
        [Xb, Yb] = meshgrid(x, y);
        % Combine the grid points for evaluation
        XY = [Xb(:) Yb(:)];
        % Evaluate the GMM PDF at the grid points
        Zb = reshape(pdf(GausB, XY), size(Xb));
    end
    
    TrimCtrans2D = TrimCtrans(:,1:2);
    if size(TrimCtrans2D,1)>2
        GausC = fitgmdist(TrimCtrans2D,1);
        % Create a grid
        [Xc, Yc] = meshgrid(x, y);
        % Combine the grid points for evaluation
        XY = [Xc(:) Yc(:)];
        % Evaluate the GMM PDF at the grid points
        Zc = reshape(pdf(GausC, XY), size(Xc));
    end
    
    % combine Z-values fo gaussians for plotting in single graph
    AllZ = Za+Zb+Zc;
    % AllZ = zeros(30);
    % AllZ = zeros+Za+Zb+Zc;
    
    % define sigma levels for optional contour plot
    sigma_levels = [0.607, 0.135, 0.011];
    
    % colorschemes
    if colorscheme == 1
        color1 =[0 0.3 0.1];
        color2 = [0 0.6 0];
        color3 = [0 1 0];
    elseif colorscheme == 2
        color1 =[1 0 1];
        color2 = [0.6 0 1];
        color3 = [1 0.5 1];
    elseif colorscheme == 3
        color2 = 'k';
        color3 = [0.7 0.7 0.7];
        color1 = [0.4 0.4 0.4];
    elseif colorscheme == 4
        color1 = [0 0.5 1];
        color2 = [0.5 0.7 0.9];
        color3 = [0.25 0.4 0.55];
    end 
    
    
    % plot 3D trimer - optional
    if PlotTrimer3D
        figure('Position',[100 400 220 180]);
         set(gcf,'renderer','Painters');
        scatter3(TrimAtrans2(:,1),TrimAtrans2(:,2),TrimAtrans2(:,3),120,color1,'filled','MarkerFaceAlpha',0.7);
        hold on;
        scatter3(TrimBtrans2(:,1),TrimBtrans2(:,2),TrimBtrans2(:,3),120,color2,'filled','MarkerFaceAlpha',0.7);
        scatter3(TrimCtrans2(:,1),TrimCtrans2(:,2),TrimCtrans2(:,3),120,color3,'filled','MarkerFaceAlpha',0.7);
        axis tight; axis equal;
        xlim([-30 30]); ylim([-30 30]); zlim ([-30 30]);
        xlabel('x (nm)'); ylabel('y (nm)'); zlabel('z (nm)');
                ax = gca;
                ax.XAxis.FontSize = 14; ax.YAxis.FontSize = 14; ax.ZAxis.FontSize = 14;
                ax.GridColor = [0.8 0.8 0.8];
                ax.GridAlpha = 1;
                ax.GridLineWidth = 0.25;
    end
    
    %% plot results
    figure('Position',[100 100 650 200]);
    tiledlayout(1,3);
    nexttile
    scatter3(TrimAtrans2(:,1),TrimAtrans2(:,2),TrimAtrans2(:,3),120,color1,'filled','MarkerFaceAlpha',0.7);
        hold on;
        scatter3(TrimBtrans2(:,1),TrimBtrans2(:,2),TrimBtrans2(:,3),120,color2,'filled','MarkerFaceAlpha',0.7);
        scatter3(TrimCtrans2(:,1),TrimCtrans2(:,2),TrimCtrans2(:,3),120,color3,'filled','MarkerFaceAlpha',0.7);
        axis tight; axis equal;
        xlim([-30 30]); ylim([-30 30]); zlim ([-30 30]);
        xlabel('x (nm)'); ylabel('y (nm)'); zlabel('z (nm)');
                ax = gca;
                ax.XAxis.FontSize = 14; ax.YAxis.FontSize = 14; ax.ZAxis.FontSize = 14;
                ax.GridColor = [0.8 0.8 0.8];
                ax.GridAlpha = 1;
                ax.GridLineWidth = 0.25;
    nexttile
    hold on;
    scatter(TrimAtrans(:,1), TrimAtrans(:,2), 50, 'MarkerFaceColor',color1,'MarkerEdgeColor','none','MarkerFaceAlpha',0.8);
    scatter(TrimBtrans(:,1), TrimBtrans(:,2),50,'MarkerFaceColor',color2,'MarkerEdgeColor','none','MarkerFaceAlpha',0.8);
    scatter(TrimCtrans(:,1), TrimCtrans(:,2),50,'MarkerFaceColor',color3,'MarkerEdgeColor','none','MarkerFaceAlpha',0.8);
    
    plot(AllCenters(:,1),AllCenters(:,2),'+','MarkerEdgeColor','red', 'MarkerSize',15,'LineWidth',1);
    
    if PlotContours
        contour(Xa, Ya, Za, [sigma_levels(1) * max(Za(:)) sigma_levels(1) * max(Za(:))], 'LineWidth', 0.5, 'LineColor', 'red');
        contour(Xa, Ya, Za, [sigma_levels(2) * max(Za(:)) sigma_levels(2) * max(Za(:))], 'LineWidth', 0.5, 'LineColor', 'red');
        % % contour(Xa, Ya, Za, [sigma_levels(3) * max(Za(:)) sigma_levels(3) * max(Za(:))], 'LineWidth', 2, 'LineColor', 'red');
        % 
        contour(Xb, Yb, Zb, [sigma_levels(1) * max(Zb(:)) sigma_levels(1) * max(Zb(:))], 'LineWidth', 0.5, 'LineColor', 'red');
        contour(Xb, Yb, Zb, [sigma_levels(2) * max(Zb(:)) sigma_levels(2) * max(Zb(:))], 'LineWidth', 0.5, 'LineColor', 'red');
        % % contour(Xb, Yb, Zb, [sigma_levels(3) * max(Zb(:)) sigma_levels(3) * max(Zb(:))], 'LineWidth', 2, 'LineColor', 'red');
        % 
        contour(Xc, Yc, Zc, [sigma_levels(1) * max(Zc(:)) sigma_levels(1) * max(Zc(:))], 'LineWidth', 0.5, 'LineColor', 'red');
        contour(Xc, Yc, Zc, [sigma_levels(2) * max(Zc(:)) sigma_levels(2) * max(Zc(:))], 'LineWidth', 0.5, 'LineColor', 'red');
        % % contour(Xc, Yc, Zc, [sigma_levels(3) * max(Zc(:)) sigma_levels(3) * max(Zc(:))], 'LineWidth', 2, 'LineColor', 'red');
    end
    
    axis tight; axis equal;
    xlim([-30 30]); ylim([-30 30]);
    xlim([-31 31]); ylim([-31 31]);
    xlabel('(nm)'); ylabel('(nm)');
    ax = gca; ax.XAxis.FontSize = 14; ax.YAxis.FontSize = 14;
    ax.Box = 'on';
    ax.BoxStyle = 'full';
    ax.TickLength = [0.03,0.03];

    nexttile
    hold on;
    s=pcolor(Xa,Ya,AllZ);%,'EdgeColor','none')%,shading interp;
    s.EdgeColor = 'none';%[1 0.7 0.3];
    colormap hot; axis tight; axis equal;
    colorbar;
    xlabel('(nm)'); ylabel('(nm)');
    ax = gca; ax.XAxis.FontSize = 14; ax.YAxis.FontSize = 14;
    c = colorbar; c.Label.String = 'probability'; c.Label.FontSize = 14; c.FontSize = 14;
        ax.Box = 'on';
    ax.BoxStyle = 'full';
end
