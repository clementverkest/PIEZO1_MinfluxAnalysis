%% Align MINFLUX localisations of PIEZO trimers to a reference PIEZO trimer using the iterative closest point algorithm to generate a superparticle for presentation purposes
% – Sep 11th 2024, Stefan G. Lechner (s.lechner@uke.de) –
% Uses the 'AllTrimerCoord' array created by AFLA_analysis_finalOPT as the input 
% The script uses the pcregistericp function of Matlab, which is part of
% the computer vision tool box that must be installed.

function [Superparticle] = PIEZO1Superparticle(result);

%% create empty result array
    result1 = [];

%% choose options
    PlotSuperparticle = false;

%% read coordinates
    AllTrimerCoord = result.TrimerCoord;
    locID = [1:3*size(AllTrimerCoord,1)]';

%% select if final rotation should be applied to align the closest points to vertex A (1 = perform additional rotation alignment; 0 = no additional alignment)
    
        option = 1;
    
%% calculate coordinates of the references trimer
    P2dist=26; % set PIEZO reference trimer interblade distance; this value does not affect the results
    hP2trimer=(sqrt(3)/2)*P2dist; %calcluate heigth of triangle
    rP2trimer=(sqrt(3)/3)*P2dist; %calculate radius of triangle (distance from centroid to vertex
    refA = [(-P2dist/2) (-(hP2trimer-rP2trimer)) 0];
    refB = [(P2dist/2) (-(hP2trimer-rP2trimer)) 0];
    refC = [0 rP2trimer 0];
    refTrimer=[refA; refB; refC];
    ptCloudFix = pointCloud(refTrimer);
    plotX=[refTrimer(:,1);refTrimer(1,1)];
    plotY=[refTrimer(:,2);refTrimer(1,2)];

%% Create rotation matrices for optional rotation alignment
    R120 = [cosd(120) -sind(120); sind(120) cosd(120)];
    R240 = [cosd(240) -sind(240); sind(240) cosd(240)];

%% loop through data  
for i=1:size(AllTrimerCoord,1)
    A = AllTrimerCoord(i,1:3);
    B = AllTrimerCoord(i,4:6);
    C = AllTrimerCoord(i,7:9);
    
    % transform to 2D – i.e. rotate trimer coordinates into X-Y plane
    lengthA=norm(B-A);
    lengthB=norm(C-A);
    lengthC=norm(B-C);
    Cx=(lengthB^2-lengthC^2+lengthA^2)/(2*lengthA);
    Cy=sqrt((lengthB^2-Cx^2));
    Trimer =[0 0 0; lengthA 0 0; Cx Cy 0];
    
    % run iterative closest point algorithm   
    ptCloudmove = pointCloud(Trimer);
    tform = pcregistericp(ptCloudmove,ptCloudFix);   
    ptCloudOut = pctransform(ptCloudmove,tform);
    tformTrimer = transpose(ptCloudOut.Location.');    
    tformA=tformTrimer(1,1:2);
    tformB=tformTrimer(2,1:2);
    tformC=tformTrimer(3,1:2);
    
    % rotate such that closest of the three points is A

    if option==1    
    distA=norm(tformTrimer(1,:)-refA);
    distB=norm(tformTrimer(2,:)-refB);
    distC=norm(tformTrimer(3,:)-refC);
        if distA<distB && distA<distC % do nothing
            result1=[result1; tformA; tformB; tformC];
        elseif distB<distA && distB<distC % rotate by 120           
            tformA=transpose(R120*tformA');
            tformB=transpose(R120*tformB');
            tformC=transpose(R120*tformC');
            result1=[result1; tformC; tformA; tformB];
        elseif distC<distB && distC<distA %rotate by 240°           
            tformA=transpose(R240*tformA');
            tformB=transpose(R240*tformB');
            tformC=transpose(R240*tformC');  
            result1=[result1; tformB; tformC; tformA];
        end
    else

    result1=[result1; tformA; tformB; tformC];
    end
   
end

% split result
    resultA=result1(mod(locID+2,3)==0,:);
    resultB=result1(mod(locID+1,3)==0,:);
    resultC=result1(mod(locID,3)==0,:);

%% plot results

% create colorscale for distance from centroid
if PlotSuperparticle
    Avg_A=mean(resultA);
    Avg_B=mean(resultB);
    Avg_C=mean(resultC);
    ColorA=resultA-Avg_A;
    ColorB=resultB-Avg_B;
    ColorC=resultC-Avg_C;

    ColorScaleA=vecnorm(ColorA,2,2);
    ColorScaleB=vecnorm(ColorB,2,2);
    ColorScaleC=vecnorm(ColorC,2,2);
    ColorScale=[ColorScaleA;ColorScaleB;ColorScaleC];

  
    [ColorScaleA,sortIdx] = sort(ColorScaleA, 'descend');
    resultA = resultA(sortIdx,:);
 
    [ColorScaleB,sortIdx] = sort(ColorScaleB, 'descend');
    resultB = resultB(sortIdx,:);

    [ColorScaleC,sortIdx] = sort(ColorScaleC, 'descend');
    resultC = resultC(sortIdx,:);

    f2=figure;
    hold on;
    scatter(resultA(:,1), resultA(:,2),400, ColorScaleA, "o", "filled");%, 'MarkerFaceColor', [1 0 0], 'MarkerFaceAlpha',[0.5]);
    scatter(resultB(:,1), resultB(:,2),400, ColorScaleB, "o", "filled");%, 'MarkerFaceColor', [0 1 0], 'MarkerFaceAlpha',[0.5]);
    scatter(resultC(:,1), resultC(:,2),400, ColorScaleC, "o", "filled");%, 'MarkerFaceColor', [0 0 1], 'MarkerFaceAlpha',[0.5]);

    plot(plotX, plotY, '--o','color','red','LineWidth',2);
    grid on; axis equal; xlim([-30 30]); ylim([-30 30]);
    colormap('parula'); colorbar;
    c = colorbar; c.Label.String = 'Distance from centroid (nm)';
    c.Label.FontSize = 14;
    title('PIEZO1 superparticle (overlay of all identified trimers'); xlabel('x (nm)'); ylabel('y (nm)');
end

    % return superparticle coordinates
    Superparticle = [resultA resultB resultC];

end
