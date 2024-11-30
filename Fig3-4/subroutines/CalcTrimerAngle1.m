%%  This function calculates the angles of a triangle and returns them as a 1 X 3 array
    % – Sep 11th 2024, Stefan G. Lechner (s.lechner@uke.de) –
    %   Input arguments: 
    %   PointA, PointB and PointC as 1 x 3 double

function [TrimerAngles] = CalcTrimerAngle1(PointA, PointB, PointC)
    AB = PointA-PointB;
    AC = PointA-PointC;
    BC = PointB-PointC;
    normAB = AB/norm(AB);
    normAC = AC/norm(AC);
    normBC = BC/norm(BC);

    dotProdA = dot(normAB, normAC);
    dotProdB = dot(-normAB, normBC);
    dotProdC = dot(-normBC, -normAC);

    angleRadA = rad2deg(acos(dotProdA));
    angleRadB = rad2deg(acos(dotProdB));
    angleRadC = rad2deg(acos(dotProdC));
    TrimerAngles = [angleRadA angleRadB angleRadC];
end
