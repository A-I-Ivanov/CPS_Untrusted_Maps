function [ R ] = rMatrix(thetaX, thetaY, thetaZ )
%RMATRIX Summary of this function goes here
%   Detailed explanation goes here
   R = (rotz(thetaZ))*(roty(thetaY))*(rotx(thetaX));
   
    function Rz = rotz(theta)
        Cz = cos(theta); Sz = sin(theta);
        Rz = [ Cz, -Sz, 0;
               Sz,  Cz, 0;
                0,   0, 1;];
        
    end

    function Rx = rotx(theta)
        Cx = cos(theta); Sx = sin(theta);
        Rx = [ 1,  0,   0;
              0, Cx, -Sx;
              0, Sx,  Cx];
    end
    function Ry = roty(theta)
        Cy = cos(theta); Sy = sin(theta);
        Ry = [  Cy, 0, Sy;
                 0, 1,  0;
               -Sy, 0, Cy;];
        
    end
end

