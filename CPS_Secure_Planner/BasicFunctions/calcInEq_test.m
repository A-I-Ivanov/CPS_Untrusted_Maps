function ineqs = calcInEq_test(xNow, x2C)
      
DeltaI = 5;
deltaT = .1;
thetaSensor = .7854;
rSensor = 2;
        if(abs(x2C(3))>pi)
            x2C(3) = wrapToPi(x2C(3)); %angle wrap
        end
        
        if(abs(xNow(3))>pi)
            xNow(3) = wrapToPi(xNow(3)); %angle wrap
        end
        
        

        %The remaining eqconst are safety visual contraints
        %They utilze eliptic transformations to take the reactive 
        %set to the unit circle. This is script F in the paper
        [a,Q] = interpReactive(x2C(4));
        COS = cos(x2C(3)); SIN = sin(x2C(3));
        rotMat = [COS -SIN; 
                SIN COS];
        aRot = rotMat*a;
        xDiff = xNow(1:2)- x2C(1:2) - aRot;
        ineqs(1) = calcDistConst(x2C); %Check distance of elipse to obstacles

        QTic = rotMat*Q*rotMat';
        
        
        
        sqrtQTic = inv(chol(QTic));
        x2Tic = sqrtQTic*xDiff;
        
        
        %D = x cross normal vector
        
        len = sqrt(abs(1-norm(x2Tic)^2));
        sinThe = 1/norm(x2Tic);
        cosThe = sqrt(abs(1-sinThe^2));
        
        xTangent = len*cosThe - norm(x2Tic);
        yTangent = len*sinThe;
        point_up   = [xTangent; yTangent];
        point_down = [xTangent; -yTangent];
        
        rot2xTic = rotz(360*cart2pol(xDiff(1), xDiff(2))/2/pi+180);
        OccludeTriangle = [x2Tic, rot2xTic(1:2,1:2)*[point_up, point_down], x2Tic];
        
        
        
% %         slopes = yTangent/abs(abs(xTangent) - norm(x2Tic));
% %         plotTri = [x2Tic, rot2xTic(1:2,1:2)*[point_up + [1; slopes;], point_down + [1; -slopes;]]];
        rescaledTri = chol(QTic)*OccludeTriangle;
        rescaledTri(1,:) = rescaledTri(1,:)+x2C(1) + aRot(1);
        rescaledTri(2,:) = rescaledTri(2,:)+x2C(2) + aRot(2);
        
%        ineqs(8) = calcOcculusionConst(rescaledTri);
        
% %         plot([rescaledTri(1,:),rescaledTri(1,1)] , [rescaledTri(2,:),rescaledTri(2,1)], 'b');
        
        
        %%Begin testing of polyhedral notation 
        
        normals = [-sin(thetaSensor) -sin(thetaSensor)  1;
           -cos(thetaSensor) cos(thetaSensor)   0;];
        b = [0;0;cos(thetaSensor)*rSensor;];

        COS = cos(xNow(3)); SIN = sin(xNow(3));
        rot_normals = [COS -SIN; 
                       SIN COS]*normals;
                   
        %%Check ellipsoidal contianment
% %         if(true)
% %         testNorms = [rot_normals(1,:);rot_normals(2,:)];
% %         plot([xNow(1),testNorms(1,1)+xNow(1)],[xNow(2),testNorms(2,1)+xNow(2)],'r')
% %         plot([xNow(1),testNorms(1,2)+xNow(1)],[xNow(2),testNorms(2,2)+xNow(2)],'r')
% %         plot([xNow(1),testNorms(1,3)+xNow(1)],[xNow(2),testNorms(2,3)+xNow(2)],'r')
% %         
% %         
% %         axis([-2,2,-2,2])
% %         end
        
        test_eqs(1:3) = rot_normals'*(x2C(1:2) + aRot - xNow(1:2)) - b;
        b_translate = rot_normals'*(xNow(1:2)-(x2C(1:2)+aRot)); 
        b_new  = b + b_translate;
        
        %%ToDo check this seems wrong
        %COS = cos(xNow(3)); SIN = sin(xNow(3));
        %rot_normals = [COS -SIN; 
        %               SIN COS]*rot_normals;
        
        scew_norm = inv(chol(inv(Q)))'*rot_normals;
        
        for p =1:length(b)
            b_new(p) = b_new(p)/(norm(scew_norm(:,p)));
        end
        
        test_eqs(4:6) = 1- b_new.^2;
        
        
        
        ineqs(2:7) = test_eqs;
        
% %          hold off
        
end
    
function distConst = calcDistConst(xNow)
        if(abs(xNow(3))>pi)
            xNow(3) = wrapToPi(xNow(3)); %angle wrap
        end
        
        %The remaining eqconst are safety distance contraints
        %They utilze eliptic transformations 
        [a,Q] = interpReactive(xNow(4));
        COS = cos(xNow(3)); SIN = sin(xNow(3));
        
        rotMat = [COS -SIN; 
                  SIN COS];
        %%Has to be in sensor range 
        arot = rotMat*a;
        QTic = rotMat*Q*rotMat';
        
        sqrtQTic = chol(inv(QTic));
        
        distConst = -obstDistance(xNow(1:2)+arot, sqrtQTic, 0); %Check distance of elipse to obstacles
        
    end