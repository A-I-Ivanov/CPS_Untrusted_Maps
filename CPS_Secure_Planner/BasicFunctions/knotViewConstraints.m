%%%%Written by Alexander I. Ivanov - 2017%%%%
function [cin, eqconst, DCIn, DCeq] = knotViewConstraints(x)
    finDiffDelta = 1e-6;
    
%This function provides dynamic defect constraints as well as 
%visibility constraints and distance constraints for the reactive 
%set
global nx nu K rSensor thetaSensor num2 num1

    E = eye(nx);
    numConst = 8; 
    numEqConst = nx*(K-1);
    deltaEye = 5;
    eqconst = zeros(numEqConst,1);
    
    DCeq = zeros(numEqConst, length(x)); 
    DCIn = zeros(numConst*(K-deltaEye)+deltaEye, length(x)); %+deltaEye
    

    
    cin = zeros(numConst*(K-deltaEye)+deltaEye,1);
    deltaT = x(1);
    
    %Define state-based constraints and visibility constraints
    for i=1:K-deltaEye-1
        xNow = x((nx+nu)*(i-1)+2:(nx+nu)*(i-1)+1+nx);
        xNext = x((nx+nu)*(i)+2:(nx+nu)*(i)+1+nx);
        uNow = x((nx+nu)*(i-1)+2+nx:(nx+nu)*(i)+1);
        eqconst((i-1)*(nx)+1:i*nx) = xNext - diffDriveKinematics(xNow, uNow, deltaT);
        
        %%Populate sparse jacobian matrix for differential drive
        DCeq((i-1)*(nx)+1:i*nx,:)=popSparseJac(i, DCeq((i-1)*(nx)+1:i*nx,:), xNow, uNow, deltaT);
        
        
        x2C = x((nx+nu)*(i-1+deltaEye)+2:(nx+nu)*(i-1+deltaEye)+1+nx);
        
        cin(numConst*(i-1)+1:numConst*(i)) = calcInEq(xNow, x2C);
        DCIn(numConst*(i-1)+1:numConst*(i),:) = popSparseInJac(i,DCIn(numConst*(i-1)+1:numConst*(i),:), cin(numConst*(i-1)+1:numConst*(i)), xNow, x2C);
        
        
    end
    
    %The following are just state-based constraints
    l = 0; %this makes other indexing easier
    for i=K-deltaEye:K-1
        l = l+1;
        xNow = x((nx+nu)*(i-1)+2:(nx+nu)*(i-1)+1+nx);
        xNext = x((nx+nu)*(i)+2:(nx+nu)*(i)+1+nx);
        uNow = x((nx+nu)*(i-1)+2+nx:(nx+nu)*(i)+1);
        eqconst((i-1)*(nx)+1:i*nx) = xNext - diffDriveKinematics(xNow, uNow, deltaT);
        
        %%Populate sparse jacobian matrix for differential drive
        DCeq((i-1)*(nx)+1:i*nx,:)=popSparseJac(i, DCeq((i-1)*(nx)+1:i*nx,:), xNow, uNow, deltaT);
        
        cVal =  calcDistConst(xNow);
        cin(numConst*(K-deltaEye)+l) = cVal;
        
        %Take care of terminal inequality constraints
        %These have no analytic form so we do finite differences 
        for k=1:nx
            vDelt = finDiffDelta*E(:,k);
            deltDistPlus = calcDistConst(xNow +vDelt); %Check distance to obstacles
            deltDistMinus = calcDistConst(xNow -vDelt); %Check distance to obstacles
            DCIn(numConst*(K-deltaEye)+l,(nx+nu)*(i-1)+1+k) = (deltDistPlus - deltDistMinus)/(2*finDiffDelta);
        end
        
        
    end
    
    DCIn = DCIn';
    DCeq = DCeq';
    
    return
    
    
    
    %%%%%%%Begin Helper Functions%%%%%%%%%%%
    
    %This function populates the sparse inequality jacobian using finite
    %differences
    function rows = popSparseInJac(i,DinRows, cin, xNow, x2C)
        rows = DinRows;
        for j=1:nx
            vectDelt = finDiffDelta*E(:,j);
            cDeltPlus = calcInEq(xNow + vectDelt,x2C)';
            cDeltMinus = calcInEq(xNow - vectDelt,x2C)';
            rows(:,1+j+(i-1)*(nx+nu))= (cDeltPlus-cDeltMinus)/(2*finDiffDelta);
            cDeltPlus = calcInEq(xNow,x2C  + vectDelt)';
            cDeltMinus = calcInEq(xNow, -vectDelt+x2C)';
            rows(:,1+j+(i-1+deltaEye)*(nx+nu))= (cDeltPlus-cDeltMinus)/(2*finDiffDelta);
        end

    end
    
    %This function populates rows of the sparse equality jacobian analytically
    function rows = popSparseJac(i, DeqRows, xNow, uNow, deltaT)
        C = cos(xNow(3));
        S = sin(xNow(3));
        V = xNow(4);
        rows = DeqRows;
        
        rows(1,1) = -V*C;  
        rows(1,(i-1)*(nx+nu)+2) = -1;
        rows(1,(i-1)*(nx+nu)+4) = V*S*deltaT;
        rows(1,(i-1)*(nx+nu)+5) = -C*deltaT;
        rows(1,(i-1)*(nx+nu)+9) = 1;
        
        rows(2,1) = -V*S;
        rows(2,(i-1)*(nx+nu)+3) = -1;
        rows(2,(i-1)*(nx+nu)+4) = -V*C*deltaT;
        rows(2,(i-1)*(nx+nu)+5) = -S*deltaT;
        rows(2,(i-1)*(nx+nu)+10) = 1;
        
        rows(3,1) = -xNow(5);  
        rows(3,(i-1)*(nx+nu)+4) = -1;
        rows(3,(i-1)*(nx+nu)+6) = -deltaT;
        rows(3,(i-1)*(nx+nu)+11) = 1;
        
        rows(4,1) = -uNow(1);
        rows(4,(i-1)*(nx+nu)+5) = -1;
        rows(4,(i-1)*(nx+nu)+7) = -deltaT;
        rows(4,(i-1)*(nx+nu)+12) = 1;
        
        rows(5,1) = -uNow(2);
        rows(5,(i-1)*(nx+nu)+6) = -1;
        rows(5,(i-1)*(nx+nu)+8) = -deltaT;
        rows(5,(i-1)*(nx+nu)+13) = 1;
    end


    function ineqs = calcInEq(xNow, x2C)
      
        if(abs(x2C(3))>pi)
            x2C(3) = wrapToPi(x2C(3)); %angle wrap
        end
        
        if(abs(xNow(3))>pi)
            xNow(3) = wrapToPi(xNow(3)); %angle wrap
        end
        
        ineqs(1) = calcDistConst(xNow); %Check distance of elipse to obstacles

        %The remaining eqconst are safety visual contraints
        %They utilze eliptic transformations to take the reactive 
        %set to the unit circle. This is script F in the paper
        [a,Q] = interpReactive(x2C(4));
        COS = cos(x2C(3)); SIN = sin(x2C(3));
        rotMat = [COS -SIN; 
                SIN COS];
        aRot = rotMat*a;
        xDiff = xNow(1:2)- x2C(1:2) - aRot;
        
        if(true)
        plot(xNow(1), xNow(2), 'rx');
        hold on
        [m,n] = pol2cart(x2C(3),.5);
        quiver(x2C(1)+aRot(1), x2C(2)+aRot(2),m,n);
        end
        
        QTic = rotMat*Q*rotMat';
        
        
        %%debug
        if(true)
         t = linspace(-pi,pi,100);
        circle = [cos(t); sin(t)];
        COS = cos(x2C(3)); SIN = sin(x2C(3));
        rotMat = [COS -SIN; 
                  SIN COS];
        ellipse = rotMat*chol(Q)*circle;
        ellipse(1,:) = ellipse(1,:)+x2C(1)+aRot(1);
        ellipse(2,:) = ellipse(2,:)+x2C(2)+aRot(2);
        plot(ellipse(1,: ), ellipse(2,:),'r');
        plot([xNow(1), xNow(1) + rSensor*cos(xNow(3)+thetaSensor)],[xNow(2), xNow(2) + rSensor*sin(xNow(3)+thetaSensor)],'g');
        plot([xNow(1), xNow(1) + rSensor*cos(xNow(3)-thetaSensor)],[xNow(2), xNow(2) + rSensor*sin(xNow(3)-thetaSensor)],'g');
        plot([xNow(1) + rSensor*cos(xNow(3)+thetaSensor), xNow(1) + rSensor*cos(xNow(3)-thetaSensor)],[xNow(2) + rSensor*sin(xNow(3)+thetaSensor), xNow(2) + rSensor*sin(xNow(3)-thetaSensor)],'g');
        end
        
        
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
        
        
        
        slopes = yTangent/abs(abs(xTangent) - norm(x2Tic));
        plotTri = [x2Tic, rot2xTic(1:2,1:2)*[point_up + [1; slopes;], point_down + [1; -slopes;]]];
        rescaledTri = chol(QTic)*OccludeTriangle;
        rescaledTri(1,:) = rescaledTri(1,:)+x2C(1) + aRot(1);
        rescaledTri(2,:) = rescaledTri(2,:)+x2C(2) + aRot(2);
        
        ineqs(8) = calcOcculusionConst(rescaledTri)
        
        plot([rescaledTri(1,:),rescaledTri(1,1)] , [rescaledTri(2,:),rescaledTri(2,1)], 'b');
        
        
        %%Begin testing of polyhedral notation 
        
        normals = [-sin(thetaSensor) -sin(thetaSensor)  1;
           -cos(thetaSensor) cos(thetaSensor)   0;];
        b = [0;0;cos(thetaSensor)*rSensor;];

        COS = cos(xNow(3)); SIN = sin(xNow(3));
        rot_normals = [COS -SIN; 
                       SIN COS]*normals;
                   
        %%Check ellipsoidal contianment
        if(true)
        testNorms = [rot_normals(1,:);rot_normals(2,:)];
        plot([xNow(1),testNorms(1,1)+xNow(1)],[xNow(2),testNorms(2,1)+xNow(2)],'r')
        plot([xNow(1),testNorms(1,2)+xNow(1)],[xNow(2),testNorms(2,2)+xNow(2)],'r')
        plot([xNow(1),testNorms(1,3)+xNow(1)],[xNow(2),testNorms(2,3)+xNow(2)],'r')
        
        
        axis([-2,2,-2,2])
        end
        
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
        
         hold off
        
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

     function distConst = calcOcculusionConst(triangle)
            
            distConst =  -Dist2OccTri(triangle); %Check distance of elipse to obstacles

        end

   
end



