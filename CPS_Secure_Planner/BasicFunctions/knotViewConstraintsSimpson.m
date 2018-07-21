%%%%Written by Alexander I. Ivanov - 2017%%%%
function [cin, eqconst, DCIn, DCeq] = knotViewConstraintsSimpson(x)
    finDiffDelta = 1e-6;
    
%This function provides dynamic defect constraints as well as 
%visibility constraints and distance constraints for the reactive 
%set
global nx nu K rSensor thetaSensor 

    blkSz = 2*nu+nx;
    E = eye(nx);
    numConst = 8; 
    deltaEye = 5;
    numVar = length(x);
 
    DCIn = zeros(numConst*(K-deltaEye)+deltaEye, length(x)); %+deltaEye
   
    cin = zeros(numConst*(K-deltaEye)+deltaEye,1);
    deltaT = x(1);
    
    [eqconst,DCeq] = basicDynamicsSimpsonEqonly(x);
    
    ineqOffset =1;
    colOffset =2;
    %Define state-based constraints and visibility constraints
    for i =1:deltaEye-1
        xNow = x(colOffset:colOffset+nx-1);

        cin(ineqOffset) = -obstDistance(xNow, [],0);
        DCIn(ineqOffset,:) = popDistJac(colOffset, xNow);
        ineqOffset = ineqOffset+1;
        colOffset = colOffset+blkSz;
    end
    
    colOffset =2;
    
    for i=1:K-deltaEye
        xNow = x(colOffset:colOffset+nx-1);
        x2C = x(colOffset+deltaEye*blkSz:colOffset+deltaEye*blkSz+nx-1);
        
        cin(ineqOffset:ineqOffset+numConst-1) = calcInEq(xNow, x2C);
        DCIn(ineqOffset:ineqOffset+numConst-1,:) = popSparseInJac(colOffset, xNow, x2C);
        
        colOffset = colOffset + blkSz;
        ineqOffset = ineqOffset+numConst;
    end
    

    DCIn = DCIn';
   
    
    return
    
    
    
    %%%%%%%Begin Helper Functions%%%%%%%%%%%
    
    %This function populates the sparse inequality jacobian using finite
    %differences
    function rows = popSparseInJac(colOffset, xNow, x2C)
        rows = zeros(numConst,numVar);
        for j=1:nx
            vectDelt = finDiffDelta*E(:,j);
            cDeltPlus = calcInEq(xNow + vectDelt,x2C)';
            cDeltMinus = calcInEq(xNow - vectDelt,x2C)';
            rows(:,colOffset+j-1)= (cDeltPlus-cDeltMinus)/(2*finDiffDelta);
            cDeltPlus = calcInEq(xNow,x2C  + vectDelt)';
            cDeltMinus = calcInEq(xNow, -vectDelt+x2C)';
            rows(:,colOffset+blkSz*deltaEye+j-1)= (cDeltPlus-cDeltMinus)/(2*finDiffDelta);
        end

    end
    %This function populates the sparse inequality for the first few steps using finite
    %differences
    function rows = popDistJac(colOffset, xNow)
        rows = zeros(1, numVar);
        for j=1:nx
            vectDelt = finDiffDelta*E(:,j);
            cDeltPlus = -obstDistance(xNow + vectDelt,[],0);
            cDeltMinus = -obstDistance(xNow - vectDelt,[],0);
            rows(:,colOffset+j-1)= (cDeltPlus-cDeltMinus)/(2*finDiffDelta);
        end

    end
    
    %This function populates rows of the sparse equality jacobian analytically
   


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
        
        
        

        rescaledTri = chol(QTic)*OccludeTriangle;
        rescaledTri(1,:) = rescaledTri(1,:)+x2C(1) + aRot(1);
        rescaledTri(2,:) = rescaledTri(2,:)+x2C(2) + aRot(2);
        
        ineqs(8) = calcOcculusionConst(rescaledTri);
        
        
        normals = [-sin(thetaSensor) -sin(thetaSensor)  1;
           -cos(thetaSensor) cos(thetaSensor)   0;];
        b = [0;0;cos(thetaSensor)*rSensor;];

        COS = cos(xNow(3)); SIN = sin(xNow(3));
        rot_normals = [COS -SIN; 
                       SIN COS]*normals;
                   
        %%Check ellipsoidal contianment
        
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



