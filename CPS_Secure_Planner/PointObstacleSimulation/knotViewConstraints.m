function [cin, eqconst, DCIn, DCeq] = knotViewConstraints(x)
    finDiffDelta = 1e-6;
global nx nu K rSensor thetaSensor num2 num1

    E = eye(nx);
    numConst = 6; 
    numEqConst = nx*(K-1);
    deltaEye = 4;
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
        
        
        %The remaining eqconst are safety visual contraints
        %They utilze eliptic transformations 
        [a,Q] = interpReactive(x2C(4));
        COS = cos(x2C(3)); SIN = sin(x2C(3));
        aRot = [COS -SIN; 
                SIN COS]*a;
        xDiff = xNow(1:2)- x2C(1:2) - aRot;
        
        %%Has to be in sensor range
        ineqs(2) = (norm(xDiff - aRot) - rSensor)/rSensor; 
        theta12 = -atan2(xDiff(2), xDiff(1)) + pi;
        
        COS = cos(theta12); SIN = sin(theta12);
        rotMat = [COS -SIN; 
                  SIN COS];
        xTic = rotMat*xDiff;
        
        
        %Make sure we actually have intercepts
        angUp = thetaSensor + xNow(3) + theta12;
        angDown = -thetaSensor + xNow(3) + theta12;
        if(angUp/pi-floor(angUp/pi)>(1/2)) %Fast angle wrap
            m1 = 100;
            num1 = num1+1;
        else
            m1 = tan(angUp);
        end
        
        if(angDown/pi-ceil(angDown/pi)<-(1/2)) %fast angle wrap
            m2 = -100;
            num2 = num2+1;
        else
            m2 = tan(angDown);
        end
        

        yTic1 = -[0 m1*xTic(1)]';
        yTic2 = -[0 m2*xTic(1)]';
        
        
        QTic = rotMat*Q*rotMat';
        
        sqrtQTic = chol(inv(QTic));
        
        %This is technically wrong. Should be using a different Q and a
        ineqs(1) = calcDistConst(xNow); %1-obstDistance(xNow(1:2)+aRot, sqrtQTic); %Check distance of elipse to obstacles

        
        x2Tic = sqrtQTic*xTic;
        y2Tic1 = sqrtQTic*yTic1;
        y2Tic2 = sqrtQTic*yTic2;
        
        A1 = (x2Tic-y2Tic1);
        A2 = (x2Tic - y2Tic2);
        scale = norm(x2Tic);
        
        %D = x cross normal vector
        
        %Also these should be in linear inequality constriants
        %We normalize to one to get good constraint conditioning
        ineqs(3) = (-abs(prod(A1)/norm(A1))+ 1)/scale;
        ineqs(4) = (-abs(prod(A2)/norm(A2))+ 1)/scale;
        ineqs(5) = -y2Tic1(2)/scale;
        ineqs(6) = y2Tic2(2)/scale;
    end

    function distConst = calcDistConst(xNow)
        if(abs(xNow(3))>pi)
            xNow(3) = wrapToPi(xNow(3)); %angle wrap
        end
        
        %The remaining eqconst are safety visual contraints
        %They utilze eliptic transformations 
        [a,Q] = interpReactive(xNow(4));
        COS = cos(xNow(3)); SIN = sin(xNow(3));
        
        rotMat = [COS -SIN; 
                  SIN COS];
        %%Has to be in sensor range 
        arot = rotMat*a;
        QTic = rotMat*Q*rotMat';
        
        sqrtQTic = chol(inv(QTic));
        
        distConst = -obstDistance(xNow(1:2)+arot, sqrtQTic, xNow(4)); %Check distance of elipse to obstacles
        
    end

   
end



