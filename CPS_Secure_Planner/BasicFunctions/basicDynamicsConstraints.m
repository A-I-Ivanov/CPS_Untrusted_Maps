function [cin, eqconst, DCIn, DCeq] = basicDynamicsConstraints(x)
    finDiffDelta = 1e-6;
global nx nu K

    E = eye(nx);
    numConst = 1; 
    numEqConst = nx*(K-1);
    eqconst = zeros(numEqConst,1);
    
    DCeq = zeros(numEqConst, length(x)); 
    DCIn = zeros(numConst*(K), length(x)); 
    

    
    cin = zeros(numConst*(K),1);
    deltaT = x(1);
    
    %Define state-based constraints and visibility constraints
    for i=1:K-1
        xNow = x((nx+nu)*(i-1)+2:(nx+nu)*(i-1)+1+nx);
        xNext = x((nx+nu)*(i)+2:(nx+nu)*(i)+1+nx);
        uNow = x((nx+nu)*(i-1)+2+nx:(nx+nu)*(i)+1);
        eqconst((i-1)*(nx)+1:i*nx) = xNext - diffDriveKinematics(xNow, uNow, deltaT);
        
        %%Populate sparse jacobian matrix for differential drive
        DCeq((i-1)*(nx)+1:i*nx,:)=popSparseJac(i, DCeq((i-1)*(nx)+1:i*nx,:), xNow, uNow, deltaT);
        
        
        cin(numConst*(i-1)+1:numConst*(i)) = calcDistConst(xNow);
        DCIn(numConst*(i-1)+1:numConst*(i),:) = popSparseInJac(i,DCIn(numConst*(i-1)+1:numConst*(i),:), xNow);
        
        
    end
    
   
    
    DCIn = DCIn';
    DCeq = DCeq';
    
    return
    
    
    
    %%%%%%%Begin Helper Functions%%%%%%%%%%%
    
    %This function populates the sparse inequality jacobian using finite
    %differences
    function rows = popSparseInJac(i,DinRows, xNow)
        rows = DinRows;
        for j=1:nx
            vectDelt = finDiffDelta*E(:,j);
            cDeltPlus = calcDistConst(xNow + vectDelt)';
            cDeltMinus = calcDistConst(xNow - vectDelt)';
            rows(:,1+j+(i-1)*(nx+nu))= (cDeltPlus-cDeltMinus)/(2*finDiffDelta);
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


    function distConst = calcDistConst(xNow)
        
        distConst = -obstDistance(xNow(1:2), [], xNow(4)); 
        
    end

   
end



