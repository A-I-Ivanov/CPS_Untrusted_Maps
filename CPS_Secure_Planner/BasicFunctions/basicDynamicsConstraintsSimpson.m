%%%%Written by Alexander I. Ivanov - 2017%%%%
function [cin, eqconst, DCIn, DCeq] = basicDynamicsConstraintsSimpson(x, safetyDistance)
    finDiffDelta = 1e-6;
    
%These constraints use a Euler approximation to the defect constraints 
%of a differential drive robot. More accurate defects are required 
%for robots with fast or complex dynamics.
global K

persistent A B blkSz nx nu

    if(isempty(A))
        nx =5;
        nu =2;
       blkSz = nx+2*nu;
       %%Following bettes 142
       A = zeros((K-1)*nx, length(x));
       B=-1/6*eye((K-1)*nx,(K-1)*nx);
       offset = 2;
       rowOffset =1;

       for i=1:K-1
          A(rowOffset:rowOffset+nx-1,offset:offset+nx-1) = -1*eye(nx);
          A(rowOffset:rowOffset+nx-1,offset+blkSz:offset+blkSz+nx-1) = eye(nx);
          offset = offset+blkSz;
          rowOffset =rowOffset+nx;
       end
       
    end

   
    numConst = 1; 
    numEqConst = nx*(K-1);
    eqconst = zeros(numEqConst,1);
    
    %DCeq = zeros(numEqConst, length(x)); 
    %DCIn = zeros(numConst*(K), length(x)); 
    

    
    cin = zeros(numConst*(K),1);
    deltaT = x(1);

    offset = 2;
    qOffset=1;
    rowIneqOffset=1;
    %Define state-based defects. q follows Bettes pg 142
    q = zeros((K-1)*nx,1);
    fnow = diffDriveDynamics(x(offset:offset+nx-1), x(offset+nx:offset+nx+nu-1));
    for i=1:K-1
        xNow = x(offset:offset+nx-1);
        uNow = x(offset+nx:offset+nx+nu-1);
        uMid = x(offset+nx+nu:offset+blkSz-1);
        xNext = x(offset+blkSz:offset+blkSz+nx-1);
        uNext = x(offset+blkSz+nx:offset+blkSz+nx+nu-1);
        
        
        %We can make the foloowing line faster by doing
        [q(qOffset:qOffset+nx-1), flast] =  simpsonUpdate(xNow, xNext, uNow, uMid, uNext, deltaT, fnow);
        fnow = flast; %Save the evaluation of dynamics function
        qOffset = qOffset+nx;
       
        cin(rowIneqOffset) = calcDistConst(xNow);
        %DCIn(rowIneqOffset,:) = popSparseInJac(i,DCIn(rowIneqOffset,:), xNow);
        rowIneqOffset = rowIneqOffset+1;
        
        offset = offset+blkSz;
        
    end
    
      %%Populate sparse jacobian matrix for differential drive
    ceq = A*x+B*q;
    %DCeq = popSparseJacSimpson(q, A, B, x);

    

    
return
    
    
    
    %%%%%%%Begin Helper Functions%%%%%%%%%%%
    
    %This function populates the sparse inequality jacobian using finite
    %differences
    function rows = popSparseInJac(i,DinRows, xNow)
        rows = DinRows;
        Idt = eye(nx);
        for j=1:nx
            vectDelt = finDiffDelta*Idt(:,j);
            cDeltPlus = calcDistConst(xNow + vectDelt)';
            cDeltMinus = calcDistConst(xNow - vectDelt)';
            rows(:,1+j+(i-1)*(nx+nu))= (cDeltPlus-cDeltMinus)/(2*finDiffDelta);
        end

    end
    
    %This function populates rows of the sparse equality jacobian analytically
    function jac = popSparseJacSimpson(q,A,B,x)
        jac = zeroes(nx*(K-1),length(x));
        
        
        jac= A+B*jac;
    end


    function distConst = calcDistConst(xNow)
        
        distConst = -obstDistance(xNow(1:2), [], safetyDistance); 
        
    end


   
end



