%%%%Written by Alexander I. Ivanov - 2017%%%%
function [eqconst, DCeq] = basicDynamicsSimpsonEqonly(x)
    finDiffDelta = 1e-6;
    
%These constraints use a Euler approximation to the defect constraints 
%of a differential drive robot. More accurate defects are required 
%for robots with fast or complex dynamics.
global K

persistent A B blkSz nx nu numVar

    if(isempty(A))
        nx =5;
        nu =2;
        numVar = length(x);
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

    
    deltaT = x(1);

    offset = 2;
    qOffset=1;
    rowIneqOffset=1;
    %Define state-based defects. q follows Bettes pg 142
    q = zeros((K-1)*nx,1);
    fNow = diffDriveDynamics(x(offset:offset+nx-1), x(offset+nx:offset+nx+nu-1));
    
    
    for l=1:K-1
        xNow = x(offset:offset+nx-1);
        uNow = x(offset+nx:offset+nx+nu-1);
        uMid = x(offset+nx+nu:offset+blkSz-1);
        xNext = x(offset+blkSz:offset+blkSz+nx-1);
        uNext = x(offset+blkSz+nx:offset+blkSz+nx+nu-1);
        
        
        %We can make the foloowing line faster by doing
        [q(qOffset:qOffset+nx-1), fNext, fMid, xMid] =  simpsonUpdate(xNow, xNext, uNow, uMid, uNext, deltaT, fNow);
        D(qOffset:qOffset+nx-1,:) = popSparseJacSimpson(q(qOffset:qOffset+nx-1),xNow, xMid,xNext, uNow,uMid,uNext, fNow, fNext, offset);
        
        fNow = fNext; %Save the evaluation of dynamics function
        
        qOffset = qOffset+nx;

       
        rowIneqOffset = rowIneqOffset+1;
        
        offset = offset+blkSz;
        
    end
    %cin(rowIneqOffset) = calcDistConst(x(offset+blkSz:offset+blkSz+nx-1));
      %%Populate sparse jacobian matrix for differential drive
    eqconst = A*x+B*q;
    DCeq = A+(B*D);
    DCeq = DCeq';

    

    
return
    
    
    
    
    %This function populates rows of the sparse equality jacobian analytically
    function Drows = popSparseJacSimpson(qk,xNow, xMid,xNext, uNow,uMid,uNext, fNow, fNext,colOffset)
        Drows = zeros(nx, numVar);
        dfMid = diffDriveDynamicsDeriv(xMid, uMid);
        dfNow = diffDriveDynamicsDeriv(xNow, uNow);
        dfNext = diffDriveDynamicsDeriv(xNext, uNext);
        Inx = eye(nx);
        temp1 = [Inx/2 zeros(nx,2)]+deltaT/8*dfNow;
        
        temp2 = [Inx/2 zeros(nx,2)]-deltaT/8*dfNext;
      
        
        chainXdfMid1 = 4*dfMid*[temp1; zeros(2,7)];
        chainXdfMid2 = 4*dfMid*[temp2; zeros(2,7)];
        Drows(:,1) = qk/deltaT+deltaT*(dfMid*[fNow-fNext;0;0;])/2;
        Drows(:,colOffset:colOffset+blkSz+nx+nu-1) = deltaT*[dfNow+chainXdfMid1 4*dfMid(:,nx+1:end) dfNext+chainXdfMid2];

        
    end


    function distConst = calcDistConst(xNow)
        
        distConst = -obstDistance(xNow(1:2), [], safetyDistance); 
        
    end


   
end



