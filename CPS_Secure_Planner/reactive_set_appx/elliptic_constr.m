function [cin, constraints] = elliptic_constr(x)
global vectorized 
    
    constraints = [];
    xRep = repmat(x(1:2),length(vectorized),1);
    baseScew = zeros(2);
    baseScew(1,1) = x(3);
    baseScew(2,2) = x(4);
    scewRep= inv(baseScew);
    
  
    
    distances = xRep-vectorized;
    
    for i=1:length(vectorized)
       cin(i) = distances(i,:)*scewRep*distances(i,:)'- 1; 
    end
    
    
    
end
