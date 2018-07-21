%%%%Written by Alexander I. Ivanov - 2017%%%%
function [update, fnext, fmid, xMid] = simpsonUpdate(xNow, xNext,uNow, uMid,uNext, deltaT, fnow)
fnext = diffDriveDynamics(xNext, uNext);
xMid = (xNow+xNext)/2+deltaT/8*(fnow-fnext);
fmid =  diffDriveDynamics(xMid,uMid);

%Note, may have to return the un-scaled update to increase efficiency
update = deltaT*(fnext+4*fmid+fnow);

    %Time doesnt enter explicity into dynamics   
end

