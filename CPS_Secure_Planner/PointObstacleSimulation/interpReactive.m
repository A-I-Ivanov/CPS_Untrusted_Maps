function [a,Q] = interpReactive(velocity)

%scaling = mod(velocity, .1); %We discretized velocity into .1m/s previously
scaling = (velocity)/1.9;
if(scaling<0)
    scaling = 0;
end
%if(scaling>1)
%    scaling =1;
%end

start = floor((velocity-scaling)/.1);

%if(start<1)
%   start
%   velocity
%   scaling
%end
fin = start+1;

if(start<1)
    start=1;
    fin = 1;
end

if(fin>18)
    start = 18;
    fin = 18;
end

%Data from simulation. Should read from file instead 
velData = [0.0786498261593462, 0,0.0213969775637027,0.0109216275074815 ;
    0.496269572116941,0,0.0643095906324830, 0.0484791498039827];

%%debugging
start=1; fin =2;

%a1 = velData(start, 1:2);
%a2 = velData(fin, 1:2);

%q1 = velData(start, 3:4);
%q2 = velData(fin, 3:4);

a = velData(start, 1:2)*(1-scaling)+velData(fin, 1:2)*scaling;
q = velData(start, 3:4)*(1-scaling)+velData(fin, 3:4)*scaling;

Q = diag(q);


end