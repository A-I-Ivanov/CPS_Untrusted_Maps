%%%%Written by Alexander I. Ivanov - 2017%%%%
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

a1 = velData(1, 1:2);
a2 = velData(2, 1:2);

q1 = velData(1, 3:4);
q2 = velData(2, 3:4);

a = a1*(1-scaling)+a2*scaling;
q = q1*(1-scaling)+q2*scaling;

Q = diag(q);

    %a = [0.1978 0];
    %Q = [0.0691 0; 0 0.0351];

end