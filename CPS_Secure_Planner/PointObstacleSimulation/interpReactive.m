function [a,Q] = interpReactive(velocity)
persistent velData indexVels pp
%scaling = mod(velocity, .1); %We discretized velocity into .1m/s previously

if isempty(indexVels)
   indexVels = linspace(.1,1.9,19);
   velData = [0,0.000359216275074815,0.000196504917537513,0.000359216275074815;
    0,0.000943716823467047,0.000196504917537513,0.000943716823467047;
    0,0.00254285533931621,0.000196504917537513,0.00254285533931621;
    0,0.00582479718882070,0.000196504917537513,0.00582479718882070;
    0,0.0110712226331336,0.000196504917537513,0.0110712226331336;
    0,0.0198164538931231,0.000196504917537513,0.0198164538931231;
    0,0.0326262945366763,0.000196504917537513,0.0326262945366763;
    0,0.0450236300128605,0.000196504917537513,0.0508902388392306;
    0,0.0786498261593462,0.00213969775637027,0.0776507759547518;
    0,0.138126372920431,0.00617618268587879,0.114334432291792;
    0,0.220479277216832,0.0494064747829253,0.176807383272560;
    0,0.268057475038376,0.202148193699846,0.219081691624916;
    0,0.334266259375488,0.268700658453067,0.285445284358851;
    0,0.369356152219121,0.343288038586703,0.324166618236438;
    0,0.406139231222650,0.412846065905222,0.367405708963829;
    0,0.434009312481233,0.490523301399586,0.401968738236303;
    0,0.464297203401451,0.570909961570336,0.441290958251584;
    0,0.496269572116941,0.643095906324830,0.484791498039827;
    0,0.5,0.65,0.5]; 

    pAy = zeros(1,4);
    pAx = sigm_fit(indexVels,  velData(:,2)', [.000359216275074815, NaN,NaN,NaN],[],0);
    pQx = sigm_fit(indexVels,  velData(:,3)', [.000196504917537513,  NaN,NaN,NaN],[],0);
    pQy = sigm_fit(indexVels,  velData(:,4)', [.000359216275074815, NaN,NaN,NaN],[],0);
    
    pAy(1) = max(pAy(1),0);
    pAx(1) = max(pAx(1),0);
    pQy(1) = max(pQy(1),0);
    pQx(1) = max(pQx(1),0);
    pp = [pAx; pAy; pQx; pQy];
    
 

end


%ix = find(indexVels>=velocity,1);
%Data from simulation. Should read from file instead 

%velData = [0.0786498261593462, 0,0.0213969775637027,0.0109216275074815 ;
   % 0.496269572116941,0,0.0643095906324830, 0.0484791498039827];

%debug 
%if velocity<0
%    ix =-1 %error out
%end
%if velocity>1.9
%   ix = 19;
%end
   
%if ix==1
%    fin = 1;
%    start =1;
%else if ix==19
%        start = ix;
%        fin = ix;
%    else
%        fin = ix;
%        start = ix-1;
%    end
%end




%a = velData(start, 1:2)*scaling+velData(fin, 1:2)*(1-scaling);
%q = velData(start, 3:4)*scaling+velData(fin, 3:4)*(1-scaling);
if(velocity<0)
    velocity =0;
end
if velocity > 1.9
   velocity = 1.9; 
end

%interps = ppval(pp, velocity);
interps(1) = sigmoid(velocity, pp(1,:));
interps(2) = sigmoid(velocity, pp(2,:));
interps(3) = sigmoid(velocity, pp(3,:));
interps(4) = sigmoid(velocity, pp(4,:));
a = interps(1:2)';
Q = diag(interps(3:4));

if (interps(3) <0 || interps(4)<0)
    gay =1;
    
end


end