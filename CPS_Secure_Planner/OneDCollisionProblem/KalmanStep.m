function xhat = KalmanStep(x, u, z, Sigma0, SigmaW, SigmaV)
global Adisc H Bdisc

nx = size(SigmaW,1);
ny = size(SigmaV,1);
SigmaK = Sigma0;
SigmaBar = Adisc*SigmaK*Adisc' + SigmaW;
Sk = SigmaV + H*SigmaBar*H';
Kk = SigmaBar*H'*inv(Sk);
SigmaK = SigmaBar - Kk*Sk*Kk';

xbar = Adisc*x+Bdisc*u;
xhat = xbar + Kk*(z - H*xbar);

end