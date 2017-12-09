function Sigmas = PreComputeSigmas(Sigma0, SigmaV, SigmaW, A, H, K)
    nx = size(SigmaW,1);
    ny = size(SigmaV,1);
    Sigmas = cell(K,1);
    SigmaK = Sigma0;
    SigmaBar = A*SigmaK*A' + SigmaW;
    for i =1:K
        Sigmas{i} = SigmaK;
        Sk = SigmaV + H*SigmaBar*H';
        Kk = SigmaBar*H'*inv(Sk);
        SigmaK = SigmaBar - Kk*Sk*Kk';
        SigmaBar = A*SigmaK*A' + SigmaW;
    end
end