% Computing the Hessian by finite difference methods

LearnTAP;
P = [KTrue; JTrueVec];

N       = length(P);
HMat    = zeros(N,N);
h       = 1e-6; 

for mm = 1:N
    for nn = 1:N
        DelP1 = zeros(N,1); DelP1([mm,nn]) = h;
        DelP2 = zeros(N,1); DelP2(mm) = h; DelP2(nn) = -h;
        DelP3 = zeros(N,1); DelP3(nn) = h; DelP3(mm) = -h;
        DelP4 = zeros(N,1); DelP4([mm,nn]) = -h;
        HMat(mm,nn) = (TAPCostCombined(P+DelP1) - TAPCostCombined(P+DelP2) - TAPCostCombined(P+DelP3) + TAPCostCombined(P+DelP4))/(4*h^2);   
    end
end
