function y_return = simu_data_generate_random_approx(param, ini, totday, startday, dt)

% function to simulate the epidemic data
% the generation of the data does not follow the exactly same probability
% law as the stochastic process described in the paper, while it is an
% approximation of the ground truth stochastic process aiming to accelerate
% the data generation process

%--------------------------------------------------------------------------
% Inputs:
% param: param(1: n)        : H(0) in the n regions
%        param((n+1):(2*n)) : E(0) in the n regions
%        param((2*n+1): end): transmission parameters in the n regions
%    * remark: We assume that R(0) = 0

% ini: ini.totPop  : the initial total population N(0) in the n regions
%      ini.traMat  : the n * n transportation matrix, whose (i,j)-th entry
%                    denotes the population travelling from region i to region j. 
%      ini.gamma   : the inverse of average time of individuals getting
%                    hospitalized
%      ini.delta   : the inverse of average time of individuals getting
%                    removed (recovered or deceased)
%    * remark: For the simulated data, we assume that the transportation
%              volumes keep constant over time.

% totday: total days that the epidemic data are generated

% startday: the starting day of storing the results

% dt: the time step of moving forward

%--------------------------------------------------------------------------
% Outputs:
% yC: the accumulated confirmed cases
% yCR: the accumulated removed cases
% yE: the population in state E
% yH: the population in state H
% yR: the population in state R

%    * remark: the output are collected on a daily basis
%--------------------------------------------------------------------------


n      = length(ini.totPop);
N_days = totday;

yC = zeros(N_days, n);
yE = zeros(N_days, n);
yCR = zeros(N_days, n);
yH = zeros(N_days, n);
yR = zeros(N_days, n);

H  = round(param(1: n));
E  = round(param((n+1):(2*n)));


tot = round(ini.totPop);
S = tot - H - E;
R = zeros(1, n);

traMat    = ini.traMat; 
lambdavec = param((2*n+1): end); 

gamma  = ini.gamma; % value
delta  = ini.delta; % value


day  = 1; % migration once a day after the transmission
deltaC_temp  = zeros(1, n);
deltaCR_temp = zeros(1, n);

traMat_now = traMat;
rowsum_now  = sum(traMat_now, 2)';


for t_step=dt:dt:N_days
        
    r_mat = dt * [lambdavec .* S .* E ./ tot; delta * E; gamma .* H];
    
    r_mat_sample = poissrnd(r_mat);
    
    % transmission
    
    
    S = S - r_mat_sample(1, :);
    E = E + r_mat_sample(1, :) - r_mat_sample(2, :);
    H = H + r_mat_sample(2, :) - r_mat_sample(3, :);
    R = R + r_mat_sample(3, :);
    
    % migration
        
    S = S + (S ./ (tot-H) * traMat_now - rowsum_now .* S ./ (tot-H)) * dt;
    E = E + (E ./ (tot-H) * traMat_now - rowsum_now .* E ./ (tot-H)) * dt;
    R = R + (R ./ (tot-H) * traMat_now - rowsum_now .* R ./ (tot-H)) * dt;
        
    S_idx = (S < 1);
    I_idx = (H < 1);
    E_idx = (E < 1);
    
    S(S_idx) = poissrnd(S(S_idx));
    H(I_idx) = poissrnd(H(I_idx));
    E(E_idx) = poissrnd(E(E_idx));
    
    
    deltatot = zeros(1, n);
    tot = tot + deltatot;
    

    %R = tot - (I + E + S);
    
    deltaC_temp  = deltaC_temp + r_mat_sample(2, :);
    deltaCR_temp = deltaCR_temp + r_mat_sample(3, :);
    
    if t_step >= day
        
        yC(day, :)  = deltaC_temp;
        yCR(day, :) = deltaCR_temp;
        yE(day, :)  = E;
        yR(day, :)  = R;
        
        if day == 1
            yH(day, :) = deltaC_temp;            
        else
            yH(day, :) = yH(day-1, :) + deltaC_temp;
        end
        
        if day ~= N_days
            day = day + 1;
            deltaC_temp  = zeros(1, n);
            deltaCR_temp = zeros(1, n);
        end
        
        
        
    end
    
end



yC  = yC(startday: end, :);
yE  = yE(startday: end, :);
yCR = yCR(startday: end, :);
yH  = yH(startday: end, :);
yR  = yR(startday: end, :);

y_return = [yC, yCR, yE, yH, yR];