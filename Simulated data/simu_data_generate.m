function [yC, yCR, yE, yH, yR] = simu_data_generate(param, ini, totday, startday)

% function to simulate the epidemic data
% the generation of the data follow the exactly same probability law as the
% stochastic process described in the paper.
% the code might be slow

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

traVal    = ini.traVal; % value

lambdavec = param((2*n+1): end); % before

gamma  = ini.gamma; % value
delta  = ini.delta; % value


t_temp = 0;
t_day  = 1; % migration once a day after the transmission
deltaC_temp  = zeros(1, n);
deltaCR_temp = zeros(1, n);

traMat_now = traVal * ones(n) - diag(traVal * ones(1, n));
rowsum_now  = sum(traMat_now, 2)';

while 1    
    
    r_mat = [lambdavec .* S .* E ./ tot; delta * E; gamma .* H];
    
    delta_t_mat = exprnd(1 ./ r_mat);
    delta_t     = min(min(delta_t_mat));
    
    t_temp = t_temp + delta_t;
    
    [ind_row, ind_col] = find(delta_t_mat == delta_t);
    
    if ind_row == 1
        S(ind_col) = S(ind_col) - 1;
        E(ind_col) = E(ind_col) + 1;
    elseif ind_row == 2
        E(ind_col) = E(ind_col) - 1;
        H(ind_col) = H(ind_col) + 1;
        deltaC_temp(ind_col) = deltaC_temp(ind_col) + 1;
    else
        H(ind_col) = H(ind_col) - 1;
        R(ind_col) = R(ind_col) + 1;
        deltaCR_temp(ind_col) = deltaCR_temp(ind_col) + 1;
    end
    
    
    
    if t_temp >= t_day
                
        % migration multinomial chain method        
        
        mnrate_mat = traVal * ones(n) ./ (tot' * ones(1, n));
        mnrate_mat = mnrate_mat - diag(sum(mnrate_mat, 2)) + eye(n);
        
        %S_mig_mat = mnrnd(S', mnrate_mat);
        E_mig_mat = mnrnd(E', mnrate_mat);
        H_mig_mat = mnrnd(H', mnrate_mat);
        
        %S = S + sum(S_mig_mat) - sum(S_mig_mat, 2)';
        S = round(S + (S ./ tot * traMat_now - rowsum_now .* S ./ tot));
        E = E + sum(E_mig_mat) - sum(E_mig_mat, 2)';
        H = H + sum(H_mig_mat) - sum(H_mig_mat, 2)';
        R = tot - S - E - H;
                
        
        yC(t_day, :)  = deltaC_temp + sum(H_mig_mat) - diag(H_mig_mat)';
        yE(t_day, :)  = E;
        %yI(t_day, :)  = I;
        yR(t_day, :)  = R;
        yCR(t_day, :) = deltaCR_temp;
        
        if t_day == 1
            yH(t_day, :) = deltaC_temp;
        else
            yH(t_day, :) = yH(t_day-1, :) + deltaC_temp;
        end
            
        
        t_day = t_day + 1;
        deltaC_temp  = zeros(1, n);
        deltaCR_temp = zeros(1, n);
        
        
    end
    
    
    if t_temp > N_days
        break;
    end
        
    
end

yC  = yC(startday: end, :);
yE  = yE(startday: end, :);
yCR = yCR(startday: end, :);
yH  = yH(startday: end, :);
yR  = yR(startday: end, :);