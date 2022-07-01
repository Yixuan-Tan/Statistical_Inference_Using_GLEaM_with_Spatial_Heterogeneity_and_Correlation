function loss = fmin_simu_data(x, ytrue, ini, loss_type, mu, beta, dt, predday, startday, sigma)

% function to compute the loss of the given parameters x
% this function is for the models that allow the parameter heterogeneity
%   * remark: loss is -log p(\theta|Y)

%--------------------------------------------------------------------------
% Inputs:
% x: x(1: n)        : H(0) in the n regions
%    x((n+1):(2*n)) : E(0) in the n regions
%    x((2*n+1): end): transmission parameters in the n regions
%   * remark: We assume that R(0) = 0

% ytrue: the observed newly confirmed cases & newly removed cases

% ini: ini.totPop  : the initial total population N(0) in the n regions
%      ini.traMat  : the n * n transportation matrix, whose (i,j)-th entry
%                    denotes the population travelling from region i to region j. 
%      ini.gamma   : the inverse of average time of individuals getting
%                    hospitalized
%      ini.delta   : the inverse of average time of individuals getting
%                    removed (recovered or deceased)
%    * remark: For the simulated data, we assume that the transportation
%              volumes keep constant over time.

% loss_type: either 'pois' or 'norm'. 'pois' means that the observed data
%            follow a Poisson distribution whose mean is the deterministic 
%            trajectory determined by the parameters. 'norm' means that the
%            observed data follow a multivariate normal distribution
%            instead.

% mu: the parameter of the graph Laplacian penalty

% beta: the parameter taking values in (0, 1) that reduces the correlation 
%       between inter-group regions, see details in the paper

% dt: timestep size for the forward Euler method to integrate the ODE
%     system

% predday: the total days that data are to be fitted

% startday: the starting day that data are to be fitted

% sigma: the factor for regularization on l2 norm of the vector(s) of 
%        transmission parameter. sigma = 1e-6 for experiments in the paper

%--------------------------------------------------------------------------
% Outputs:
% loss: loss of the current parameters given the observed data. The best
%       parameters are chosen to minimize the loss
%--------------------------------------------------------------------------


n      = length(ini.totPop);
N_days = predday; 

yC  = zeros(N_days, n);
yCR = zeros(N_days, n);

H  = x(1: n);
E  = x((n+1):(2*n));

tot = round(ini.totPop);
S = tot - H - E;
R = zeros(1, n);

traMat    = ini.traMat; % value
part      = ini.part;
partNum   = length(part);

lambdavec = x((2*n+1): end); % before

gamma  = ini.gamma; % value
delta  = ini.delta; % value
%--------------------------------------------------------
if max(traMat,[],'all') > 0
    W_reg  = traMat / max(traMat,[],'all');
elseif max(traMat,[],'all') == 0
    W_reg  = traMat;
end
mu1 = mu;
mu0 = mu*beta;
for i = 1: partNum
    for j = 1: partNum
        if j == i
            W_reg(part{i}, part{j}) = mu1 * W_reg(part{i}, part{j}); % within block, mu1 large
        else
            W_reg(part{i}, part{j}) = mu0 * W_reg(part{i}, part{j}); % between block, mu0 small
        end
    end
end
Lap_reg = diag(sum(W_reg)) - W_reg;


day = 1;
deltaC_temp  = zeros(1, n);
deltaCR_temp = zeros(1, n);

traMat_now = traMat;
rowsum_now  = sum(traMat_now, 2)';

%--------------------------------------------------------
for t_step=dt:dt:N_days
       
    r_mat = dt * [lambdavec .* S .* E ./ tot; delta * E; gamma .* H];
    
    S = S - r_mat(1, :) + (S ./ (tot-H) * traMat_now - rowsum_now .* S ./ (tot-H)) * dt;
    E = E + r_mat(1, :) - r_mat(2, :) + (E ./ (tot-H) * traMat_now - rowsum_now .* E ./ (tot-H)) * dt;
    H = H + r_mat(2, :) - r_mat(3, :);
    
    %deltatot = zeros(1, n);
    %tot = tot + deltatot;
    
    %R = tot - (I + E + S);
    
    deltaC_temp  = deltaC_temp + r_mat(2, :);
    deltaCR_temp = deltaCR_temp + r_mat(3, :);
    
    if t_step >= day
        
        yC(day, :)  = deltaC_temp;
        yCR(day, :) = deltaCR_temp;
        
        if day ~= N_days
            
            day = day + 1;
            deltaC_temp  = zeros(1, n);
            deltaCR_temp = zeros(1, n);
        end
    end
    
end

ytrue = ytrue(1:(predday-startday+1), 1:n);
%--------------------------------------------------------

if strcmp(loss_type, 'pois')
    y    = max(yC(startday: predday, :), 1);
    loss = -sum(ytrue(:) .* log(y(:)) - y(:)) + lambdavec * Lap_reg * lambdavec' + sigma * norm(lambdavec)^2;

elseif strcmp(loss_type, 'norm')
    y    = max([yC(startday: predday, :), yCR(startday: predday, :)], 1);
    %    loss = norm(y(:) - ytrue(:), 2)^2 / temp;
    loss = sum(diag(((y - ytrue)./(max(ytrue,1)).^(1/2)) * ((y - ytrue)./(max(ytrue,1)).^(1/2))'))...
        + lambdavec * Lap_reg * lambdavec';
end