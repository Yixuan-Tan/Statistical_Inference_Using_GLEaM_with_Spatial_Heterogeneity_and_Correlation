function [newly_conf_inf, newly_recov_inf, accu_conf_inf, I,E,R] = real_world_data_generate_determ(param, ini, dt, startday,endday)

% generate predicted trajectories for the countries based on the estimated
% parameters `param`

% param: estimated parameters:
% param(1:n)           : I0
% param((n+1):(2*n))   : E0
% param((2*n+1):(3*n)) : lambda_before
% param((3*n+1):(4*n)) : lambda_later

% ini                  : initial values that are prefixed, delta, gamma, etc.
% dt                   : \Delta_t
% period_ind           : should be 1/0, if 1, then this period is from day
% 1:infer_days, otherwise it is from day (infer_days+1):end

% change_days would be used only if period_ind = 1

N_days    = endday;
n_country = ini.n_country;

newly_conf_inf  = zeros(N_days, n_country);
newly_recov_inf = zeros(N_days, n_country);
accu_conf_inf   = zeros(N_days, n_country);

delta               = ini.delta;
gamma               = ini.gamma;
lambda_before       = param((2*n_country+1):(3*n_country));
lambda_later        = param((3*n_country+1):(4*n_country));

I   = param(1:n_country);
E   = param((n_country+1):(2*n_country));
tot = ini.totPop;
S   = tot - I - E - ini.R0;
R   = ini.R0;

change_day        = ini.change_day;
change_day_unique = unique(change_day); % unique change_day, sorted in ascending order
num_change_day    = length(change_day_unique);
change_day_ind    = cell(num_change_day, 1);
for i = 1: num_change_day
    change_day_ind{i} = find(change_day == change_day_unique(i));
end

change_day_temp_ind = 1;

day_now      = 1;
deltaC_temp  = zeros(1, n_country); % newly confirmed cases in [t, t+dt]
deltaCR_temp = zeros(1, n_country); % newly recovered cases in [t, t+dt]

lambdavec    = lambda_before;


for t_step = dt: dt: N_days
    
    
    if (change_day_temp_ind <= num_change_day) && (t_step >= change_day_unique(change_day_temp_ind))
        lambdavec(change_day_ind{change_day_temp_ind}) = lambda_later(change_day_ind{change_day_temp_ind});
        change_day_temp_ind = change_day_temp_ind + 1;
    end
    
    
    r_mat = dt * [lambdavec .* S .* E ./ tot; delta * E; gamma * I];
    
    S = S - r_mat(1, :);
    E = E + r_mat(1, :) - r_mat(2, :);
    I = I + r_mat(2, :) - r_mat(3, :);
    R = R + r_mat(3, :);
    
    deltaC_temp = deltaC_temp + r_mat(2, :);
    deltaCR_temp = deltaCR_temp + r_mat(3, :);
    
    if t_step >= day_now
        newly_conf_inf(day_now, :)  = deltaC_temp;
        newly_recov_inf(day_now, :) = deltaCR_temp;
        if day_now == 1
            accu_conf_inf(day_now, :)= ini.I_accu + deltaC_temp;
        else
            accu_conf_inf(day_now, :)= accu_conf_inf(day_now-1, :) + deltaC_temp;
        end
        day_now                     = day_now + 1;
        deltaC_temp                 = zeros(1, n_country); % set to 0
        deltaCR_temp                = zeros(1, n_country);
    end
end


newly_conf_inf = newly_conf_inf(startday:endday, :);
newly_recov_inf = newly_recov_inf(startday:endday, :);
accu_conf_inf = accu_conf_inf(startday:endday, :);


end
















