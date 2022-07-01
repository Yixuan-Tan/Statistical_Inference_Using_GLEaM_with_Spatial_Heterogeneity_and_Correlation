function [newly_conf_inf, newly_recov_inf, accu_conf_inf] = simu_europe_data_generate_random(param, ini, startday, endday)

% generate predicted trajectories for the countries based on the estimated
% (random version)
% parameters `param`

% param: estimated parameters:
% param(1:n)           : I0
% param((n+1):(2*n))   : E0
% param((2*n+1):(3*n)) : lambda_before
% param((3*n+1):(4*n)) : lambda_later

% ini                  : initial values that are prefixed, delta, gamma, etc.
% dt                   : \Delta_t

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

change_day        = ini.change_day;
change_day_unique = unique(change_day); % unique change_day, sorted in ascending order
num_change_day    = length(change_day_unique);
change_day_ind    = cell(num_change_day, 1);
for i = 1: num_change_day
    change_day_ind{i} = find(change_day == change_day_unique(i));
end

change_day_temp_ind = 1;


day_now      = 1;
t_now        = 0;
deltaC_temp  = zeros(1, n_country); % newly confirmed cases in [t, t+dt]
deltaCR_temp = zeros(1, n_country); % newly recovered cases in [t, t+dt]
lambdavec    = lambda_before;


while 1
    
    if (change_day_temp_ind <= num_change_day) && (t_now >= change_day_unique(change_day_temp_ind))
        lambdavec(change_day_ind{change_day_temp_ind}) = lambda_later(change_day_ind{change_day_temp_ind});
        change_day_temp_ind = change_day_temp_ind + 1;
    end
    
    r_mat = [lambdavec .* S .* E ./ tot; delta * E; gamma * I];
    
    delta_t_mat = exprnd(1 ./ r_mat);
    delta_t     = min(min(delta_t_mat));
    
    t_now = t_now + delta_t;
    [ind_row, ind_col] = find(delta_t_mat == delta_t);
    
    if ind_row == 1
        S(ind_col) = S(ind_col) - 1;
        E(ind_col) = E(ind_col) + 1;
    elseif ind_row == 2
        E(ind_col) = E(ind_col) - 1;
        I(ind_col) = I(ind_col) + 1;
        deltaC_temp(ind_col)  = deltaC_temp(ind_col) + 1;
    else
        I(ind_col) = I(ind_col) - 1;
        deltaCR_temp(ind_col) = deltaCR_temp(ind_col) + 1;
    end
    
    if t_now >= day_now
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
    
    if t_now > N_days
        break;
    end

end

newly_conf_inf = newly_conf_inf(startday:endday, :);
newly_recov_inf = newly_recov_inf(startday:endday, :);
accu_conf_inf = accu_conf_inf(startday:endday, :);


end














