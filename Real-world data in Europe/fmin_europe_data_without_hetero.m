function loss = fmin_europe_data_without_hetero(x, newly_conf_true, ini, loss_type, dt, startday,endday,sigma)

% change_day could be different for countries, which should be stored in
% ini

% x: to be predicted
% x(1:n)           : I0
% x((n+1):(2*n))   : E0
% x((2*n+1):(3*n)) : lambda_before
% x((3*n+1):(4*n)) : lambda_later

% newly_conf_true  : N * n, the real newly confirmed cases in the n countries
% ini              : initial values that are prefixed, delta, gamma, etc.
% change_day       : the threshold when lambda would be changed
% loss_type        : could be 'Pois' (Poisson loss) or 'Norm' (Normal loss)
% dt               : \Delta_t

[~, n_country] = size(newly_conf_true);
N_days         = endday;
yCstartday = ini.yCstartday;
newly_conf_est = zeros(N_days, n_country);

delta               = ini.delta;
gamma               = ini.gamma;
lambda_before       = x((2*n_country+1)); % scalar
lambda_later        = x((2*n_country+2)); % scalar


I   = x(1:n_country)*1e6;
E   = x((n_country+1):(2*n_country))*1e6;
tot = ini.totPop; 
S   = tot - I - E - ini.R0;

day_now     = 1;
deltaC_temp = zeros(1, n_country);

change_day        = ini.change_day;
change_day_unique = unique(change_day); % unique change_day, sorted in ascending order
num_change_day    = length(change_day_unique);
change_day_ind    = cell(num_change_day, 1);
for i = 1: num_change_day
    change_day_ind{i} = find(change_day == change_day_unique(i));
end

change_day_temp_ind = 1;

lambdavec   = lambda_before * ones(1, n_country);

for t_step = dt: dt: N_days
%     if t_step >= change_day
%         lambdavec = lambda_later;
%     end
    
    if (change_day_temp_ind <= num_change_day) && (t_step >= change_day_unique(change_day_temp_ind))
        lambdavec(change_day_ind{change_day_temp_ind}) = lambda_later;
        change_day_temp_ind = change_day_temp_ind + 1;
    end
    
    r_mat = dt * [lambdavec .* S .* E ./ tot; delta * E; gamma * I];
    
    S = S - r_mat(1, :); 
    E = E + r_mat(1, :) - r_mat(2, :);
    I = I + r_mat(2, :) - r_mat(3, :);
    
    deltaC_temp = deltaC_temp + r_mat(2, :);
    
    if t_step >= day_now
        newly_conf_est(day_now, :) = deltaC_temp;
        day_now                    = day_now + 1;
        deltaC_temp                = zeros(1, n_country);
    end
end


newly_conf_est  = newly_conf_est(startday:endday, :);
newly_conf_true = newly_conf_true((startday-yCstartday+1):(endday-yCstartday+1), :);


if strcmp(loss_type, 'pois')
    newly_conf_est    = max(newly_conf_est, 1);
    %loss = -sum(newly_conf_true(:) .* log(newly_conf_est(:)) - newly_conf_est(:))+ 1e5*lambda_later * GL_mat * lambda_later' + 1e5*lambda_before * GL_mat * lambda_before'...
    %    + sigma*norm(lambda_later,2)^2 + sigma*norm(lambda_before,2)^2;
    loss = -sum(newly_conf_true(:) .* log(newly_conf_est(:)) - newly_conf_est(:))...
         + sigma*norm(lambda_later,2)^2 + sigma*norm(lambda_before,2)^2;
    
elseif strcmp(loss_type, 'norm') % to be revised
    y    = [yC1, yCR1];
    %    loss = norm(y(:) - ytrue(:), 2)^2 / temp;
    loss = sum(diag(((y - ytrue1)./(max(ytrue1,1)).^(1/2)) * ((y - ytrue1)./(max(ytrue1,1)).^(1/2))')) ;
end



















