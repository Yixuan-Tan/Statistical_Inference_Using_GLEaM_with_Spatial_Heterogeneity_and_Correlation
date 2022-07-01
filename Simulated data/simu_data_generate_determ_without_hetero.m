function y = simu_data_generate_determ_without_hetero(param, ini, dt, totday, startday, R)


n      = length(ini.totPop);
N_days = totday;

yC = zeros(N_days, n);
yCR = zeros(N_days, n);
yI = zeros(N_days, n);


I  = round(param(1: n));
E  = round(param((n+1):(2*n)));



tot = round(ini.totPop);
S = tot - I - E- R;


traMat    = ini.traMat; % value

lambdavec = param(2*n+1) * ones(1, n); % before

gamma  = ini.gamma; % value
delta  = ini.delta; % value


day = 1;
deltaC_temp  = zeros(1, n);
deltaCR_temp = zeros(1, n);
traMat_now   = traMat;
rowsum_now   = sum(traMat_now, 2)';


for t_step=dt:dt:N_days
    
    r_mat = dt * [lambdavec .* S .* E ./ tot; delta * E; gamma .* I];
    N_active = S+E+R;%S+E+R+I;
    
    S = S - r_mat(1, :) + (S ./ N_active * traMat_now - rowsum_now .* S ./ N_active) * dt;
    E = E + r_mat(1, :) - r_mat(2, :) + (E ./ N_active * traMat_now - rowsum_now .* E ./ N_active) * dt;
    I = I + r_mat(2, :) - r_mat(3, :);
    %I = I + r_mat(2, :) - r_mat(3, :) + (I ./ N_active * traMat_now - rowsum_now .* I ./ N_active) * dt;
    R = R + r_mat(3, :) + (R ./ N_active * traMat_now - rowsum_now .* R ./ N_active) * dt;
    
    deltatot = zeros(1, n);
    tot = tot + deltatot;
    
    %R = tot - (I + E + S);
    
    deltaC_temp  = deltaC_temp + r_mat(2, :);
    %deltaC_temp  = deltaC_temp + r_mat(2, :) + I ./ N_active * traMat_now * dt;
    deltaCR_temp = deltaCR_temp + r_mat(3, :);
    
    if t_step >= day
        
        yC(day, :)  = deltaC_temp;
        yCR(day, :) = deltaCR_temp;
        
        if day == 1
            yI(day, :) = deltaC_temp;
        else
            yI(day, :) = yI(day-1, :) + deltaC_temp;
        end
        
        if day ~= N_days
                        
            day = day + 1;
            deltaC_temp  = zeros(1, n);
            deltaCR_temp = zeros(1, n);
        end
    end
    
end

yC1  = yC(startday: end, :);
yCR1 = yCR(startday: end, :);
yI1  = yI(startday: end, :);

y = [yC1, yCR1, yI1];