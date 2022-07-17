clear all;
load("data_averaged.mat");
whos

% for example, Italy 80, France 60, Norway 75

%% using Germany's data to estimate gamma (recovery rate)

German_ind = map_country_to_ind("Germany");

newly_recov_ave_Germany  = newly_recovered_data_ave(:, German_ind);
accu_inHosp_ave_Germany = accu_confirmed_data_ave(1: end-1, German_ind) - accu_recovered_data_ave(1: end-1, German_ind);

figure(1), clf;
subplot(1,2,1);
plot(1:(tot_days-7), newly_recov_ave_Germany, 'x-')
title('newly confirmed cases in Germany');
subplot(1,2,2);
plot(1:(tot_days-7), accu_inHosp_ave_Germany, 'x-')
title('cases in hospital in Germany');

f_gamma = @(gamma) gamma_fmin(gamma, accu_inHosp_ave_Germany, newly_recov_ave_Germany);
opts_gamma  = optimoptions('fminunc');
opts_gamma.MaxIterations =  2000;
opts_gamma.MaxFunctionEvaluations = 10^6;
opts_gamma.Display = 'off';
gamma_Germ = fminunc(f_gamma, 0.14, opts_gamma) % will be used as the uniform gamma for all European countries


%%

data = mean(newly_confirmed_data_ave(1:100,:),2);


%------------------------------------------------------
[~,change_day] = min(mean(newly_confirmed_data_ave(1:105,:),2))


%% model 3' with different mu


% suppose mu1 = beta mu0 for fixed beta, and we let mu0
% change in a range

beta = 0.1; % beta = mu0 / mu1, which should be a constant smaller than or equal to 1
muvec = [0, 10.^(2:0.1:7.5)]

sigma = 1e-6;
  
%---------------------------------------
% initialize settings

dt         = 0.1;
delta      = 0.14;
n_country  = length(country_names);

ini.delta  = delta;
ini.gamma  = gamma_Germ;
ini.totPop = [5806081, 5521773, 5328212, 8859992, 82979100, 8526932, 60377663, 46733038, 11455358, 66992699, 4857000];
ini.R0     = accu_recovered_data_ave(1,:);
ini.I_accu = accu_confirmed_data_ave(1,:);
ini.n_country = n_country;

group_num = 4;
group_cell = cell(group_num,1);

group_cell{1} = [1;2;3];
group_cell{2} = [4;5;6];
group_cell{3} = [7;8];
group_cell{4} = [9;10;11];

% % group_cell{1} = [2;3]
% % group_cell{2} = [1;4;5]
% % group_cell{3} = [8;9;10;11]
% % group_cell{4} = [6;7]

commu_mat_within_clust = zeros(n_country);
for i = 1:group_num
    ind_temp = group_cell{i};
    commu_mat_within_clust(ind_temp,ind_temp) = ones(length(ind_temp));
end

commu_mat_within_clust = commu_mat_within_clust - eye(n_country);

commu_mat_between_clust = ones(11) - commu_mat_within_clust;
commu_mat_between_clust = commu_mat_between_clust - eye(11);
ini.commu_mat_within_clust  = commu_mat_within_clust;
ini.commu_mat_between_clust = commu_mat_between_clust;

% change_day = 50;
ini.change_day = change_day * ones(n_country, 1); % same for every country

yCstartday = 1;
ini.yCstartday = yCstartday;

%---------------------------------------
infer_days = 102;%100;

traindays_ratio = 0.07;

totDays  = tot_days-7; % after averaging the data
% used to choose mu
trainDays_begin1 = 1;
trainDays_end1   = 90;
validDays_begin1 = trainDays_end1+1;
validDays_end1   = 105;
% validDays_end1   = change_day;

trainDays_begin2 = change_day+1;
trainDays_end2   = change_day+6;
validDays_begin2 = trainDays_end2+1;

% used to compute lambda in the second period after choosing mu
trainDays_begin3 = 1;
trainDays_end3   = 105;
testDays_begin   = trainDays_end3+1;
testDays_end     = totDays;

%
%---------------------------------------
% inputs for fmincon

lb2 = [0.0* ones(1, 2*n_country), 0.08*ones(1,2*n_country)];
ub2 = [0.5* ones(1, n_country), 0.5 * ones(1, n_country), 0.5*ones(1,2*n_country)];
param_pre2 = [0.1 * ones(1, n_country), 0.5 * ones(1, n_country), 0.14 * ones(1,n_country),0.14 * ones(1,n_country)];

opts_wGL = optimoptions('fmincon');
opts_wGL.MaxIterations =  2000;
opts_wGL.MaxFunctionEvaluations = 10^6;
opts_wGL.Display = 'off';%'final-detailed';%'off';%'iter-detailed';
opts_wGL.ConstraintTolerance = 1e-8;


% weights for computing errors
weights2_train = 1/infer_days * ones(infer_days,n_country);
weights2_test  = 1/(tot_days - 7 - infer_days) * ones(tot_days-7-infer_days,n_country);

weights1_train = newly_confirmed_data_ave(1:infer_days,:) ./ (ones(infer_days,1)*accu_confirmed_data_ave(infer_days,:));
weights1_test  = newly_confirmed_data_ave((infer_days+1):end,:) ./ (ones(tot_days - 7 - infer_days,1)...
    *(accu_confirmed_data_ave(end,:) - accu_confirmed_data_ave(infer_days+1,:) ) );

%---------------------------------------
% validation weights

% weight2: simple average

% tot_valid_days = (validDays_end2-validDays_begin2+1) + (validDays_end1-validDays_begin1+1);
% weights2_valid = ones( tot_valid_days , n_country) / tot_valid_days;
tot_valid_days = (validDays_end1-validDays_begin1+1);
weights2_valid = ones( tot_valid_days , n_country) / tot_valid_days;

% weight1: propotional to newly confirmed cases

% weights1_valid  = newly_confirmed_data_ave([validDays_begin1:validDays_end1, validDays_begin2: validDays_end2],:);
% weights1_valid  = weights1_valid ./ (ones(tot_valid_days,1) * sum(weights1_valid,1));
weights1_valid  = newly_confirmed_data_ave(validDays_begin1:validDays_end1,:);
weights1_valid  = weights1_valid ./ (ones(tot_valid_days,1) * sum(weights1_valid,1));

%%
%---------------------------------------
t_list = zeros(length(muvec), 1);

paramsmat1 = zeros(length(muvec), 4 * n_country);

validRelErrmat1 = zeros(length(muvec), 1);
validRelErrmat2 = zeros(length(muvec), 1);
validRelMSEmat1 = zeros(length(muvec), 1);
validRelMSEmat2  = zeros(length(muvec), 1);

traj_newly_conf_valid  = zeros(length(muvec), validDays_end1, n_country);
validRelErrmat1_country = zeros(length(muvec), n_country);

tic
for i = 1:length(muvec)
    
    mu = muvec(i);
    mu1 = mu;
    mu0 = beta * mu1;
    fprintf('the %d-th mu, mu = 10^%.2f\n', i,log10(mu));

    %----------------------------------------------------------------------
    % compute validation error using a part of total training set as
    % "training set" and another part of total training set as "validation
    % set"
    
    % optimization process
    if mu == 0
        f_wGL = @(x) fmin_europe_data(x, newly_confirmed_data_ave, ini, 'pois', dt, mu, beta, trainDays_begin1, trainDays_end1, 0);
    else
        f_wGL = @(x) fmin_europe_data(x, newly_confirmed_data_ave, ini, 'pois', dt, mu, beta, trainDays_begin1, trainDays_end1, sigma);
    end
    
    [params1, ~] = fmincon(f_wGL, param_pre2, [], [], [], [], lb2, ub2, [], opts_wGL);
    params1(1:2*n_country) = 1e6 * params1(1:2*n_country);
    % store the inferred parameters
    paramsmat1(i,:) = params1;
    %
    y_reg = simu_europe_data_generate_determ(params1, ini, dt, trainDays_begin1, validDays_end1);
    y_reg_valid = y_reg(validDays_begin1:validDays_end1,:);
    y_true_valid = newly_confirmed_data_ave(validDays_begin1:validDays_end1,:);
    traj_newly_conf_valid(i,:,:) = y_reg(1:validDays_end1,:);
    
    
    validAbsErr  = abs(y_true_valid - y_reg_valid);
    validRelErr  = validAbsErr ./ y_true_valid;
    
    %
    validRelErr_vec_1 = sum(validRelErr .* weights1_valid, 1);
    validRelErr_vec_2 = sum(validRelErr .* weights2_valid, 1);
    validRelMSE_vec_1 = sqrt(sum(validRelErr.^2 .* weights1_valid, 1));
    validRelMSE_vec_2 = sqrt(sum(validRelErr.^2 .* weights2_valid, 1));

    validRelErrmat2(i) = sum(validRelErr_vec_2);
    validRelErrmat1(i) = sum(validRelErr_vec_1);
    validRelMSEmat2(i) = sum(validRelMSE_vec_2);
    validRelMSEmat1(i) = sum(validRelMSE_vec_1);
    
    validRelErrmat1_country(i,:) = validRelErr_vec_1;
    
    fprintf('total validation err: %f\n\n', validRelErrmat1(i)/n_country);
    
end

toc;

return

%%
%---------------------------------------
% store the inferred parameters, and training &
% testing errors


trainRelErrmat_2periods1 = zeros(length(muvec), 1); 
testRelErrmat1  = zeros(length(muvec), 1);

trainRelErrmat_2periods2 = zeros(length(muvec), 1); 
testRelErrmat2  = zeros(length(muvec), 1);

trainRelMSEmat_2periods1 = zeros(length(muvec), 1); 
testRelMSEmat1  = zeros(length(muvec), 1);

trainRelMSEmat_2periods2 = zeros(length(muvec), 1); 
testRelMSEmat2  = zeros(length(muvec), 1);

paramsmat2 = zeros(length(muvec), 4 * n_country);

traj_newly_conf_train = zeros(length(muvec), trainDays_end3, n_country);
traj_newly_conf_test  = zeros(length(muvec), totDays - trainDays_end3, n_country);


testRelErrmat1_country  = zeros(length(muvec), n_country);


tic
for i = 1:length(muvec)
    
    mu = muvec(i);
    mu1 = mu;
    mu0 = beta * mu1;
    fprintf('the %d-th mu, mu = %6.2e\n', i,mu);
    %----------------------------------------------------------------------
    % compute total training error & testing error
    if mu == 0
        f_wGL = @(x) fmin_europe_data(x, newly_confirmed_data_ave, ini, 'pois', dt, mu, beta, trainDays_begin3, trainDays_end3, 0);
    else
        f_wGL = @(x) fmin_europe_data(x, newly_confirmed_data_ave, ini, 'pois', dt, mu, beta, trainDays_begin3, trainDays_end3, sigma);
    end
   
    [params1, ~] = fmincon(f_wGL, param_pre2, [], [], [], [], lb2, ub2, [], opts_wGL);
    params1(1:2*n_country) = 1e6 * params1(1:2*n_country);
    % store the inferred parameters
    paramsmat2(i,:) = params1;
    y_reg = simu_europe_data_generate_determ(params1, ini, dt, trainDays_begin3, testDays_end);
    
    traj_newly_conf_train(i,:,:) = y_reg(1:trainDays_end3,:);
    traj_newly_conf_test(i,:,:)  = y_reg((trainDays_end3+1):end,:);
    
    trainAbsErr = abs(y_reg(1:(trainDays_end3-trainDays_begin3+1), 1: n_country) - newly_confirmed_data_ave((trainDays_begin3-yCstartday+1):(trainDays_end3-yCstartday+1), :));
    trainRelErr_vec_1 = (sum(trainAbsErr))' ./ (sum(newly_confirmed_data_ave((trainDays_begin3-yCstartday+1):(trainDays_end3-yCstartday+1), :)))' ;
    trainRelErr_vec_2 = sum( trainAbsErr ./ newly_confirmed_data_ave((trainDays_begin3-yCstartday+1):(trainDays_end3-yCstartday+1), :) )' / (trainDays_end3-trainDays_begin3+1);
    trainRelMSE_vec_1 = (sum( (trainAbsErr ./ newly_confirmed_data_ave((trainDays_begin3-yCstartday+1):(trainDays_end3-yCstartday+1), :)).^2 .* newly_confirmed_data_ave((trainDays_begin3-yCstartday+1):(trainDays_end3-yCstartday+1), :)...
        ./ (ones(trainDays_end3-trainDays_begin3+1,1) * (sum(newly_confirmed_data_ave((trainDays_begin3-yCstartday+1):(trainDays_end3-yCstartday+1), :))) )  ) ).^(1/2) ;
    trainRelMSE_vec_2 = (sum( (trainAbsErr ./ newly_confirmed_data_ave((trainDays_begin3-yCstartday+1):(trainDays_end3-yCstartday+1), :)).^2 / (trainDays_end3-trainDays_begin3+1) )).^(1/2) ;

    testAbsErr  = abs(y_reg((testDays_begin-trainDays_begin3+1):(testDays_end-trainDays_begin3+1), 1: n_country) - newly_confirmed_data_ave((testDays_begin-yCstartday+1):(testDays_end-yCstartday+1), :));
    testRelErr_vec_1 = (sum(testAbsErr))' ./ (sum(newly_confirmed_data_ave((testDays_begin-yCstartday+1):(testDays_end-yCstartday+1), :)))' ;
    testRelErr_vec_2 = sum( testAbsErr ./ newly_confirmed_data_ave((testDays_begin-yCstartday+1):(testDays_end-yCstartday+1), :) )' / (testDays_end - testDays_begin + 1);
    testRelMSE_vec_1 = (sum( (testAbsErr ./ newly_confirmed_data_ave((testDays_begin-yCstartday+1):(testDays_end-yCstartday+1), :)).^2 .* newly_confirmed_data_ave((testDays_begin-yCstartday+1):(testDays_end-yCstartday+1), :)...
        ./ (ones(testDays_end-testDays_begin+1,1) * sum(newly_confirmed_data_ave((testDays_begin-yCstartday+1):(testDays_end-yCstartday+1), :))  ) )).^(1/2) ;
    testRelMSE_vec_2 = (sum( (testAbsErr ./  newly_confirmed_data_ave((testDays_begin-yCstartday+1):(testDays_end-yCstartday+1), :)).^2 / (testDays_end-testDays_begin+1) )).^(1/2);
    
    trainRelErrmat_2periods2(i) = sum(trainRelErr_vec_2);
    trainRelErrmat_2periods1(i) = sum(trainRelErr_vec_1);
    trainRelMSEmat_2periods2(i) = sum(trainRelMSE_vec_2);
    trainRelMSEmat_2periods1(i) = sum(trainRelMSE_vec_1);    
    
    testRelErrmat2(i) = sum(testRelErr_vec_2);
    testRelErrmat1(i) = sum(testRelErr_vec_1);
    testRelMSEmat2(i) = sum(testRelMSE_vec_2);
    testRelMSEmat1(i) = sum(testRelMSE_vec_1);
    fprintf('total training err: %f, total testing err: %f\n', trainRelErrmat_2periods1(i)/n_country, testRelErrmat1(i)/n_country);
    fprintf('total training mse: %f, total testing mse: %f\n', trainRelMSEmat_2periods1(i)/n_country, testRelMSEmat1(i)/n_country);

    testRelErrmat1_country(i,:) = testRelErr_vec_1;
end
toc;
%%


[~,ind2_plot] = min(validRelErrmat1)
ind1_plot = ind2_plot-3%12
ind3_plot = ind2_plot+3
% [~,ind3] = min(testRelErrmat1)%42


figure(5), clf;
subplot(1,2,1);
hold on;
plot(log10(muvec(2:end)), validRelErrmat1(2:end)/n_country, 'x-', 'LineWidth', 2, 'MarkerSize', 12);
plot(log10(muvec(2:end)), testRelErrmat1(2:end)/n_country, '--', 'LineWidth', 2, 'Marker', 'o', 'MarkerSize', 7);
plot(log10(muvec(ind1_plot)), validRelErrmat1(ind1_plot)/n_country, 'p', 'MarkerSize', 28, 'MarkerFaceColor', '#EDB120', 'MarkerEdgeColor', '#EDB120');
plot(log10(muvec(ind2_plot)), validRelErrmat1(ind2_plot)/n_country, 'bs', 'MarkerSize', 24, 'MarkerFaceColor', 'b');
plot(log10(muvec(ind3_plot)), validRelErrmat1(ind3_plot)/n_country, 'gd', 'MarkerSize', 24, 'MarkerFaceColor', 'g');
yline(validRelErrmat1(1)/n_country, '-', '$\mu = 0$', 'LineWidth', 2, 'FontSize', 30, 'LabelHorizontalAlignment', 'left', 'Interpreter', 'latex', 'color', '#0072BD');
yline(testRelErrmat1(1)/n_country, '--', '$\mu = 0$', 'LineWidth', 2, 'FontSize', 30, 'LabelHorizontalAlignment', 'left', 'Interpreter', 'latex', 'color', '#D95319');
xlabel('$\log_{10} \mu$', 'Interpreter', 'latex');
ylabel('prediction error', 'Interpreter', 'latex');
set(gca, 'FontSize', 36);
legend('validation error $\textnormal{MAE}^{[\rm Val]}_{(w)}$', ...
    'testing error $\textnormal{MAE}^{[\rm Te]}_{(w)}$', 'FontSize', 32, 'Interpreter', 'latex',...
    'location', 'best')
title('$\textnormal{MAE}_{(w)}$', 'FontSize', 36, 'Interpreter', 'latex');
grid on;
% xlim([0.5 4.3])
ylim([0.29 0.64])
subplot(1,2,2);
hold on;
plot(log10(muvec(2:end)), validRelMSEmat1(2:end)/n_country, 'x-', 'LineWidth', 2, 'MarkerSize', 12);
plot(log10(muvec(2:end)), testRelMSEmat1(2:end)/n_country, '--', 'LineWidth', 2, 'Marker', 'o', 'MarkerSize', 7);
plot(log10(muvec(ind1_plot)), validRelMSEmat1(ind1_plot)/n_country, 'p', 'MarkerSize', 28, 'MarkerFaceColor', '#EDB120', 'MarkerEdgeColor', '#EDB120');
plot(log10(muvec(ind2_plot)), validRelMSEmat1(ind2_plot)/n_country, 'bs', 'MarkerSize', 24, 'MarkerFaceColor', 'b');
plot(log10(muvec(ind3_plot)), validRelMSEmat1(ind3_plot)/n_country, 'gd', 'MarkerSize', 24, 'MarkerFaceColor', 'g');
yline(validRelMSEmat1(1)/n_country, '-', '$\mu = 0$', 'LineWidth', 2, 'FontSize', 30, 'LabelHorizontalAlignment', 'left', 'Interpreter', 'latex', 'color', '#0072BD', 'LabelVerticalAlignment', 'bottom');
yline(testRelMSEmat1(1)/n_country, '--', '$\mu = 0$', 'LineWidth', 2, 'FontSize', 30, 'LabelHorizontalAlignment', 'left', 'Interpreter', 'latex', 'color', '#D95319');
xlabel('$\log_{10} \mu$', 'Interpreter', 'latex');
ylabel('prediction error', 'Interpreter', 'latex');
set(gca, 'FontSize', 36);
legend('validation error $\textnormal{MSE}^{[\rm Val]}_{(w)}$', ...
    'testing error $\textnormal{MSE}^{[\rm Te]}_{(w)}$', 'FontSize', 32, 'Interpreter', 'latex',...
    'location', 'best')
title('$\textnormal{MSE}_{(w)}$', 'FontSize', 36, 'Interpreter', 'latex');
grid on;
% xlim([0.5 4.3])
ylim([0.32 0.69])

sgtitle({'Total weighted relative prediction errors (real-world data in Europe)',''}, 'Fontsize', 42, 'Interpreter', 'latex');%


%--------------------------------------------------------------------------

figure(7), clf;
subplot(1,2,1);
plot(log10(muvec(2:end)), validRelErrmat2(2:end)/n_country, 'x-', 'LineWidth', 2, 'MarkerSize', 12);
hold on;
plot(log10(muvec(2:end)), testRelErrmat2(2:end)/n_country, '--', 'LineWidth', 2, 'Marker', 'o', 'MarkerSize', 7);
plot(log10(muvec(ind1_plot)), validRelErrmat2(ind1_plot)/n_country, 'p', 'MarkerSize', 28, 'MarkerFaceColor', '#EDB120', 'MarkerEdgeColor', '#EDB120');
plot(log10(muvec(ind2_plot)), validRelErrmat2(ind2_plot)/n_country, 'bs', 'MarkerSize', 24, 'MarkerFaceColor', 'b');
plot(log10(muvec(ind3_plot)), validRelErrmat2(ind3_plot)/n_country, 'gd', 'MarkerSize', 24, 'MarkerFaceColor', 'g');
yline(validRelErrmat2(1)/n_country, '-', '$\mu = 0$', 'LineWidth', 2, 'FontSize', 30, 'LabelHorizontalAlignment', 'left', 'Interpreter', 'latex', 'color', '#0072BD', 'LabelVerticalAlignment', 'bottom');
yline(testRelErrmat2(1)/n_country, '--', '$\mu = 0$', 'LineWidth', 2, 'FontSize', 30, 'LabelHorizontalAlignment', 'left', 'Interpreter', 'latex', 'color', '#D95319');
xlabel('$\log_{10} \mu$', 'Interpreter', 'latex');
ylabel('prediction error', 'Interpreter', 'latex');
set(gca, 'FontSize', 36);
legend('validation error $\textnormal{MAE}^{[\rm Val]}_{(s)}$', ...
    'testing error $\textnormal{MAE}^{[\rm Te]}_{(s)}$', 'FontSize', 32, 'Interpreter', 'latex',...
    'location', 'best')
title('$\textnormal{MAE}_{(s)}$', 'FontSize', 36, 'Interpreter', 'latex');
grid on;
% xlim([0.5 4.3])
ylim([0.29 0.65])
subplot(1,2,2);
plot(log10(muvec(2:end)), validRelMSEmat2(2:end)/n_country, 'x-', 'LineWidth', 2, 'MarkerSize', 12);
hold on;
plot(log10(muvec(2:end)), testRelMSEmat2(2:end)/n_country, '--', 'LineWidth', 2, 'Marker', 'o', 'MarkerSize', 7);
plot(log10(muvec(ind1_plot)), validRelMSEmat2(ind1_plot)/n_country, 'p', 'MarkerSize', 28, 'MarkerFaceColor', '#EDB120', 'MarkerEdgeColor', '#EDB120');
plot(log10(muvec(ind2_plot)), validRelMSEmat2(ind2_plot)/n_country, 'bs', 'MarkerSize', 24, 'MarkerFaceColor', 'b');
plot(log10(muvec(ind3_plot)), validRelMSEmat2(ind3_plot)/n_country, 'gd', 'MarkerSize', 24, 'MarkerFaceColor', 'g');
yline(validRelMSEmat2(1)/n_country, '-', '$\mu = 0$', 'LineWidth', 2, 'FontSize', 30, 'LabelHorizontalAlignment', 'left', 'Interpreter', 'latex', 'color', '#0072BD');
yline(testRelMSEmat2(1)/n_country, '--', '$\mu = 0$', 'LineWidth', 2, 'FontSize', 30, 'LabelHorizontalAlignment', 'left', 'Interpreter', 'latex', 'color', '#D95319');
xlabel('$\log_{10} \mu$', 'Interpreter', 'latex');
ylabel('prediction error', 'Interpreter', 'latex');
set(gca, 'FontSize', 36);
legend('validation error $\textnormal{MSE}^{[\rm Val]}_{(s)}$', ...
    'testing error $\textnormal{MSE}^{[\rm Te]}_{(s)}$', 'FontSize', 32, 'Interpreter', 'latex',...
    'location', 'best')
title('$\textnormal{MSE}_{(s)}$', 'FontSize', 36, 'Interpreter', 'latex');
grid on;
% xlim([0.5 4.3])
ylim([0.32 0.71])

sgtitle({'Total simply averaged relative prediction errors (real-world data in Europe)',''}, 'Fontsize', 42, 'Interpreter', 'latex');


%%

fprintf('Print errors:\n\n')

for i = 1:length(muvec)
    fprintf('the %d-th mu0 = 10^%.1f:\n', i,log10(muvec(i)));
%     fprintf('total weighted rela training error     : %5.3f, total weighted rela testing error      : %5.3f\n',...
%         trainRelErrmat_2periods1(i), testRelErrmat1(i));  
%     fprintf('total weighted rela training error(MSE): %5.3f, total weighted rela testing error (MSE): %5.3f\n\n',...
%         trainRelMSEmat_2periods1(i), testRelMSEmat1(i));
    fprintf('total weighted rela training error     : %5.3f, total weighted rela testing error      : %5.3f\n',...
        trainRelErrmat_2periods1(i)/n_country, testRelErrmat1(i)/n_country);  
    fprintf('total weighted rela training error(MSE): %5.3f, total weighted rela testing error (MSE): %5.3f\n\n',...
        trainRelMSEmat_2periods1(i)/n_country, testRelMSEmat1(i)/n_country);

end

%% errors for model 1'


f_woGL = @(x) fmin_europe_data_without_hetero(x, newly_confirmed_data_ave, ini, 'pois', dt, trainDays_begin3, trainDays_end3, 0);
opts_wGL = optimoptions('fmincon');
opts_wGL.MaxIterations =  2000;
opts_wGL.MaxFunctionEvaluations = 10^6;
opts_wGL.Display = 'off';%'off';%'iter-detailed';

% optimization

% lb = zeros(1, 2*n_country+2);
% ub = [5e5 * ones(1, 2*n_country), ones(1,2)];
% param_pre = [1000 * ones(1, n_country), 500 * ones(1, n_country), 0.14,0.14];

lb = zeros(1, 2*n_country+2);
ub = [5 * ones(1, 2*n_country), ones(1,2)];
param_pre = [0.1 * ones(1, n_country), 0.5 * ones(1, n_country), 0.14,0.14];

tic;
[params1, ~] = fmincon(f_woGL, param_pre, [], [], [], [], lb, ub, [], opts_wGL);
params1(1:2*n_country) = 1e5 * params1(1:2*n_country);
%print the results

separate_ind = [3, 6, 7, 8];

fprintf('lambda_before: %6.4f\n:', params1(2*n_country+1));
fprintf('lambda_later : %6.4f\n:', params1(2*n_country+2));

for i = 1: n_country
    fprintf("%11s I_0: %6d, E_0: %6d \n", ...
        country_names(i), round(params1(i)), round(params1(n_country+i)));
    if ismember(i, separate_ind)
        fprintf('\n');
    end
end

toc;
param2 = zeros(1,4*n_country);
param2(1:2*n_country) = params1(1:2*n_country);
param2((2*n_country+1):3*n_country) = params1(1+2*n_country)*ones(1,n_country);
param2((3*n_country+1):4*n_country) = params1(2+2*n_country)*ones(1,n_country);
params1 = param2;
y_reg = simu_europe_data_generate_determ(params1, ini, dt, trainDays_begin3, testDays_end);
y_reg_woHwoM = y_reg;

trainAbsErr = abs(y_reg(1:(trainDays_end3-trainDays_begin3+1), 1: n_country) - newly_confirmed_data_ave((trainDays_begin3-yCstartday+1):(trainDays_end3-yCstartday+1), :));
trainRelErr_vec_1 = (sum(trainAbsErr))' ./ (sum(newly_confirmed_data_ave((trainDays_begin3-yCstartday+1):(trainDays_end3-yCstartday+1), :)))' ;
trainRelErr_vec_2 = sum( trainAbsErr ./ newly_confirmed_data_ave((trainDays_begin3-yCstartday+1):(trainDays_end3-yCstartday+1), :) )' / (trainDays_end3-trainDays_begin3+1);
trainRelMSE_vec_1 = (sum( (trainAbsErr ./ newly_confirmed_data_ave((trainDays_begin3-yCstartday+1):(trainDays_end3-yCstartday+1), :)).^2 .* newly_confirmed_data_ave((trainDays_begin3-yCstartday+1):(trainDays_end3-yCstartday+1), :)...
    ./ (ones(trainDays_end3-trainDays_begin3+1,1) * (sum(newly_confirmed_data_ave((trainDays_begin3-yCstartday+1):(trainDays_end3-yCstartday+1), :))) )  ) ).^(1/2) ;
trainRelMSE_vec_2 = (sum( (trainAbsErr ./ newly_confirmed_data_ave((trainDays_begin3-yCstartday+1):(trainDays_end3-yCstartday+1), :)).^2 / (trainDays_end3-trainDays_begin3+1) )).^(1/2) ;

testAbsErr  = abs(y_reg((testDays_begin-trainDays_begin3+1):(testDays_end-trainDays_begin3+1), 1: n_country) - newly_confirmed_data_ave((testDays_begin-yCstartday+1):(testDays_end-yCstartday+1), :));
testRelErr_vec_1 = (sum(testAbsErr))' ./ (sum(newly_confirmed_data_ave((testDays_begin-yCstartday+1):(testDays_end-yCstartday+1), :)))' ;
testRelErr_vec_2 = sum( testAbsErr ./ newly_confirmed_data_ave((testDays_begin-yCstartday+1):(testDays_end-yCstartday+1), :) )' / (testDays_end - testDays_begin + 1);
testRelMSE_vec_1 = (sum( (testAbsErr ./ newly_confirmed_data_ave((testDays_begin-yCstartday+1):(testDays_end-yCstartday+1), :)).^2 .* newly_confirmed_data_ave((testDays_begin-yCstartday+1):(testDays_end-yCstartday+1), :)...
    ./ (ones(testDays_end-testDays_begin+1,1) * sum(newly_confirmed_data_ave((testDays_begin-yCstartday+1):(testDays_end-yCstartday+1), :))  ) )).^(1/2) ;
testRelMSE_vec_2 = (sum( (testAbsErr ./  newly_confirmed_data_ave((testDays_begin-yCstartday+1):(testDays_end-yCstartday+1), :)).^2 / (testDays_end-testDays_begin+1) )).^(1/2);


trainRelErrmat_2periods2_woHwoM = sum(trainRelErr_vec_2);
trainRelErrmat_2periods1_woHwoM = sum(trainRelErr_vec_1);
trainRelMSEmat_2periods2_woHwoM = sum(trainRelMSE_vec_2);
trainRelMSEmat_2periods1_woHwoM = sum(trainRelMSE_vec_1);    

testRelErrmat2_woHwoM = sum(testRelErr_vec_2);
testRelErrmat1_woHwoM = sum(testRelErr_vec_1);
testRelMSEmat2_woHwoM = sum(testRelMSE_vec_2);
testRelMSEmat1_woHwoM = sum(testRelMSE_vec_1);


%
fprintf('\n\n');
% fprintf('total training relative error: %7.3f\n\n', trainRelErrmat_2periods1_woHwoM);
% fprintf('total testing  relative error: %7.3f\n\n', testRelErrmat1_woHwoM);
% 
% 
% fprintf('total training relative error (MSE): %7.3f\n\n', trainRelMSEmat_2periods1_woHwoM);
% fprintf('total testing  relative error (MSE): %7.3f\n\n', testRelMSEmat1_woHwoM);
fprintf('total training relative error: %7.3f\n\n', trainRelErrmat_2periods1_woHwoM/n_country);
fprintf('total testing  relative error: %7.3f\n\n', testRelErrmat1_woHwoM/n_country);


fprintf('total training relative error (MSE): %7.3f\n\n', trainRelMSEmat_2periods1_woHwoM/n_country);
fprintf('total testing  relative error (MSE): %7.3f\n\n', testRelMSEmat1_woHwoM/n_country);





%% random trajectories sampled from one theta
rng(1218); 

muind = ind2_plot;
muvec(muind)
numsamp    = 1e2;
trajec_mat_one_Theta = zeros(numsamp, tot_days-7, 2*n_country);
nan_count  = 0;
ini.I_accu = accu_confirmed_data_ave(1,:);
params_temp1 = paramsmat2(muind,:);

tic;
for i = 1: numsamp
    [traject_temp_newly, ~, traject_temp_accu] = simu_europe_data_generate_random(params_temp1, ini, trainDays_begin1, totDays);
    trajec_mat_one_Theta(i, :, :) = [traject_temp_newly, traject_temp_accu];
    if max(traject_temp_newly(:, 1: n_country)) == 0
        nan_count = nan_count + 1;
        fprintf('\n the %d the iteration, the trajectory hits zero\n', i);
    end
    
    if mod(i, 20) == 0
        fprintf('i=%d\n',i);
    end
end
toc;
ini.N_days = infer_days; % reset N_days for further inference

%% MCMC

%rng('shuffle');
rng(1218);

% muind = 21;
muind = ind2_plot;
mu   = muvec(muind);

params1        = paramsmat2(muind,:);
params1(1:2*n_country) = 1e-6 * params1(1:2*n_country);
proposal_dist  = 1e-8 * ones(1, length(params1));
proposal_dist(1:(1+2*n_country)) = zeros(1, 1+2*n_country);
iterations    = 5e5; % total number of iterations

numtrasamp    = 1e2; % total number of sampled trajectories
step          = iterations / numtrasamp
%---------------------------------------------------------
% store history
params_post  = zeros(iterations, length(params1));
traject_post = zeros(numtrasamp, tot_days-7, 2*n_country);
loss_post    = zeros(iterations, 1);
Avec = zeros(iterations, 1);

%---------------------------------------------------------

params_temp1 = params1; % starting point
sigma        = 1;
ini.I_accu   = accu_confirmed_data_ave(1,:);
%---------------------------------------------------------


tic;
for i = 1: iterations
    params_post(i, :)     = params_temp1;
    
    loss_post(i) = fmin_europe_data(params_temp1, newly_confirmed_data_ave, ini, 'pois', dt, mu, beta, trainDays_begin3, trainDays_end3, sigma);
    
    if mod(i, step) == 0
        ini.N_days = tot_days - 7;
        params_temp1_samp = params_temp1;
        params_temp1_samp(1:2*n_country) = 1e6 * params_temp1_samp(1:2*n_country);
        [traject_temp_newly, ~, traject_temp_accu] = simu_europe_data_generate_random(params_temp1_samp, ini, trainDays_begin1, totDays);
        traject_post(i/step, :, :) = [traject_temp_newly, traject_temp_accu];
        ini.N_days = infer_days; % reset N_days
    end
    
    params_temp2 = mvnrnd(params_temp1, diag(proposal_dist), 1);
    if min(params_temp2((1+2*n_country) : end)) >= 0       
        A = fmin_europe_data(params_temp1, newly_confirmed_data_ave, ini, 'pois', dt, mu, beta, trainDays_begin3, trainDays_end3, sigma) - ...
            fmin_europe_data(params_temp2, newly_confirmed_data_ave, ini, 'pois', dt, mu, beta, trainDays_begin3, trainDays_end3, sigma);
        Avec(i) = A;    
        A = exp(A);
    else
        A = 0;
    end
    if A >= 1
        params_temp1 = params_temp2;
    else
        U = unifrnd(0, 1);
        if U <= A
            params_temp1 = params_temp2;
        end
    end
    
     if mod(i, 5000) == 0
         fprintf('iteration: %d, A = %f, lambda_{Italy,1}: %f, lambda_{Italy,2}: %f, \n', i, Avec(i), params_temp1(2*n_country+7), ...
             params_temp1(3*n_country+7));
     end
end
toc;

%% plotting section for specific country



country_ind_temp = 4;
% country_ind_temp = 5;
% country_ind_temp = 7;
country_names(country_ind_temp)

figure(10), clf;

FontSize = 30;

muind1 = ind2_plot;
muind2 = ind2_plot-3;
muind3 = ind2_plot+3;

y_reg_wHwoMwGL1 = simu_europe_data_generate_determ(paramsmat2(muind1,:), ini, dt, trainDays_begin3, testDays_end);
y_reg_wHwoMwGL2 = simu_europe_data_generate_determ(paramsmat2(muind2,:), ini, dt, trainDays_begin3, testDays_end);
y_reg_wHwoMwGL3 = simu_europe_data_generate_determ(paramsmat2(muind3,:), ini, dt, trainDays_begin3, testDays_end);
y_reg_wHwoM     = simu_europe_data_generate_determ(paramsmat2(1,:), ini, dt, trainDays_begin3, testDays_end);

subplot(2,2,1)
grid on
hold on
plot(1: totDays, y_reg_wHwoMwGL1(:, country_ind_temp), 'x-', 'LineWidth', 2, 'Color', '#000099', 'MarkerSize', 5);
plot(1: totDays, y_reg_wHwoMwGL2(:, country_ind_temp), 'x-', 'LineWidth', 2, 'Color', '#0066ff', 'MarkerSize', 5);
plot(1: totDays, y_reg_wHwoMwGL3(:, country_ind_temp), 'x-', 'LineWidth', 2, 'Color', '#99ccff', 'MarkerSize', 5);
plot(1: totDays, newly_confirmed_data_ave(:, country_ind_temp), 'o-', 'LineWidth', 2, 'Color', '#ff9900', 'MarkerSize', 5);
xline(trainDays_end3, 'Linewidth', 2);
legend_string = strings(4,1);
legend_string(1) = strcat("Model 5, $\mu = 10^{", num2str(log10(muvec(muind1)),'%.1f'), "}$");
legend_string(2) = strcat("Model 5, $\mu = 10^{", num2str(log10(muvec(muind2)),'%.1f'), "}$");
legend_string(3) = strcat("Model 5, $\mu = 10^{", num2str(log10(muvec(muind3)),'%.1f'), "}$");
legend_string(4) = "True";
legend(legend_string,'Location','best', 'Interpreter', 'latex', 'FontSize', 24);
set(gca, 'FontSize', FontSize);
xlabel('day', 'Interpreter', 'latex')
ylabel('Newly Confirmed Cases', 'Interpreter', 'latex')
title(sprintf('True and predicted trajectories with Model 3'''), 'Interpreter', 'latex', 'FontSize', 33)


subplot(2,2,2);
grid on
hold on
plot(1: totDays, y_reg_woHwoM(:, country_ind_temp), '+-', 'LineWidth', 2, 'Color', '#cccc00', 'MarkerSize', 5);
plot(1: totDays, y_reg_wHwoM(:, country_ind_temp), 'd-', 'LineWidth', 2, 'Color', '#009933', 'MarkerSize', 5);
plot(1: totDays, y_reg_wHwoMwGL1(:, country_ind_temp), 'x-', 'LineWidth', 2, 'Color', '#000099', 'MarkerSize', 5);
plot(1: totDays, newly_confirmed_data_ave(:, country_ind_temp), 'o-', 'LineWidth', 2, 'Color', '#ff9900', 'MarkerSize', 5);
xline(trainDays_end3, 'Linewidth', 2);
% legend('woHwoM', 'woHwM', 'wHwoM', 'wHwM', 'wHwM,wGL','True','Location','best');
legend('Model 1''', 'Model 2''', ...
    strcat("Model 3'', $\mu = 10^{", num2str(log10(muvec(muind1)),'%.1f'), "}$"),...
    'True','Location','best', 'Interpreter', 'latex', 'FontSize', 24);
set(gca, 'FontSize', FontSize);
xlabel('day', 'Interpreter', 'latex')
ylabel('Newly Confirmed Cases', 'Interpreter', 'latex')
title(sprintf('True and predicted trajectories for all models'), 'Interpreter', 'latex', 'FontSize', 33)

subplot(2,2,3);
grid on
hold on
for j = 1: numtrasamp
    s = scatter(1: tot_days-7, reshape(traject_post(j,:,country_ind_temp),totDays,1), 'filled', 'MarkerFaceColor', '#EDB120', 'SizeData', 30);
    alpha(s, .1);
end
for j = 1: numsamp
    t = scatter(1: tot_days-7, reshape(trajec_mat_one_Theta(j,:,country_ind_temp),totDays,1), 'filled', 'c', 'SizeData', 30);
    alpha(t, .1);
end
plot(1: totDays, y_reg_wHwoMwGL1(:, country_ind_temp), 'x-', 'LineWidth', 2, 'Color', '#000099', 'MarkerSize', 5);
plot(1: totDays, newly_confirmed_data_ave(:, country_ind_temp), 'o-', 'LineWidth', 2, 'Color', '#ff9900', 'MarkerSize', 5);
xline(trainDays_end3, 'Linewidth', 2);
set(gca, 'FontSize', FontSize);
xlabel('day', 'Interpreter', 'latex')
ylabel('Newly Confirmed Cases', 'Interpreter', 'latex')
% title(sprintf('Sampled trajectories'), 'Interpreter', 'latex', 'FontSize', 30)
title(strcat("Sampled trajectories with Model 3', $\mu = 10^{", num2str(log10(muvec(muind1)),'%.1f'), "}$"), 'Interpreter', 'latex', 'FontSize', 33)


subplot(2,2,4);
grid on
hold on
plot(1: totDays, abs(y_reg_woHwoM(:, country_ind_temp)-newly_confirmed_data_ave(:, country_ind_temp)), '+-', 'LineWidth', 2, 'Color', '#cccc00', 'MarkerSize', 5);
plot(1: totDays, abs(y_reg_wHwoM(:, country_ind_temp)-newly_confirmed_data_ave(:, country_ind_temp)), 'd-', 'LineWidth', 2, 'Color', '#009933', 'MarkerSize', 5);
plot(1: totDays, abs(y_reg_wHwoMwGL1(:, country_ind_temp)-newly_confirmed_data_ave(:, country_ind_temp)), 'x-', 'LineWidth', 2, 'Color', '#000099', 'MarkerSize', 5);
plot(1: totDays, abs(y_reg_wHwoMwGL2(:, country_ind_temp)-newly_confirmed_data_ave(:, country_ind_temp)), 'x-', 'LineWidth', 2, 'Color', '#0066ff', 'MarkerSize', 5);
plot(1: totDays, abs(y_reg_wHwoMwGL3(:, country_ind_temp)-newly_confirmed_data_ave(:, country_ind_temp)), 'x-', 'LineWidth', 2, 'Color', '#99ccff', 'MarkerSize', 5);
xline(trainDays_end3, 'Linewidth', 2);
% legend('woHwoM', 'woHwM', 'wHwoM', 'wHwM', 'wHwM,wGL','Location','best');
legend('Model 1''', 'Model 2''', ...
    strcat("Model 3'', $\mu = 10^{", num2str(log10(muvec(muind1)),'%.1f'), "}$"),...
    strcat("Model 3'', $\mu = 10^{", num2str(log10(muvec(muind2)),'%.1f'), "}$"),...
    strcat("Model 3'', $\mu = 10^{", num2str(log10(muvec(muind3)),'%.1f'), "}$"),...
    'Location','best', 'Interpreter', 'latex', 'FontSize', 24);
set(gca, 'FontSize', FontSize);
xlabel('day', 'Interpreter', 'latex')
ylabel('Absolute Errors', 'Interpreter', 'latex')
title(sprintf('Error of predicted trajectories for all models'), 'Interpreter', 'latex', 'FontSize', 33)

sgtitle(sprintf("True and predicted newly confirmed cases in %s\n", country_names(country_ind_temp)), 'Interpreter', 'latex', "FontSize", 42);



%%

figure(35), clf;

id_country = 7;
country_names(id_country)
muind  = ind2_plot;

subplot(1,2,1);
lambda_post = params_post(1:iterations, 2*n_country+id_country);
mean(lambda_post)
std(lambda_post)
paramsmat2(muind,2*n_country+id_country)
h2 = histogram(lambda_post, 'Normalization', 'pdf');
xline(paramsmat2(muind,2*n_country+id_country), 'black', 'Linewidth', 4);
set(gca, 'FontSize', 30);
title('$\lambda_{\rm Italy}^{(1)}$', 'FontSize', 36, 'Interpreter', 'latex');
yticks([0 5e2 1e3 1.5e3 2e3 2.5e3 3e3 3.5e3])
yticklabels({'0', '500', '1,000', '1,500', '2,000', '2,500', '3,000', '3,500'})
grid on


subplot(1,2,2);
lambda_post = params_post(1:iterations, 3*n_country+id_country);
mean(lambda_post)
std(lambda_post)
paramsmat2(muind,3*n_country+id_country)
h2 = histogram(lambda_post, 'Normalization', 'pdf');
xline(paramsmat2(muind,3*n_country+id_country), 'black', 'Linewidth', 4);
set(gca, 'FontSize', 30);
title('$\lambda_{\rm Italy}^{(1)}$', 'FontSize', 36, 'Interpreter', 'latex');
grid on

sgtitle("Estimated posterior distribution for $\lambda$ in Italy", "FontSize", 42, 'Interpreter', 'latex');



%%

































