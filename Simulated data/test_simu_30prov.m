% thirty provinces

% validation set last for the same time period as the testing set, for example:
% training   set: 1-10days
% validation set: 11-30days
% testing    set: 11-30 days
% trajectories from which these sets are from are all independent

% predicted testing trajectories are computed in a "hidden way":
% using the fitted [I0, E0, lambda], compute the trajectory from the 1st
% day instead of using the ground truth information from the 10th day
% namely, the only known facts are newly confirmed cases


% fix some bugs 


clear all;

% rng(621);
rng(100000);

%%

provNum = 30;
provPop = 1e6 * ones(1, provNum);



% randomly generate the group numbers
group_num_min = 3; group_num_max = 3;
group_num = randi([group_num_min group_num_max], 1);
fprintf('group number = %d\n\n', group_num);
% randomly divide the provinces into groups
groupIndex = randi(group_num,provNum,1); % Generate an index that assigns each row to one of the groups
group_cell = cell(group_num,1);
group_ind_temp = 1:provNum;
for i = 1:group_num
    group_cell{i} = group_ind_temp(groupIndex == i);
end
group_cell

%------------------------------------

% design travelling matrix
traVal_within  = 5e3; % travelling volume between provinces in the same group
traVal_between  = 5e3; % travelling volume between provinces not in the same group
traMat = traVal_between * ones(provNum, provNum);
for i = 1:group_num
    ind_temp = group_cell{i};
    traMat(ind_temp, ind_temp) = traMat(ind_temp, ind_temp) + (traVal_within-traVal_between)*ones(length(ind_temp));
end
traMat  = traMat - diag(diag(traMat)); % traval from prov i to prov i = 0

%------------------------------------

delta     = 0.14;
gamma     = 0.14;


%------------------------------------
% first randomly generate "the transmission parameter for each group". It
% is the mean/median of the transmission parameters for each province in that group (which are also random)

% group_lambda_min = 0.2; group_lambda_max = 0.5;
% group_lambda_vec = group_lambda_min + (group_lambda_max-group_lambda_min) * rand(group_num,1);
group_lambda_vec = [0.3 0.4 0.5];
group_lambda_vec 
%------------------------------------
% randomly generate the transmission parameters for provinces in each
% group, which is "centered" around "the transmission parameter sampled for
% this group"
lambdavec = zeros(1,provNum);
lambda_sigma_norm_random = 0.01;
for i = 1:group_num
    ind_temp = group_cell{i};
    lambdavec(ind_temp) = normrnd(group_lambda_vec(i), lambda_sigma_norm_random, length(ind_temp), 1);
end

%------------------------------------

% print the generated transmission parameters
for i = 1:group_num
    fprintf('group %d, the "mean of transmission parameters" = %.3f\n\n', i, group_lambda_vec(i));
    ind_temp = group_cell{i};
    for j = 1:length(ind_temp)
        fprintf('prov %d, transmission parameter %.3f\n', ind_temp(j), lambdavec(ind_temp(j)));
    end
    fprintf('\n\n');
end

%------------------------------------

E0      = 30 * ones(1, provNum);
I0      = 10 * ones(1, provNum);
param_true   = [I0, E0, lambdavec];

totDays  = 20;
startDay = 1;
trainDays_end   = 10;
validDays_start = trainDays_end+1;
validDays_end   = totDays;
testDays_start  = trainDays_end+1;
testDays_end    = totDays;

ini.totPop = provPop;
%ini.traVal = traVal;

%------------------------------------
% % mismatch of group
% groupIndex2 = groupIndex;
% group_ind_temp = 1:provNum;
% groupIndex2(min(group_ind_temp(groupIndex==1))) = 2;
% groupIndex2(min(group_ind_temp(groupIndex==2))) = 3;
% groupIndex2(min(group_ind_temp(groupIndex==3))) = 1;
% group_cell2 = cell(group_num,1);
% for i = 1:group_num
%     group_cell2{i} = group_ind_temp(groupIndex2 == i);
% end
% 
% ini.part   = group_cell2;
% 
% group_cell2


ini.part   = group_cell;

%------------------------------------
ini.delta  = delta;
ini.gamma  = gamma;
%ini.lambda = lambda;

dt      = 0.1;
%day_end = 25;

%%
numrep = 1e2;
sigma  = 10^(-6);

traj_obs_mat = zeros(totDays, 5*provNum, numrep); % store observable trajectories
traj_obs_test_mat  = zeros(totDays, 5*provNum, numrep); % store observable trajectories
traj_obs_valid_mat = zeros(totDays, 5*provNum, numrep); % store observable trajectories

traj_train_woHwoM_mat = zeros(trainDays_end-startDay+1, provNum, numrep); % store fitted trajectories, wo Hetereo wo Migrat
traj_train_woHwM_mat  = zeros(trainDays_end-startDay+1, provNum, numrep); % store fitted trajectories, wo Hetereo w/ Migrat
traj_train_wHwoM_mat  = zeros(trainDays_end-startDay+1, provNum, numrep); % store fitted trajectories, w/ Hetereo wo Migrat
traj_train_wHwM_mat   = zeros(trainDays_end-startDay+1, provNum, numrep); % store fitted trajectories, w/ Hetereo w/ Migrat

traj_valid_woHwoM_mat = zeros(validDays_end-validDays_start+1, provNum, numrep); % store fitted trajectories, wo Hetereo wo Migrat
traj_valid_woHwM_mat  = zeros(validDays_end-validDays_start+1, provNum, numrep); % store fitted trajectories, wo Hetereo w/ Migrat
traj_valid_wHwoM_mat  = zeros(validDays_end-validDays_start+1, provNum, numrep); % store fitted trajectories, w/ Hetereo wo Migrat
traj_valid_wHwM_mat   = zeros(validDays_end-validDays_start+1, provNum, numrep); % store fitted trajectories, w/ Hetereo w/ Migrat

traj_test_woHwoM_mat = zeros(testDays_end-testDays_start+1, provNum, numrep); % store fitted trajectories, wo Hetereo wo Migrat
traj_test_woHwM_mat  = zeros(testDays_end-testDays_start+1, provNum, numrep); % store fitted trajectories, wo Hetereo w/ Migrat
traj_test_wHwoM_mat  = zeros(testDays_end-testDays_start+1, provNum, numrep); % store fitted trajectories, w/ Hetereo wo Migrat
traj_test_wHwM_mat   = zeros(testDays_end-testDays_start+1, provNum, numrep); % store fitted trajectories, w/ Hetereo w/ Migrat

param_woHwoM_mat = zeros(2*provNum+1, numrep); % store inferred params, wo Hetereo wo Migrat
param_woHwM_mat  = zeros(2*provNum+1, numrep); % store inferred params, wo Hetereo w/ Migrat
param_wHwoM_mat  = zeros(3*provNum, numrep);   % store inferred params, w/ Hetereo wo Migrat
param_wHwM_mat   = zeros(3*provNum, numrep);   % store inferred params, w/ Hetereo w/ Migrat, no GL

%--------------------------------------------------------------------------
% muvec = 10.^ (-3:0.1:-1)*5000;%10.^ (-2:0.2:2);
% muvec = 10.^ (1.5:0.1:3.5)
% muvec = 10.^ (0:0.1:3.5)
muvec = 10.^ (0:0.1:3.5)
% muvec = 10.^(-2);% * 5000;
beta  = 0.2; % mu0 = beta * mu1
%--------------------------------------------------------------------------

%param_wHwMwGL_mat = zeros(length(mu0vec), length(mu1vec), 3*provNum, numrep);
param_wHwMwGL_mat = zeros(length(muvec), 3*provNum, numrep);
% store inferred params, w/ Hetereo w/ Migrat, w/ GL

trainErr_woHwoM_mat1  = zeros(numrep,1); % relative error with weighted average
trainErr_woHwM_mat1   = zeros(numrep,1);
trainErr_wHwoM_mat1   = zeros(numrep,1);
trainErr_wHwM_mat1    = zeros(numrep,1);
trainErr_wHwMwGL_mat1 = zeros(length(muvec), numrep);

validErr_woHwoM_mat1  = zeros(numrep,1);
validErr_woHwM_mat1   = zeros(numrep,1);
validErr_wHwoM_mat1   = zeros(numrep,1);
validErr_wHwM_mat1    = zeros(numrep,1);
validErr_wHwMwGL_mat1 = zeros(length(muvec), numrep);

testErr_woHwoM_mat1  = zeros(numrep,1);
testErr_woHwM_mat1   = zeros(numrep,1);
testErr_wHwoM_mat1   = zeros(numrep,1);
testErr_wHwM_mat1    = zeros(numrep,1);
testErr_wHwMwGL_mat1 = zeros(length(muvec), numrep);

trainErr_woHwoM_mat2  = zeros(numrep,1);
trainErr_woHwM_mat2   = zeros(numrep,1);
trainErr_wHwoM_mat2   = zeros(numrep,1);
trainErr_wHwM_mat2    = zeros(numrep,1);
trainErr_wHwMwGL_mat2 = zeros(length(muvec), numrep);

validErr_woHwoM_mat2  = zeros(numrep,1); % relative error with weighted average
validErr_woHwM_mat2   = zeros(numrep,1);
validErr_wHwoM_mat2   = zeros(numrep,1);
validErr_wHwM_mat2    = zeros(numrep,1);
validErr_wHwMwGL_mat2 = zeros(length(muvec), numrep);

testErr_woHwoM_mat2  = zeros(numrep,1);
testErr_woHwM_mat2   = zeros(numrep,1);
testErr_wHwoM_mat2   = zeros(numrep,1);
testErr_wHwM_mat2    = zeros(numrep,1);
testErr_wHwMwGL_mat2 = zeros(length(muvec), numrep);

trainErr_woHwoM_mat3  = zeros(numrep,1); % mean square relative error with weighted average
trainErr_woHwM_mat3   = zeros(numrep,1);
trainErr_wHwoM_mat3   = zeros(numrep,1);
trainErr_wHwM_mat3    = zeros(numrep,1);
trainErr_wHwMwGL_mat3 = zeros(length(muvec), numrep);

validErr_woHwoM_mat3  = zeros(numrep,1);
validErr_woHwM_mat3   = zeros(numrep,1);
validErr_wHwoM_mat3   = zeros(numrep,1);
validErr_wHwM_mat3    = zeros(numrep,1);
validErr_wHwMwGL_mat3 = zeros(length(muvec), numrep);

testErr_woHwoM_mat3  = zeros(numrep,1);
testErr_woHwM_mat3   = zeros(numrep,1);
testErr_wHwoM_mat3   = zeros(numrep,1);
testErr_wHwM_mat3    = zeros(numrep,1);
testErr_wHwMwGL_mat3 = zeros(length(muvec), numrep);

trainErr_woHwoM_mat4  = zeros(numrep,1); % mean square relative error with simple average
trainErr_woHwM_mat4   = zeros(numrep,1);
trainErr_wHwoM_mat4   = zeros(numrep,1);
trainErr_wHwM_mat4    = zeros(numrep,1);
trainErr_wHwMwGL_mat4 = zeros(length(muvec), numrep);

validErr_woHwoM_mat4  = zeros(numrep,1);
validErr_woHwM_mat4   = zeros(numrep,1);
validErr_wHwoM_mat4   = zeros(numrep,1);
validErr_wHwM_mat4    = zeros(numrep,1);
validErr_wHwMwGL_mat4 = zeros(length(muvec), numrep);

testErr_woHwoM_mat4  = zeros(numrep,1);
testErr_woHwM_mat4   = zeros(numrep,1);
testErr_wHwoM_mat4   = zeros(numrep,1);
testErr_wHwM_mat4    = zeros(numrep,1);
testErr_wHwMwGL_mat4 = zeros(length(muvec), numrep);


%%

% rng(2022)

E_pre = 50;
I_pre = 10;
lambda_pre = 0.2;

tic;
for rep = 1: numrep
    
    fprintf('the %d th replica\n', rep);
    
    fprintf('generating data\n');
    ini.traMat = traMat;
    %y_true  = max(SEIR_data_gen_test_pred_inf_random_2(param, ini, totDays, startDay),0.01);
    
    y_true_test  = max(SEIR_data_gen_test_pred_inf_random_dt_2(param_true, ini, totDays, startDay, dt),1);
    y_true_valid = max(SEIR_data_gen_test_pred_inf_random_dt_2(param_true, ini, totDays, startDay, dt),1);%y_true;%max(SEIR_data_gen_test_pred_inf_random_dt_2(param, ini, totDays, startDay, dt),1);
    y_true       = max(SEIR_data_gen_test_pred_inf_random_dt_2(param_true, ini, totDays, startDay, dt),1);
    
    y_obs   = y_true(:, 1: (5*provNum)); % deltaC, deltaCR,  E, I, R
    traj_obs_mat(:, :, rep) = y_obs;
    traj_obs_test_mat(:, :, rep) = y_true_test(:, 1: (5*provNum));
    traj_obs_valid_mat(:, :, rep) = y_true_valid(:, 1: (5*provNum));
    
    %%
    %-----------------------------------------------------------------------------------------------------
    % wo Hetero, wo Migrat
    fprintf('wo Hetero, wo Migrat\n');
    ini.traMat = zeros(provNum);
    
    f1 = @(x)SEIR_fminunc_test_pred_inf_wo_h_2(x, y_obs, ini, 'pois', dt, trainDays_end, startDay );
    opts1 = optimoptions('fmincon');
    opts1.MaxIterations =  2000;
    opts1.MaxFunctionEvaluations = 10^6;
    opts1.Display = 'off';%'iter-detailed';
    opts1.OptimalityTolerance = 1e-10;
    %     lb = [ones(1,2*provNum), 0];
    %     ub = ones(1,2*provNum+1) * 1000;
    lb = [ones(1,2*provNum), 0];
    ub = [ones(1,2*provNum) * 100, 1];
    params_pre = [I_pre  * ones(1, length(provPop)), E_pre * ones(1, length(provPop)), lambda_pre];
    [params1, ~] = fmincon(f1, params_pre, [],[],[],[],lb,ub,[],opts1);
    param_woHwoM_mat(:, rep) = params1';
    
    
    R_train = zeros(1, provNum);
    params1_test = params1;
    R_test  = zeros(1, provNum);
    params1_valid = params1;
    R_valid = zeros(1, provNum);
    y_reg_train = SEIR_data_gen_test_pred_inf_determ_wo_h_2(params1, ini, dt, trainDays_end-startDay+1, 1, R_train);
    y_reg_valid = SEIR_data_gen_test_pred_inf_determ_wo_h_2(params1_valid, ini, dt, validDays_end-startDay+1, 1, R_valid);
    y_reg_valid = y_reg_valid(validDays_start:validDays_end,:);
    y_reg_test  = SEIR_data_gen_test_pred_inf_determ_wo_h_2(params1_test, ini, dt, testDays_end-startDay+1, 1, R_test);
    y_reg_test  = y_reg_test(testDays_start:testDays_end,:);
    traj_train_woHwoM_mat(:, :, rep) = y_reg_train(:, 1: provNum);
    traj_valid_woHwoM_mat(:, :, rep) = y_reg_valid(:, 1: provNum);
    traj_test_woHwoM_mat(:, :, rep)  = y_reg_test(:, 1: provNum);
    %-------------------------
    trainErr_1 = abs(y_reg_train(:, 1: provNum) - y_true(1: trainDays_end-startDay+1, 1: provNum));
    trainErr_vec_1 = (sum(trainErr_1))' ./ (y_true(trainDays_end-startDay+1, (3*provNum+1):(4*provNum))  ...
        -   y_true(1, (3*provNum+1):(4*provNum))   )' ;
    validErr_1  = abs(y_reg_valid(:, 1: provNum) ...
        - y_true_valid(validDays_start: validDays_end, 1: provNum));
    validErr_vec_1 = (sum(validErr_1))' ./ (y_true_valid(validDays_end, (3*provNum+1):(4*provNum))  ...
        -   y_true_valid(validDays_start, (3*provNum+1):(4*provNum))   )' ;
    testErr_1  = abs(y_reg_test(:, 1: provNum) ...
        - y_true_test(testDays_start: testDays_end, 1: provNum));
    testErr_vec_1 = (sum(testErr_1))' ./ (y_true_test(testDays_end, (3*provNum+1):(4*provNum))  ...
        -   y_true_test(testDays_start, (3*provNum+1):(4*provNum))   )' ;
    trainErr_woHwoM_mat1(rep) = sum(trainErr_vec_1);
    validErr_woHwoM_mat1(rep)  = sum(validErr_vec_1);
    testErr_woHwoM_mat1(rep)  = sum(testErr_vec_1);
    %--------------------------
    trainErr_2  = abs(y_reg_train(:, 1: provNum) ...
        - y_true(1: trainDays_end-startDay+1, 1: provNum)) ...
        ./ y_true(1: trainDays_end-startDay+1, 1: provNum);
    trainErr_vec_2 = (sum(trainErr_2))' / trainDays_end ;
    validErr_2  = abs(y_reg_valid(:, 1: provNum) ...
        - y_true_valid(validDays_start: validDays_end, 1: provNum)) ...
        ./ y_true_valid(validDays_start: validDays_end, 1: provNum);
    validErr_vec_2 = (sum(validErr_2))' / (validDays_end - validDays_start + 1) ;
    testErr_2  = abs(y_reg_test(:, 1: provNum) ...
        - y_true_test(testDays_start: testDays_end, 1: provNum)) ...
        ./ y_true_test(testDays_start: testDays_end, 1: provNum);
    testErr_vec_2 = (sum(testErr_2))' / (testDays_end - testDays_start + 1) ;
    trainErr_woHwoM_mat2(rep) = sum(trainErr_vec_2);
    validErr_woHwoM_mat2(rep) = sum(validErr_vec_2);
    testErr_woHwoM_mat2(rep)  = sum(testErr_vec_2);
    %--------------------------
    trainErr_3  = ( abs(y_reg_train(:, 1: provNum) ...
        - y_true(1: trainDays_end-startDay+1, 1: provNum)) ...
        ./ y_true(1: trainDays_end-startDay+1, 1: provNum) ).^2;
    trainErr_vec_3 = sqrt((sum(trainErr_3 .* y_true(1: trainDays_end-startDay+1, 1: provNum)...
        ./ (ones(trainDays_end,1)*(y_true(trainDays_end-startDay+1, (3*provNum+1):(4*provNum))  ...
        -   y_true(1, (3*provNum+1):(4*provNum)))) ))' ) ;
    validErr_3  = ( abs(y_reg_valid(:, 1: provNum) ...
        - y_true_valid(validDays_start: validDays_end, 1: provNum)) ...
        ./ y_true_valid(validDays_start: validDays_end, 1: provNum) ).^2;
    validErr_vec_3 = sqrt((sum(validErr_3 .* y_true_valid(validDays_start: validDays_end, 1: provNum)...
        ./ (ones(validDays_end-validDays_start+1,1)*(y_true_valid(validDays_end, (3*provNum+1):(4*provNum))  ...
        -   y_true_valid(validDays_start, (3*provNum+1):(4*provNum)))) ))' ) ;
    testErr_3  = ( abs(y_reg_test(:, 1: provNum) ...
        - y_true_test(testDays_start: testDays_end, 1: provNum)) ...
        ./ y_true_test(testDays_start: testDays_end, 1: provNum) ).^2;
    testErr_vec_3 = sqrt((sum(testErr_3 .* y_true_test(testDays_start: testDays_end, 1: provNum)...
        ./ (ones(testDays_end-testDays_start+1,1)*(y_true_test(testDays_end, (3*provNum+1):(4*provNum))  ...
        -   y_true_test(testDays_start, (3*provNum+1):(4*provNum)))) ))' ) ;
    trainErr_woHwoM_mat3(rep) = sum(trainErr_vec_3);
    validErr_woHwoM_mat3(rep) = sum(validErr_vec_3);
    testErr_woHwoM_mat3(rep)  = sum(testErr_vec_3);
    %--------------------------
    trainErr_vec_4 = sqrt( (sum(trainErr_3))' / (trainDays_end-startDay+1) );
    validErr_vec_4 = sqrt( (sum(validErr_3))' / (validDays_end - validDays_start+1) );
    testErr_vec_4 = sqrt( (sum(testErr_3))' / (testDays_end-testDays_start+1) );
    trainErr_woHwoM_mat4(rep) = sum(trainErr_vec_4);
    validErr_woHwoM_mat4(rep) = sum(validErr_vec_4);
    testErr_woHwoM_mat4(rep) = sum(testErr_vec_4);
    
    %%
    %-----------------------------------------------------------------------------------------------------
    % wo Hetero, w/ Migrat
    fprintf('wo Hetero, w/ Migrat\n');
    ini.traMat = traMat;
    
    f1 = @(x)SEIR_fminunc_test_pred_inf_wo_h_2(x, y_obs, ini, 'pois', dt, trainDays_end, startDay );
    opts1 = optimoptions('fmincon');
    opts1.MaxIterations =  2000;
    opts1.MaxFunctionEvaluations = 10^6;
    opts1.Display = 'off';%'iter-detailed';
    opts1.OptimalityTolerance = 1e-10;
    %     lb = [ones(1,2*provNum), 0];
    %     ub = ones(1,2*provNum+1) * 1000;
    lb = [ones(1,2*provNum), 0];
    ub = [ones(1,2*provNum) * 100, 1];
    params_pre = [I_pre  * ones(1, length(provPop)), E_pre * ones(1, length(provPop)), lambda_pre];
    [params1, ~] = fmincon(f1, params_pre, [],[],[],[],lb,ub,[],opts1);
    param_woHwM_mat(:, rep) = params1';
    
    R_train = zeros(1, provNum);
    params1_test = params1;
    R_test  = zeros(1, provNum);
    params1_valid = params1;
    R_valid = zeros(1, provNum);
    y_reg_train = SEIR_data_gen_test_pred_inf_determ_wo_h_2(params1, ini, dt, trainDays_end-startDay+1, 1, R_train);
    y_reg_valid = SEIR_data_gen_test_pred_inf_determ_wo_h_2(params1_valid, ini, dt, validDays_end-startDay+1, 1, R_valid);
    y_reg_valid = y_reg_valid(validDays_start:validDays_end,:);
    y_reg_test  = SEIR_data_gen_test_pred_inf_determ_wo_h_2(params1_test, ini, dt, testDays_end-startDay+1, 1, R_test);
    y_reg_test  = y_reg_test(testDays_start:testDays_end,:);
    traj_train_woHwM_mat(:, :, rep) = y_reg_train(:, 1: provNum);
    traj_valid_woHwM_mat(:, :, rep) = y_reg_valid(:, 1: provNum);
    traj_test_woHwM_mat(:, :, rep)  = y_reg_test(:, 1: provNum);
    %-------------------------
    trainErr_1 = abs(y_reg_train(:, 1: provNum) - y_true(1: trainDays_end-startDay+1, 1: provNum));
    trainErr_vec_1 = (sum(trainErr_1))' ./ (y_true(trainDays_end-startDay+1, (3*provNum+1):(4*provNum))  ...
        -   y_true(1, (3*provNum+1):(4*provNum))   )' ;
    validErr_1  = abs(y_reg_valid(:, 1: provNum) ...
        - y_true_valid(validDays_start: validDays_end, 1: provNum));
    validErr_vec_1 = (sum(validErr_1))' ./ (y_true_valid(validDays_end, (3*provNum+1):(4*provNum))  ...
        -   y_true_valid(validDays_start, (3*provNum+1):(4*provNum))   )' ;
    testErr_1  = abs(y_reg_test(:, 1: provNum) ...
        - y_true_test(testDays_start: testDays_end, 1: provNum));
    testErr_vec_1 = (sum(testErr_1))' ./ (y_true_test(testDays_end, (3*provNum+1):(4*provNum))  ...
        -   y_true_test(testDays_start, (3*provNum+1):(4*provNum))   )' ;
    trainErr_woHwM_mat1(rep) = sum(trainErr_vec_1);
    validErr_woHwM_mat1(rep)  = sum(validErr_vec_1);
    testErr_woHwM_mat1(rep)  = sum(testErr_vec_1);
    %--------------------------
    trainErr_2  = abs(y_reg_train(:, 1: provNum) ...
        - y_true(1: trainDays_end-startDay+1, 1: provNum)) ...
        ./ y_true(1: trainDays_end-startDay+1, 1: provNum);
    trainErr_vec_2 = (sum(trainErr_2))' / trainDays_end ;
    validErr_2  = abs(y_reg_valid(:, 1: provNum) ...
        - y_true_valid(validDays_start: validDays_end, 1: provNum)) ...
        ./ y_true_valid(validDays_start: validDays_end, 1: provNum);
    validErr_vec_2 = (sum(validErr_2))' / (validDays_end - validDays_start + 1) ;
    testErr_2  = abs(y_reg_test(:, 1: provNum) ...
        - y_true_test(testDays_start: testDays_end, 1: provNum)) ...
        ./ y_true_test(testDays_start: testDays_end, 1: provNum);
    testErr_vec_2 = (sum(testErr_2))' / (testDays_end - testDays_start + 1) ;
    trainErr_woHwM_mat2(rep) = sum(trainErr_vec_2);
    validErr_woHwM_mat2(rep) = sum(validErr_vec_2);
    testErr_woHwM_mat2(rep)  = sum(testErr_vec_2);
    %--------------------------
    trainErr_3  = ( abs(y_reg_train(:, 1: provNum) ...
        - y_true(1: trainDays_end-startDay+1, 1: provNum)) ...
        ./ y_true(1: trainDays_end-startDay+1, 1: provNum) ).^2;
    trainErr_vec_3 = sqrt((sum(trainErr_3 .* y_true(1: trainDays_end-startDay+1, 1: provNum)...
        ./ (ones(trainDays_end,1)*(y_true(trainDays_end-startDay+1, (3*provNum+1):(4*provNum))  ...
        -   y_true(1, (3*provNum+1):(4*provNum)))) ))' ) ;
    validErr_3  = ( abs(y_reg_valid(:, 1: provNum) ...
        - y_true_valid(validDays_start: validDays_end, 1: provNum)) ...
        ./ y_true_valid(validDays_start: validDays_end, 1: provNum) ).^2;
    validErr_vec_3 = sqrt((sum(validErr_3 .* y_true_valid(validDays_start: validDays_end, 1: provNum)...
        ./ (ones(validDays_end-validDays_start+1,1)*(y_true_valid(validDays_end, (3*provNum+1):(4*provNum))  ...
        -   y_true_valid(validDays_start, (3*provNum+1):(4*provNum)))) ))' ) ;
    testErr_3  = ( abs(y_reg_test(:, 1: provNum) ...
        - y_true_test(testDays_start: testDays_end, 1: provNum)) ...
        ./ y_true_test(testDays_start: testDays_end, 1: provNum) ).^2;
    testErr_vec_3 = sqrt((sum(testErr_3 .* y_true_test(testDays_start: testDays_end, 1: provNum)...
        ./ (ones(testDays_end-testDays_start+1,1)*(y_true_test(testDays_end, (3*provNum+1):(4*provNum))  ...
        -   y_true_test(testDays_start, (3*provNum+1):(4*provNum)))) ))' ) ;
    trainErr_woHwM_mat3(rep) = sum(trainErr_vec_3);
    validErr_woHwM_mat3(rep) = sum(validErr_vec_3);
    testErr_woHwM_mat3(rep)  = sum(testErr_vec_3);
    %--------------------------
    trainErr_vec_4 = sqrt( (sum(trainErr_3))' / (trainDays_end-startDay+1) );
    validErr_vec_4 = sqrt( (sum(validErr_3))' / (validDays_end - validDays_start+1) );
    testErr_vec_4 = sqrt( (sum(testErr_3))' / (testDays_end-testDays_start+1) );
    trainErr_woHwM_mat4(rep) = sum(trainErr_vec_4);
    validErr_woHwM_mat4(rep) = sum(validErr_vec_4);
    testErr_woHwM_mat4(rep) = sum(testErr_vec_4);
    
    %%
    %-----------------------------------------------------------------------------------------------------
    % w/ Hetero, wo Migrat
    fprintf('w/ Hetero, wo Migrat\n');
    ini.traMat = zeros(provNum); mu0 = 0; mu1 = 0;
    
    f1 = @(x)SEIR_fminunc_test_pred_inf_multi_2(x, y_obs, ini, 'pois', mu0, mu1, dt, trainDays_end, startDay, 0 );
    opts1 = optimoptions('fmincon');
    opts1.MaxIterations =  2000;
    opts1.MaxFunctionEvaluations = 10^6;
    opts1.Display = 'off';%'iter-detailed';
    opts1.OptimalityTolerance = 1e-10;
    lb = [ones(1,2*provNum), zeros(1,provNum)];
    ub = [ones(1,2*provNum) * 100, ones(1,provNum) * 1];
    params_pre = [I_pre  * ones(1, length(provPop)), E_pre * ones(1, length(provPop)), lambda_pre*ones(1,length(provPop))];
    [params1, ~] = fmincon(f1, params_pre, [],[],[],[],lb,ub,[],opts1);
    param_wHwoM_mat(:, rep) = params1';
    %     params1(9:12)
    
    R_train = zeros(1, provNum);
    params1_test = params1;
    R_test  = zeros(1, provNum);
    params1_valid = params1;
    R_valid = zeros(1, provNum);
    y_reg_train = SEIR_data_gen_test_pred_inf_determ_2(params1, ini, dt, trainDays_end-startDay+1, 1, R_train);
    y_reg_valid = SEIR_data_gen_test_pred_inf_determ_2(params1_valid, ini, dt, validDays_end-startDay+1, 1, R_valid);
    y_reg_valid = y_reg_valid(validDays_start:validDays_end,:);
    y_reg_test  = SEIR_data_gen_test_pred_inf_determ_2(params1_test, ini, dt, testDays_end-startDay+1, 1, R_test);
    y_reg_test  = y_reg_test(testDays_start:testDays_end,:);
    traj_train_wHwoM_mat(:, :, rep) = y_reg_train(:, 1: provNum);
    traj_valid_wHwoM_mat(:, :, rep) = y_reg_valid(:, 1: provNum);
    traj_test_wHwoM_mat(:, :, rep)  = y_reg_test(:, 1: provNum);
    %-------------------------
    trainErr_1 = abs(y_reg_train(:, 1: provNum) - y_true(1: trainDays_end-startDay+1, 1: provNum));
    trainErr_vec_1 = (sum(trainErr_1))' ./ (y_true(trainDays_end-startDay+1, (3*provNum+1):(4*provNum))  ...
        -   y_true(1, (3*provNum+1):(4*provNum))   )' ;
    validErr_1  = abs(y_reg_valid(:, 1: provNum) ...
        - y_true_valid(validDays_start: validDays_end, 1: provNum));
    validErr_vec_1 = (sum(validErr_1))' ./ (y_true_valid(validDays_end, (3*provNum+1):(4*provNum))  ...
        -   y_true_valid(validDays_start, (3*provNum+1):(4*provNum))   )' ;
    testErr_1  = abs(y_reg_test(:, 1: provNum) ...
        - y_true_test(testDays_start: testDays_end, 1: provNum));
    testErr_vec_1 = (sum(testErr_1))' ./ (y_true_test(testDays_end, (3*provNum+1):(4*provNum))  ...
        -   y_true_test(testDays_start, (3*provNum+1):(4*provNum))   )' ;
    trainErr_wHwoM_mat1(rep) = sum(trainErr_vec_1);
    validErr_wHwoM_mat1(rep)  = sum(validErr_vec_1);
    testErr_wHwoM_mat1(rep)  = sum(testErr_vec_1);
    %--------------------------
    trainErr_2  = abs(y_reg_train(:, 1: provNum) ...
        - y_true(1: trainDays_end-startDay+1, 1: provNum)) ...
        ./ y_true(1: trainDays_end-startDay+1, 1: provNum);
    trainErr_vec_2 = (sum(trainErr_2))' / trainDays_end ;
    validErr_2  = abs(y_reg_valid(:, 1: provNum) ...
        - y_true_valid(validDays_start: validDays_end, 1: provNum)) ...
        ./ y_true_valid(validDays_start: validDays_end, 1: provNum);
    validErr_vec_2 = (sum(validErr_2))' / (validDays_end - validDays_start + 1) ;
    testErr_2  = abs(y_reg_test(:, 1: provNum) ...
        - y_true_test(testDays_start: testDays_end, 1: provNum)) ...
        ./ y_true_test(testDays_start: testDays_end, 1: provNum);
    testErr_vec_2 = (sum(testErr_2))' / (testDays_end - testDays_start + 1) ;
    trainErr_wHwoM_mat2(rep) = sum(trainErr_vec_2);
    validErr_wHwoM_mat2(rep) = sum(validErr_vec_2);
    testErr_wHwoM_mat2(rep)  = sum(testErr_vec_2);
    %--------------------------
    trainErr_3  = ( abs(y_reg_train(:, 1: provNum) ...
        - y_true(1: trainDays_end-startDay+1, 1: provNum)) ...
        ./ y_true(1: trainDays_end-startDay+1, 1: provNum) ).^2;
    trainErr_vec_3 = sqrt((sum(trainErr_3 .* y_true(1: trainDays_end-startDay+1, 1: provNum)...
        ./ (ones(trainDays_end,1)*(y_true(trainDays_end-startDay+1, (3*provNum+1):(4*provNum))  ...
        -   y_true(1, (3*provNum+1):(4*provNum)))) ))' ) ;
    validErr_3  = ( abs(y_reg_valid(:, 1: provNum) ...
        - y_true_valid(validDays_start: validDays_end, 1: provNum)) ...
        ./ y_true_valid(validDays_start: validDays_end, 1: provNum) ).^2;
    validErr_vec_3 = sqrt((sum(validErr_3 .* y_true_valid(validDays_start: validDays_end, 1: provNum)...
        ./ (ones(validDays_end-validDays_start+1,1)*(y_true_valid(validDays_end, (3*provNum+1):(4*provNum))  ...
        -   y_true_valid(validDays_start, (3*provNum+1):(4*provNum)))) ))' ) ;
    testErr_3  = ( abs(y_reg_test(:, 1: provNum) ...
        - y_true_test(testDays_start: testDays_end, 1: provNum)) ...
        ./ y_true_test(testDays_start: testDays_end, 1: provNum) ).^2;
    testErr_vec_3 = sqrt((sum(testErr_3 .* y_true_test(testDays_start: testDays_end, 1: provNum)...
        ./ (ones(testDays_end-testDays_start+1,1)*(y_true_test(testDays_end, (3*provNum+1):(4*provNum))  ...
        -   y_true_test(testDays_start, (3*provNum+1):(4*provNum)))) ))' ) ;
    trainErr_wHwoM_mat3(rep) = sum(trainErr_vec_3);
    validErr_wHwoM_mat3(rep) = sum(validErr_vec_3);
    testErr_wHwoM_mat3(rep)  = sum(testErr_vec_3);
    %--------------------------
    trainErr_vec_4 = sqrt( (sum(trainErr_3))' / (trainDays_end-startDay+1) );
    validErr_vec_4 = sqrt( (sum(validErr_3))' / (validDays_end - validDays_start+1) );
    testErr_vec_4 = sqrt( (sum(testErr_3))' / (testDays_end-testDays_start+1) );
    trainErr_wHwoM_mat4(rep) = sum(trainErr_vec_4);
    validErr_wHwoM_mat4(rep) = sum(validErr_vec_4);
    testErr_wHwoM_mat4(rep) = sum(testErr_vec_4);
    
    %%
    %-----------------------------------------------------------------------------------------------------
    % w/ Hetero, w/ Migrat, no GL
    fprintf('w/ Hetero, w/ Migrat, no GL\n');
    ini.traMat = traMat; mu0 = 0; mu1 = 0;
    
    f1 = @(x)SEIR_fminunc_test_pred_inf_multi_2(x, y_obs, ini, 'pois', mu0, mu1, dt, trainDays_end, startDay, 0);
    opts1 = optimoptions('fmincon');
    opts1.MaxIterations =  2000;
    opts1.MaxFunctionEvaluations = 10^6;
    opts1.Display = 'off';%'iter-detailed';
    opts1.OptimalityTolerance = 1e-10;
    lb = [ones(1,2*provNum), zeros(1,provNum)];
    ub = [ones(1,2*provNum) * 100, ones(1,provNum) * 1];
    params_pre = [I_pre  * ones(1, length(provPop)), E_pre * ones(1, length(provPop)), lambda_pre*ones(1,length(provPop))];
    [params1, ~] = fmincon(f1, params_pre, [],[],[],[],lb,ub,[],opts1);
    param_wHwM_mat(:, rep) = params1';
    
    R_train = zeros(1, provNum);
    params1_test = params1;
    R_test  = zeros(1, provNum);
    params1_valid = params1;
    R_valid = zeros(1, provNum);
    y_reg_train = SEIR_data_gen_test_pred_inf_determ_2(params1, ini, dt, trainDays_end-startDay+1, 1, R_train);
    y_reg_valid = SEIR_data_gen_test_pred_inf_determ_2(params1_valid, ini, dt, validDays_end-startDay+1, 1, R_valid);
    y_reg_valid = y_reg_valid(validDays_start:validDays_end,:);
    y_reg_test  = SEIR_data_gen_test_pred_inf_determ_2(params1_test, ini, dt, testDays_end-startDay+1, 1, R_test);
    y_reg_test  = y_reg_test(testDays_start:testDays_end,:);
    traj_train_wHwM_mat(:, :, rep) = y_reg_train(:, 1: provNum);
    traj_valid_wHwM_mat(:, :, rep) = y_reg_valid(:, 1: provNum);
    traj_test_wHwM_mat(:, :, rep)  = y_reg_test(:, 1: provNum);
    %-------------------------
    trainErr_1 = abs(y_reg_train(:, 1: provNum) - y_true(1: trainDays_end-startDay+1, 1: provNum));
    trainErr_vec_1 = (sum(trainErr_1))' ./ (y_true(trainDays_end-startDay+1, (3*provNum+1):(4*provNum))  ...
        -   y_true(1, (3*provNum+1):(4*provNum))   )' ;
    validErr_1  = abs(y_reg_valid(:, 1: provNum) ...
        - y_true_valid(validDays_start: validDays_end, 1: provNum));
    validErr_vec_1 = (sum(validErr_1))' ./ (y_true_valid(validDays_end, (3*provNum+1):(4*provNum))  ...
        -   y_true_valid(validDays_start, (3*provNum+1):(4*provNum))   )' ;
    testErr_1  = abs(y_reg_test(:, 1: provNum) ...
        - y_true_test(testDays_start: testDays_end, 1: provNum));
    testErr_vec_1 = (sum(testErr_1))' ./ (y_true_test(testDays_end, (3*provNum+1):(4*provNum))  ...
        -   y_true_test(testDays_start, (3*provNum+1):(4*provNum))   )' ;
    trainErr_wHwM_mat1(rep) = sum(trainErr_vec_1);
    validErr_wHwM_mat1(rep)  = sum(validErr_vec_1);
    testErr_wHwM_mat1(rep)  = sum(testErr_vec_1);
    %--------------------------
    trainErr_2  = abs(y_reg_train(:, 1: provNum) ...
        - y_true(1: trainDays_end-startDay+1, 1: provNum)) ...
        ./ y_true(1: trainDays_end-startDay+1, 1: provNum);
    trainErr_vec_2 = (sum(trainErr_2))' / trainDays_end ;
    validErr_2  = abs(y_reg_valid(:, 1: provNum) ...
        - y_true_valid(validDays_start: validDays_end, 1: provNum)) ...
        ./ y_true_valid(validDays_start: validDays_end, 1: provNum);
    validErr_vec_2 = (sum(validErr_2))' / (validDays_end - validDays_start + 1) ;
    testErr_2  = abs(y_reg_test(:, 1: provNum) ...
        - y_true_test(testDays_start: testDays_end, 1: provNum)) ...
        ./ y_true_test(testDays_start: testDays_end, 1: provNum);
    testErr_vec_2 = (sum(testErr_2))' / (testDays_end - testDays_start + 1) ;
    trainErr_wHwM_mat2(rep) = sum(trainErr_vec_2);
    validErr_wHwM_mat2(rep) = sum(validErr_vec_2);
    testErr_wHwM_mat2(rep)  = sum(testErr_vec_2);
    %--------------------------
    trainErr_3  = ( abs(y_reg_train(:, 1: provNum) ...
        - y_true(1: trainDays_end-startDay+1, 1: provNum)) ...
        ./ y_true(1: trainDays_end-startDay+1, 1: provNum) ).^2;
    trainErr_vec_3 = sqrt((sum(trainErr_3 .* y_true(1: trainDays_end-startDay+1, 1: provNum)...
        ./ (ones(trainDays_end,1)*(y_true(trainDays_end-startDay+1, (3*provNum+1):(4*provNum))  ...
        -   y_true(1, (3*provNum+1):(4*provNum)))) ))' ) ;
    validErr_3  = ( abs(y_reg_valid(:, 1: provNum) ...
        - y_true_valid(validDays_start: validDays_end, 1: provNum)) ...
        ./ y_true_valid(validDays_start: validDays_end, 1: provNum) ).^2;
    validErr_vec_3 = sqrt((sum(validErr_3 .* y_true_valid(validDays_start: validDays_end, 1: provNum)...
        ./ (ones(validDays_end-validDays_start+1,1)*(y_true_valid(validDays_end, (3*provNum+1):(4*provNum))  ...
        -   y_true_valid(validDays_start, (3*provNum+1):(4*provNum)))) ))' ) ;
    testErr_3  = ( abs(y_reg_test(:, 1: provNum) ...
        - y_true_test(testDays_start: testDays_end, 1: provNum)) ...
        ./ y_true_test(testDays_start: testDays_end, 1: provNum) ).^2;
    testErr_vec_3 = sqrt((sum(testErr_3 .* y_true_test(testDays_start: testDays_end, 1: provNum)...
        ./ (ones(testDays_end-testDays_start+1,1)*(y_true_test(testDays_end, (3*provNum+1):(4*provNum))  ...
        -   y_true_test(testDays_start, (3*provNum+1):(4*provNum)))) ))' ) ;
    trainErr_wHwM_mat3(rep) = sum(trainErr_vec_3);
    validErr_wHwM_mat3(rep) = sum(validErr_vec_3);
    testErr_wHwM_mat3(rep)  = sum(testErr_vec_3);
    %--------------------------
    trainErr_vec_4 = sqrt( (sum(trainErr_3))' / (trainDays_end-startDay+1) );
    validErr_vec_4 = sqrt( (sum(validErr_3))' / (validDays_end - validDays_start+1) );
    testErr_vec_4 = sqrt( (sum(testErr_3))' / (testDays_end-testDays_start+1) );
    trainErr_wHwM_mat4(rep) = sum(trainErr_vec_4);
    validErr_wHwM_mat4(rep) = sum(validErr_vec_4);
    testErr_wHwM_mat4(rep) = sum(testErr_vec_4);
    
    %%
    
    %----------------------------------------------------------------------
    % w/ Hetero, w/ Migrat, w/ GL
    fprintf('w/ Hetero, w/ Migrat, w/ GL\n');
    ini.traMat = traMat;
    params_pre = [I_pre  * ones(1, length(provPop)), E_pre * ones(1, length(provPop)), lambda_pre*ones(1,length(provPop))];
    for i = 1: length(muvec)
        %for j = 1: length(mu1vec)
        fprintf('the %d-th mu1\n', i);
        
        %         mu = muvec(i); %mu1vec(j);
        %         f1 = @(x)SEIR_fminunc_test_pred_inf_multi_3(x, y_obs, ini, 'pois', mu, beta, dt, validDays, startDay, sigma );
        
        mu1 = muvec(i)/5000; %mu1vec(j);
        mu0 = mu1*beta;%mu1/10;
        f1 = @(x)SEIR_fminunc_test_pred_inf_multi_2(x, y_obs, ini, 'pois', mu0, mu1, dt, trainDays_end, startDay, sigma );
        
        opts1 = optimoptions('fmincon');
        opts1.MaxIterations =  2000;
        opts1.MaxFunctionEvaluations = 10^6;
        opts1.Display = 'off';%'off';%'iter-detailed';
        opts1.OptimalityTolerance = 1e-10;
        lb = [ones(1,2*provNum), zeros(1,provNum)];
        ub = [ones(1,2*provNum) * 100, ones(1,provNum) * 1];
        [params1, ~] = fmincon(f1, params_pre, [],[],[],[],lb,ub,[],opts1);
        %         if i == 1
        %             [params1, ~] = fmincon(f1, params_pre, [],[],[],[],lb,ub,[],opts1);
        %         else
        %             [params1, ~] = fmincon(f1, reshape(param_wHwMwGL_mat(i-1,:,rep),1,3*provNum), ...
        %                 [],[],[],[],lb,ub,[],opts1);
        %         end
        
        param_wHwMwGL_mat(i, :, rep) = params1';
        
        R_train = zeros(1, provNum);
        params1_test = params1;
        R_test  = zeros(1, provNum);
        params1_valid = params1;
        R_valid = zeros(1, provNum);
        y_reg_train = SEIR_data_gen_test_pred_inf_determ_2(params1, ini, dt, trainDays_end-startDay+1, 1, R_train);
        y_reg_valid = SEIR_data_gen_test_pred_inf_determ_2(params1_valid, ini, dt, validDays_end-startDay+1, 1, R_valid);
        y_reg_valid = y_reg_valid(validDays_start:validDays_end,:);
        y_reg_test  = SEIR_data_gen_test_pred_inf_determ_2(params1_test, ini, dt, testDays_end-startDay+1, 1, R_test);
        y_reg_test  = y_reg_test(testDays_start:testDays_end,:);
        %-------------------------
        trainErr_1 = abs(y_reg_train(:, 1: provNum) - y_true(1: trainDays_end-startDay+1, 1: provNum));
        trainErr_vec_1 = (sum(trainErr_1))' ./ (y_true(trainDays_end-startDay+1, (3*provNum+1):(4*provNum))  ...
            -   y_true(1, (3*provNum+1):(4*provNum))   )' ;
        validErr_1  = abs(y_reg_valid(:, 1: provNum) ...
            - y_true_valid(validDays_start: validDays_end, 1: provNum));
        validErr_vec_1 = (sum(validErr_1))' ./ (y_true_valid(validDays_end, (3*provNum+1):(4*provNum))  ...
            -   y_true_valid(validDays_start, (3*provNum+1):(4*provNum))   )' ;
        testErr_1  = abs(y_reg_test(:, 1: provNum) ...
            - y_true_test(testDays_start: testDays_end, 1: provNum));
        testErr_vec_1 = (sum(testErr_1))' ./ (y_true_test(testDays_end, (3*provNum+1):(4*provNum))  ...
            -   y_true_test(testDays_start, (3*provNum+1):(4*provNum))   )' ;
        trainErr_wHwMwGL_mat1(i,rep) = sum(trainErr_vec_1);
        validErr_wHwMwGL_mat1(i,rep)  = sum(validErr_vec_1);
        testErr_wHwMwGL_mat1(i,rep)  = sum(testErr_vec_1);
        %--------------------------
        trainErr_2  = abs(y_reg_train(:, 1: provNum) ...
            - y_true(1: trainDays_end-startDay+1, 1: provNum)) ...
            ./ y_true(1: trainDays_end-startDay+1, 1: provNum);
        trainErr_vec_2 = (sum(trainErr_2))' / trainDays_end ;
        validErr_2  = abs(y_reg_valid(:, 1: provNum) ...
            - y_true_valid(validDays_start: validDays_end, 1: provNum)) ...
            ./ y_true_valid(validDays_start: validDays_end, 1: provNum);
        validErr_vec_2 = (sum(validErr_2))' / (validDays_end - validDays_start + 1) ;
        testErr_2  = abs(y_reg_test(:, 1: provNum) ...
            - y_true_test(testDays_start: testDays_end, 1: provNum)) ...
            ./ y_true_test(testDays_start: testDays_end, 1: provNum);
        testErr_vec_2 = (sum(testErr_2))' / (testDays_end - testDays_start + 1) ;
        trainErr_wHwMwGL_mat2(i,rep) = sum(trainErr_vec_2);
        validErr_wHwMwGL_mat2(i,rep) = sum(validErr_vec_2);
        testErr_wHwMwGL_mat2(i,rep)  = sum(testErr_vec_2);
        %--------------------------
        trainErr_3  = ( abs(y_reg_train(:, 1: provNum) ...
            - y_true(1: trainDays_end-startDay+1, 1: provNum)) ...
            ./ y_true(1: trainDays_end-startDay+1, 1: provNum) ).^2;
        trainErr_vec_3 = sqrt((sum(trainErr_3 .* y_true(1: trainDays_end-startDay+1, 1: provNum)...
            ./ (ones(trainDays_end,1)*(y_true(trainDays_end-startDay+1, (3*provNum+1):(4*provNum))  ...
            -   y_true(1, (3*provNum+1):(4*provNum)))) ))' ) ;
        validErr_3  = ( abs(y_reg_valid(:, 1: provNum) ...
            - y_true_valid(validDays_start: validDays_end, 1: provNum)) ...
            ./ y_true_valid(validDays_start: validDays_end, 1: provNum) ).^2;
        validErr_vec_3 = sqrt((sum(validErr_3 .* y_true_valid(validDays_start: validDays_end, 1: provNum)...
            ./ (ones(validDays_end-validDays_start+1,1)*(y_true_valid(validDays_end, (3*provNum+1):(4*provNum))  ...
            -   y_true_valid(validDays_start, (3*provNum+1):(4*provNum)))) ))' ) ;
        testErr_3  = ( abs(y_reg_test(:, 1: provNum) ...
            - y_true_test(testDays_start: testDays_end, 1: provNum)) ...
            ./ y_true_test(testDays_start: testDays_end, 1: provNum) ).^2;
        testErr_vec_3 = sqrt((sum(testErr_3 .* y_true_test(testDays_start: testDays_end, 1: provNum)...
            ./ (ones(testDays_end-testDays_start+1,1)*(y_true_test(testDays_end, (3*provNum+1):(4*provNum))  ...
            -   y_true_test(testDays_start, (3*provNum+1):(4*provNum)))) ))' ) ;
        trainErr_wHwMwGL_mat3(i,rep) = sum(trainErr_vec_3);
        validErr_wHwMwGL_mat3(i,rep) = sum(validErr_vec_3);
        testErr_wHwMwGL_mat3(i,rep)  = sum(testErr_vec_3);
        %--------------------------
        trainErr_vec_4 = sqrt( (sum(trainErr_3))' / (trainDays_end-startDay+1) );
        validErr_vec_4 = sqrt( (sum(validErr_3))' / (validDays_end - validDays_start+1) );
        testErr_vec_4 = sqrt( (sum(testErr_3))' / (testDays_end-testDays_start+1) );
        trainErr_wHwMwGL_mat4(i,rep) = sum(trainErr_vec_4);
        validErr_wHwMwGL_mat4(i,rep) = sum(validErr_vec_4);
        testErr_wHwMwGL_mat4(i,rep) = sum(testErr_vec_4);
        
    end
    
    fprintf('\n')
    
end

toc;


return
%%

trainErr_woHwoM_mat1 = trainErr_woHwoM_mat1/provNum;
validErr_woHwoM_mat1 = validErr_woHwoM_mat1/provNum;
testErr_woHwoM_mat1  = testErr_woHwoM_mat1/provNum;
trainErr_woHwM_mat1  = trainErr_woHwM_mat1/provNum;
validErr_woHwM_mat1  = validErr_woHwM_mat1/provNum;
testErr_woHwM_mat1   = testErr_woHwM_mat1/provNum;
trainErr_wHwoM_mat1  = trainErr_wHwoM_mat1/provNum;
validErr_wHwoM_mat1  = validErr_wHwoM_mat1/provNum;
testErr_wHwoM_mat1   = testErr_wHwoM_mat1/provNum;
trainErr_wHwM_mat1   = trainErr_wHwM_mat1/provNum;
validErr_wHwM_mat1   = validErr_wHwM_mat1/provNum;
testErr_wHwM_mat1    = testErr_wHwM_mat1/provNum;
trainErr_wHwMwGL_mat1 = trainErr_wHwMwGL_mat1/provNum;
validErr_wHwMwGL_mat1 = validErr_wHwMwGL_mat1/provNum;
testErr_wHwMwGL_mat1  = testErr_wHwMwGL_mat1/provNum;

trainErr_woHwoM_mat2 = trainErr_woHwoM_mat2/provNum;
validErr_woHwoM_mat2 = validErr_woHwoM_mat2/provNum;
testErr_woHwoM_mat2  = testErr_woHwoM_mat2/provNum;
trainErr_woHwM_mat2  = trainErr_woHwM_mat2/provNum;
validErr_woHwM_mat2  = validErr_woHwM_mat2/provNum;
testErr_woHwM_mat2   = testErr_woHwM_mat2/provNum;
trainErr_wHwoM_mat2  = trainErr_wHwoM_mat2/provNum;
validErr_wHwoM_mat2  = validErr_wHwoM_mat2/provNum;
testErr_wHwoM_mat2   = testErr_wHwoM_mat2/provNum;
trainErr_wHwM_mat2   = trainErr_wHwM_mat2/provNum;
validErr_wHwM_mat2   = validErr_wHwM_mat2/provNum;
testErr_wHwM_mat2    = testErr_wHwM_mat2/provNum;
trainErr_wHwMwGL_mat2 = trainErr_wHwMwGL_mat2/provNum;
validErr_wHwMwGL_mat2 = validErr_wHwMwGL_mat2/provNum;
testErr_wHwMwGL_mat2  = testErr_wHwMwGL_mat2/provNum;

trainErr_woHwoM_mat3 = trainErr_woHwoM_mat3/provNum;
validErr_woHwoM_mat3 = validErr_woHwoM_mat3/provNum;
testErr_woHwoM_mat3  = testErr_woHwoM_mat3/provNum;
trainErr_woHwM_mat3  = trainErr_woHwM_mat3/provNum;
validErr_woHwM_mat3  = validErr_woHwM_mat3/provNum;
testErr_woHwM_mat3   = testErr_woHwM_mat3/provNum;
trainErr_wHwoM_mat3  = trainErr_wHwoM_mat3/provNum;
validErr_wHwoM_mat3  = validErr_wHwoM_mat3/provNum;
testErr_wHwoM_mat3   = testErr_wHwoM_mat3/provNum;
trainErr_wHwM_mat3   = trainErr_wHwM_mat3/provNum;
validErr_wHwM_mat3   = validErr_wHwM_mat3/provNum;
testErr_wHwM_mat3    = testErr_wHwM_mat3/provNum;
trainErr_wHwMwGL_mat3 = trainErr_wHwMwGL_mat3/provNum;
validErr_wHwMwGL_mat3 = validErr_wHwMwGL_mat3/provNum;
testErr_wHwMwGL_mat3  = testErr_wHwMwGL_mat3/provNum;

trainErr_woHwoM_mat4 = trainErr_woHwoM_mat4/provNum;
validErr_woHwoM_mat4 = validErr_woHwoM_mat4/provNum;
testErr_woHwoM_mat4  = testErr_woHwoM_mat4/provNum;
trainErr_woHwM_mat4  = trainErr_woHwM_mat4/provNum;
validErr_woHwM_mat4  = validErr_woHwM_mat4/provNum;
testErr_woHwM_mat4   = testErr_woHwM_mat4/provNum;
trainErr_wHwoM_mat4  = trainErr_wHwoM_mat4/provNum;
validErr_wHwoM_mat4  = validErr_wHwoM_mat4/provNum;
testErr_wHwoM_mat4   = testErr_wHwoM_mat4/provNum;
trainErr_wHwM_mat4   = trainErr_wHwM_mat4/provNum;
validErr_wHwM_mat4   = validErr_wHwM_mat4/provNum;
testErr_wHwM_mat4    = testErr_wHwM_mat4/provNum;
trainErr_wHwMwGL_mat4 = trainErr_wHwMwGL_mat4/provNum;
validErr_wHwMwGL_mat4 = validErr_wHwMwGL_mat4/provNum;
testErr_wHwMwGL_mat4  = testErr_wHwMwGL_mat4/provNum;


%%

muind_plot = 19;

figure(5), clf;
subplot(1,2,1);
hold on;
plot(log10(muvec), mean(validErr_wHwMwGL_mat1,2), 'x-', 'LineWidth', 2);
plot(log10(muvec(muind_plot)), mean(validErr_wHwMwGL_mat1(muind_plot,:),2), 'p', 'MarkerSize', 24, 'MarkerFaceColor', '#EDB120', 'MarkerEdgeColor', '#EDB120');
yline(mean(validErr_wHwM_mat1), '-', '$\mu=0$', 'LineWidth', 2, 'FontSize', 30, 'LabelHorizontalAlignment', 'left', 'Interpreter', 'latex');
xlabel('$\log_{10} \mu$', 'Interpreter', 'latex');
ylabel('validation error', 'Interpreter', 'latex');
set(gca, 'FontSize', 36);
title('$\textnormal{MAE}^{[\rm Val]}_{(w)}$', 'FontSize', 36, 'Interpreter', 'latex');
grid on;
subplot(1,2,2);
hold on;
plot(log10(muvec), mean(validErr_wHwMwGL_mat3,2), 'x-', 'LineWidth', 2);
plot(log10(muvec(muind_plot)), mean(validErr_wHwMwGL_mat3(muind_plot,:),2), 'p', 'MarkerSize', 24, 'MarkerFaceColor', '#EDB120', 'MarkerEdgeColor', '#EDB120');
yline(mean(validErr_wHwM_mat3), '-', '$\mu=0$', 'LineWidth', 2, 'FontSize', 30, 'LabelHorizontalAlignment', 'left', 'Interpreter', 'latex');
xlabel('$\log_{10} \mu$', 'Interpreter', 'latex');
ylabel('validation error', 'Interpreter', 'latex');
set(gca, 'FontSize', 36);
title('$\textnormal{MSE}^{[\rm Val]}_{(w)}$', 'FontSize', 36, 'Interpreter', 'latex');
grid on;

sgtitle({'Mean of total weighted relative validation errors over 100 replicas', '(simulated data with 30 provinces)'}, 'Fontsize', 42, 'Interpreter', 'latex');

figure(6), clf;
subplot(1,2,1);
hold on;
plot(log10(muvec), mean(testErr_wHwMwGL_mat1,2), 'x-', 'LineWidth', 2);
plot(log10(muvec(muind_plot)), mean(testErr_wHwMwGL_mat1(muind_plot,:),2), 'p', 'MarkerSize', 24, 'MarkerFaceColor', '#EDB120', 'MarkerEdgeColor', '#EDB120');
yline(mean(testErr_wHwM_mat1), '-', '$\mu=0$', 'LineWidth', 2, 'FontSize', 30, 'LabelHorizontalAlignment', 'left', 'Interpreter', 'latex');
xlabel('$\log_{10} \mu$', 'Interpreter', 'latex');
ylabel('testing error', 'Interpreter', 'latex');
set(gca, 'FontSize', 36);
title('$\textnormal{MAE}^{[\rm Te]}_{(w)}$', 'FontSize', 36, 'Interpreter', 'latex');
grid on;
subplot(1,2,2);
hold on;
plot(log10(muvec), mean(testErr_wHwMwGL_mat3,2), 'x-', 'LineWidth', 2);
plot(log10(muvec(muind_plot)), mean(testErr_wHwMwGL_mat3(muind_plot,:),2), 'p', 'MarkerSize', 24, 'MarkerFaceColor', '#EDB120', 'MarkerEdgeColor', '#EDB120');
yline(mean(testErr_wHwM_mat3), '-', '$\mu=0$', 'LineWidth', 2, 'FontSize', 30, 'LabelHorizontalAlignment', 'left', 'Interpreter', 'latex');
xlabel('$\log_{10} \mu$', 'Interpreter', 'latex');
ylabel('testing error', 'Interpreter', 'latex');
set(gca, 'FontSize', 36);
title('$\textnormal{MSE}^{[\rm Te]}_{(w)}$', 'FontSize', 36, 'Interpreter', 'latex');
grid on;

sgtitle({'Mean of total weighted relative testing errors over 100 replicas', '(simulated data with 30 provinces)'}, 'Fontsize', 42, 'Interpreter', 'latex');



figure(7), clf;
subplot(1,2,1);
hold on;
plot(log10(muvec), mean(validErr_wHwMwGL_mat2,2), 'x-', 'LineWidth', 2);
plot(log10(muvec(muind_plot)), mean(validErr_wHwMwGL_mat2(muind_plot,:),2), 'p', 'MarkerSize', 24, 'MarkerFaceColor', '#EDB120', 'MarkerEdgeColor', '#EDB120');
yline(mean(validErr_wHwM_mat2), '-', '$\mu=0$', 'LineWidth', 2, 'FontSize', 30, 'LabelHorizontalAlignment', 'left', 'Interpreter', 'latex');
xlabel('$\log_{10} \mu$', 'Interpreter', 'latex');
ylabel('validation error', 'Interpreter', 'latex');
set(gca, 'FontSize', 36);
title('$\textnormal{MAE}^{[\rm Val]}_{(s)}$', 'FontSize', 36, 'Interpreter', 'latex');
grid on;
subplot(1,2,2);
hold on;
plot(log10(muvec), mean(validErr_wHwMwGL_mat4,2), 'x-', 'LineWidth', 2);
plot(log10(muvec(muind_plot)), mean(validErr_wHwMwGL_mat4(muind_plot,:),2), 'p', 'MarkerSize', 24, 'MarkerFaceColor', '#EDB120', 'MarkerEdgeColor', '#EDB120');
yline(mean(validErr_wHwM_mat4), '-', '$\mu=0$', 'LineWidth', 2, 'FontSize', 30, 'LabelHorizontalAlignment', 'left', 'Interpreter', 'latex');
xlabel('$\log_{10} \mu$', 'Interpreter', 'latex');
ylabel('validation error', 'Interpreter', 'latex');
set(gca, 'FontSize', 36);
title('$\textnormal{MSE}^{[\rm Val]}_{(s)}$', 'FontSize', 36, 'Interpreter', 'latex');
grid on;

sgtitle({'Mean of total simple averaged relative validation errors over 100 replicas', '(simulated data with 30 provinces)'}, 'Fontsize', 42, 'Interpreter', 'latex');

figure(8), clf;
subplot(1,2,1);
hold on;
plot(log10(muvec), mean(testErr_wHwMwGL_mat2,2), 'x-', 'LineWidth', 2);
plot(log10(muvec(muind_plot)), mean(testErr_wHwMwGL_mat2(muind_plot,:),2), 'p', 'MarkerSize', 24, 'MarkerFaceColor', '#EDB120', 'MarkerEdgeColor', '#EDB120');
yline(mean(testErr_wHwM_mat2), '-', '$\mu=0$', 'LineWidth', 2, 'FontSize', 30, 'LabelHorizontalAlignment', 'left', 'Interpreter', 'latex');
xlabel('$\log_{10} \mu$', 'Interpreter', 'latex');
ylabel('testing error', 'Interpreter', 'latex');
% ylim([0.83 0.93])
set(gca, 'FontSize', 36);
title('$\textnormal{MAE}^{[\rm Te]}_{(s)}$', 'FontSize', 36, 'Interpreter', 'latex');
grid on;
subplot(1,2,2);
hold on;
plot(log10(muvec), mean(testErr_wHwMwGL_mat4,2), 'x-', 'LineWidth', 2);
plot(log10(muvec(muind_plot)), mean(testErr_wHwMwGL_mat4(muind_plot,:),2), 'p', 'MarkerSize', 24, 'MarkerFaceColor', '#EDB120', 'MarkerEdgeColor', '#EDB120');
yline(mean(testErr_wHwM_mat4), '-', '$\mu=0$', 'LineWidth', 2, 'FontSize', 30, 'LabelHorizontalAlignment', 'left', 'Interpreter', 'latex');
xlabel('$\log_{10} \mu$', 'Interpreter', 'latex');
ylabel('testing error', 'Interpreter', 'latex');
set(gca, 'FontSize', 36);
title('$\textnormal{MSE}^{[\rm Te]}_{(s)}$', 'FontSize', 36, 'Interpreter', 'latex');
grid on;

sgtitle({'Mean of total simply averaged relative testing errors over 100 replicas', '(simulated data with 30 provinces)'}, 'Fontsize', 42, 'Interpreter', 'latex');

return;


%%

trainDays_end
validDays_end
totDays

fprintf('relative error with weighted average \n\n');

fprintf('wo H, wo M, trainErr: mean %6.4f, std % 6.4f, testErr: mean %6.4f, std % 6.4f \n', mean(trainErr_woHwoM_mat1/provNum),...
    std(trainErr_woHwoM_mat1/provNum), mean(testErr_woHwoM_mat1/provNum), std(testErr_woHwoM_mat1/provNum));
fprintf('wo H, w/ M, trainErr: mean %6.4f, std % 6.4f, testErr: mean %6.4f, std % 6.4f \n', mean(trainErr_woHwM_mat1/provNum),...
    std(trainErr_woHwM_mat1/provNum), mean(testErr_woHwM_mat1/provNum), std(testErr_woHwM_mat1/provNum));
fprintf('w/ H, wo M, trainErr: mean %6.4f, std % 6.4f, testErr: mean %6.4f, std % 6.4f \n', mean(trainErr_wHwoM_mat1/provNum),...
    std(trainErr_wHwoM_mat1/provNum), mean(testErr_wHwoM_mat1/provNum), std(testErr_wHwoM_mat1/provNum));
fprintf('w/ H, w/ M, trainErr: mean %6.4f, std % 6.4f, testErr: mean %6.4f, std % 6.4f \n', mean(trainErr_wHwM_mat1/provNum),...
    std(trainErr_wHwM_mat1/provNum), mean(testErr_wHwM_mat1/provNum), std(testErr_wHwM_mat1/provNum));
fprintf('\n');
for i = 1: length(muvec)
    fprintf('w/ H, w/ M, w/ GL, mu = 10^%.1f, trainErr: mean %.4f, std %6.4f, testErr: mean %6.4f, std %.4f \n', log10(muvec(i)), mean(trainErr_wHwMwGL_mat1(i,:)/provNum),...
        std(trainErr_wHwMwGL_mat1(i,:)/provNum),  mean(testErr_wHwMwGL_mat1(i,:)/provNum), std(testErr_wHwMwGL_mat1(i,:)/provNum));
end

%%

fprintf('relative error with equal weights \n\n');

fprintf('wo H, wo M, trainErr: mean %6.4f, std % 6.4f, testErr: mean %6.4f, std % 6.4f \n', mean(trainErr_woHwoM_mat2/provNum),...
    std(trainErr_woHwoM_mat2/provNum), mean(testErr_woHwoM_mat2/provNum), std(testErr_woHwoM_mat2/provNum));
fprintf('wo H, w/ M, trainErr: mean %6.4f, std % 6.4f, testErr: mean %6.4f, std % 6.4f \n', mean(trainErr_woHwM_mat2/provNum),...
    std(trainErr_woHwM_mat2/provNum), mean(testErr_woHwM_mat2/provNum), std(testErr_woHwM_mat2/provNum));
fprintf('w/ H, wo M, trainErr: mean %6.4f, std % 6.4f, testErr: mean %6.4f, std % 6.4f \n', mean(trainErr_wHwoM_mat2/provNum),...
    std(trainErr_wHwoM_mat2/provNum), mean(testErr_wHwoM_mat2/provNum), std(testErr_wHwoM_mat2/provNum));
fprintf('w/ H, w/ M, trainErr: mean %6.4f, std % 6.4f, testErr: mean %6.4f, std % 6.4f \n', mean(trainErr_wHwM_mat2/provNum),...
    std(trainErr_wHwM_mat2/provNum), mean(testErr_wHwM_mat2/provNum), std(testErr_wHwM_mat2/provNum));
fprintf('\n');
for i = 1: length(muvec)
    %fprintf('mu0 = %6.4f: \n\n', mu0vec(i));
    %     for j = 1: length(mu1vec)
    %         fprintf('w/ H, w/ M, w/ GL, mu1 = %6.4f, trainErr: mean %6.4f, std % 6.4f, testErr: mean %6.4f, std % 6.4f \n', mu1vec(j), mean(trainErr_wHwMwGL_mat(i,j,:)),...
    %             std(trainErr_wHwMwGL_mat(i,j,:)),  mean(testErr_wHwMwGL_mat(i,j,:)), std(testErr_wHwMwGL_mat(i,j,:)));
    %     end
    
    fprintf('w/ H, w/ M, w/ GL, mu1 = %6.4f, trainErr: mean %6.4f, std % 6.4f, testErr: mean %6.4f, std % 6.4f \n', muvec(i), mean(trainErr_wHwMwGL_mat2(i,:)/provNum),...
        std(trainErr_wHwMwGL_mat2(i,:)/provNum),  mean(testErr_wHwMwGL_mat2(i,:)/provNum), std(testErr_wHwMwGL_mat2(i,:)/provNum));
    
end
%%
fprintf('\nrelative MSE with weighted average\n\n');

fprintf('wo H, wo M, trainErr: mean %6.4f, std % 6.4f, testErr: mean %6.4f, std % 6.4f \n', mean(trainErr_woHwoM_mat3/provNum),...
    std(trainErr_woHwoM_mat3/provNum), mean(testErr_woHwoM_mat3/provNum), std(testErr_woHwoM_mat3/provNum));
fprintf('wo H, w/ M, trainErr: mean %6.4f, std % 6.4f, testErr: mean %6.4f, std % 6.4f \n', mean(trainErr_woHwM_mat3/provNum),...
    std(trainErr_woHwM_mat3/provNum), mean(testErr_woHwM_mat3/provNum), std(testErr_woHwM_mat3/provNum));
fprintf('w/ H, wo M, trainErr: mean %6.4f, std % 6.4f, testErr: mean %6.4f, std % 6.4f \n', mean(trainErr_wHwoM_mat3/provNum),...
    std(trainErr_wHwoM_mat3/provNum), mean(testErr_wHwoM_mat3/provNum), std(testErr_wHwoM_mat3/provNum));
fprintf('w/ H, w/ M, trainErr: mean %6.4f, std % 6.4f, testErr: mean %6.4f, std % 6.4f \n', mean(trainErr_wHwM_mat3/provNum),...
    std(trainErr_wHwM_mat3/provNum), mean(testErr_wHwM_mat3/provNum), std(testErr_wHwM_mat3/provNum));
fprintf('\n');
for i = 1: length(muvec)
    %fprintf('mu0 = %6.4f: \n\n', mu0vec(i));
    %     for j = 1: length(mu1vec)
    %         fprintf('w/ H, w/ M, w/ GL, mu1 = %6.4f, trainErr: mean %6.4f, std % 6.4f, testErr: mean %6.4f, std % 6.4f \n', mu1vec(j), mean(trainErr_wHwMwGL_mat(i,j,:)),...
    %             std(trainErr_wHwMwGL_mat(i,j,:)),  mean(testErr_wHwMwGL_mat(i,j,:)), std(testErr_wHwMwGL_mat(i,j,:)));
    %     end
    
    fprintf('w/ H, w/ M, w/ GL, mu1 = 10^%.1f, trainErr: mean %6.4f, std % 6.4f, testErr: mean %6.4f, std % 6.4f \n', log10(muvec(i)), mean(trainErr_wHwMwGL_mat3(i,:)/provNum),...
        std(trainErr_wHwMwGL_mat3(i,:)/provNum),  mean(testErr_wHwMwGL_mat3(i,:)/provNum), std(testErr_wHwMwGL_mat3(i,:)/provNum));
    
end
%%

fprintf('\nrelative MSE with simple average\n\n');

fprintf('wo H, wo M, trainErr: mean %6.4f, std % 6.4f, testErr: mean %6.4f, std % 6.4f \n', mean(trainErr_woHwoM_mat4/provNum),...
    std(trainErr_woHwoM_mat4/provNum), mean(testErr_woHwoM_mat4/provNum), std(testErr_woHwoM_mat4/provNum));
fprintf('wo H, w/ M, trainErr: mean %6.4f, std % 6.4f, testErr: mean %6.4f, std % 6.4f \n', mean(trainErr_woHwM_mat4/provNum),...
    std(trainErr_woHwM_mat4/provNum), mean(testErr_woHwM_mat4/provNum), std(testErr_woHwM_mat4/provNum));
fprintf('w/ H, wo M, trainErr: mean %6.4f, std % 6.4f, testErr: mean %6.4f, std % 6.4f \n', mean(trainErr_wHwoM_mat4/provNum),...
    std(trainErr_wHwoM_mat4/provNum), mean(testErr_wHwoM_mat4/provNum), std(testErr_wHwoM_mat4/provNum));
fprintf('w/ H, w/ M, trainErr: mean %6.4f, std % 6.4f, testErr: mean %6.4f, std % 6.4f \n', mean(trainErr_wHwM_mat4/provNum),...
    std(trainErr_wHwM_mat4/provNum), mean(testErr_wHwM_mat4/provNum), std(testErr_wHwM_mat4/provNum));
fprintf('\n');
for i = 1: length(muvec)
    %fprintf('mu0 = %6.4f: \n\n', mu0vec(i));
    %     for j = 1: length(mu1vec)
    %         fprintf('w/ H, w/ M, w/ GL, mu1 = %6.4f, trainErr: mean %6.4f, std % 6.4f, testErr: mean %6.4f, std % 6.4f \n', mu1vec(j), mean(trainErr_wHwMwGL_mat(i,j,:)),...
    %             std(trainErr_wHwMwGL_mat(i,j,:)),  mean(testErr_wHwMwGL_mat(i,j,:)), std(testErr_wHwMwGL_mat(i,j,:)));
    %     end
    
    fprintf('w/ H, w/ M, w/ GL, mu1 = %6.4f, trainErr: mean %6.4f, std % 6.4f, testErr: mean %6.4f, std % 6.4f \n', muvec(i), mean(trainErr_wHwMwGL_mat4(i,:)/provNum),...
        std(trainErr_wHwMwGL_mat4(i,:)/provNum),  mean(testErr_wHwMwGL_mat4(i,:)/provNum), std(testErr_wHwMwGL_mat4(i,:)/provNum));
    
end


return


%%


fprintf('estimated lambda \n\n');

fprintf('wo H, wo M, lambda: mean %6.4f, std % 6.4f \n', mean(param_woHwoM_mat(2*provNum+1,:)), std(param_woHwoM_mat(2*provNum+1,:)));
fprintf('wo H, w/ M, lambda: mean %6.4f, std % 6.4f \n', mean(param_woHwM_mat(2*provNum+1,:)), std(param_woHwM_mat(2*provNum+1,:)));
%

prov_group1_ind = 4;
prov_group2_ind = 1;
prov_group3_ind = 3;


%-----------------------------------------------
fprintf('w/ H, wo M, lambda_k in prov %d in Group 1: mean %6.4f, std % 6.4f\n', ...
    prov_group1_ind,mean(param_wHwoM_mat(2*provNum+prov_group1_ind,:)), std(param_wHwoM_mat(2*provNum+prov_group1_ind,:)));
fprintf('            lambda_k in prov %d in Group 2: mean %6.4f, std % 6.4f\n', ...
    prov_group2_ind,mean(param_wHwoM_mat(2*provNum+prov_group2_ind,:)), std(param_wHwoM_mat(2*provNum+prov_group2_ind,:)));
fprintf('            lambda_k in prov %d in Group 3: mean %6.4f, std % 6.4f\n\n', ...
    prov_group3_ind,mean(param_wHwoM_mat(2*provNum+prov_group3_ind,:)), std(param_wHwoM_mat(2*provNum+prov_group3_ind,:)));

%-----------------------------------------------
fprintf('w/ H, w/ M, lambda_k in prov %d in Group 1: mean %6.4f, std % 6.4f\n', ...
    prov_group1_ind,mean(param_wHwM_mat(2*provNum+prov_group1_ind,:)), std(param_wHwM_mat(2*provNum+prov_group1_ind,:)));
fprintf('            lambda_k in prov %d in Group 2: mean %6.4f, std % 6.4f\n', ...
    prov_group2_ind,mean(param_wHwM_mat(2*provNum+prov_group2_ind,:)), std(param_wHwM_mat(2*provNum+prov_group2_ind,:)));
fprintf('            lambda_k in prov %d in Group 3: mean %6.4f, std % 6.4f\n\n', ...
    prov_group3_ind,mean(param_wHwM_mat(2*provNum+prov_group3_ind,:)), std(param_wHwM_mat(2*provNum+prov_group3_ind,:)));

%%
%-----------------------------------------------

fprintf('w/ H, w/ M, wGL: \n\n');
for i = 1: length(muvec)
    
    fprintf('%2d-th mu=10^%.1f, lambda_k in prov %d in Group 1: mean %6.4f, std % 6.4f\n', ...
        i, log10(muvec(i)),prov_group1_ind,mean(param_wHwMwGL_mat(i,2*provNum+prov_group1_ind,:)), std(param_wHwMwGL_mat(i,2*provNum+prov_group1_ind,:)));
    fprintf('                 lambda_k in prov %d in Group 2: mean %6.4f, std % 6.4f\n', ...
        prov_group2_ind,mean(param_wHwMwGL_mat(i,2*provNum+prov_group2_ind,:)), std(param_wHwMwGL_mat(i,2*provNum+prov_group2_ind,:)));
    fprintf('                 lambda_k in prov %d in Group 3: mean %6.4f, std % 6.4f\n\n', ...
        prov_group3_ind,mean(param_wHwMwGL_mat(i,2*provNum+prov_group3_ind,:)), std(param_wHwMwGL_mat(i,2*provNum+prov_group3_ind,:)));
end


%%







%           mismatched graph information










%% another group division for the same data

% mismatch of group
groupIndex2 = groupIndex;
group_ind_temp = 1:provNum;
group_ind_temp_1 = group_ind_temp(groupIndex==1);
group_ind_temp_2 = group_ind_temp(groupIndex==2);
group_ind_temp_3 = group_ind_temp(groupIndex==3);
% groupIndex2(group_ind_temp_1(2)) = 2;
% groupIndex2(group_ind_temp_2(2)) = 3;
% groupIndex2(group_ind_temp_3(2)) = 1;
groupIndex2(group_ind_temp_1(1)) = 2;
groupIndex2(group_ind_temp_2(1)) = 3;
groupIndex2(group_ind_temp_3(1)) = 1;

group_cell2 = cell(group_num,1);
for i = 1:group_num
    group_cell2{i} = group_ind_temp(groupIndex2 == i);
end

ini.part   = group_cell2;

group_cell2
%%

trainErr_wHwMwGL_mat1_group2 = zeros(length(muvec), numrep);
validErr_wHwMwGL_mat1_group2 = zeros(length(muvec), numrep);
testErr_wHwMwGL_mat1_group2  = zeros(length(muvec), numrep);

trainErr_wHwMwGL_mat2_group2 = zeros(length(muvec), numrep);
validErr_wHwMwGL_mat2_group2 = zeros(length(muvec), numrep);
testErr_wHwMwGL_mat2_group2  = zeros(length(muvec), numrep);

trainErr_wHwMwGL_mat3_group2 = zeros(length(muvec), numrep);
validErr_wHwMwGL_mat3_group2 = zeros(length(muvec), numrep);
testErr_wHwMwGL_mat3_group2  = zeros(length(muvec), numrep);

trainErr_wHwMwGL_mat4_group2 = zeros(length(muvec), numrep);
validErr_wHwMwGL_mat4_group2 = zeros(length(muvec), numrep);
testErr_wHwMwGL_mat4_group2  = zeros(length(muvec), numrep);

param_wHwMwGL_mat_group2 = zeros(length(muvec), 3*provNum, numrep);

% for rep = 1: numrep
for rep = 1: 100
    fprintf('the %d th replica\n', rep);
    
    
    y_true       = reshape(traj_obs_mat(:, :, rep), totDays, 5*provNum);
    y_true_test  = reshape(traj_obs_test_mat(:, :, rep), totDays, 5*provNum);
    y_true_valid = reshape(traj_obs_valid_mat(:, :, rep), totDays, 5*provNum);
    
    %%
    
    %----------------------------------------------------------------------
    % w/ Hetero, w/ Migrat, w/ GL
    fprintf('w/ Hetero, w/ Migrat, w/ GL\n');
    ini.traMat = traMat;
    params_pre = [I_pre  * ones(1, length(provPop)), E_pre * ones(1, length(provPop)), lambda_pre*ones(1,length(provPop))];
%     for i = 1: length(muvec)
    for i = 19
        %for j = 1: length(mu1vec)
        fprintf('the %d-th mu1 = 10^%.1f\n', i, log10(muvec(i)));
        
        %         mu = muvec(i); %mu1vec(j);
        %         f1 = @(x)SEIR_fminunc_test_pred_inf_multi_3(x, y_obs, ini, 'pois', mu, beta, dt, validDays, startDay, sigma );
        
        mu1 = muvec(i)/5000; %mu1vec(j);
        mu0 = mu1*beta;%mu1/10;
        f1 = @(x)SEIR_fminunc_test_pred_inf_multi_2(x, y_true, ini, 'pois', mu0, mu1, dt, trainDays_end, startDay, sigma );
        
        opts1 = optimoptions('fmincon');
        opts1.MaxIterations =  2000;
        opts1.MaxFunctionEvaluations = 10^6;
        opts1.Display = 'off';%'off';%'iter-detailed';
        opts1.OptimalityTolerance = 1e-10;
        lb = [ones(1,2*provNum), zeros(1,provNum)];
        ub = [ones(1,2*provNum) * 100, ones(1,provNum) * 1];
        tic
        [params1, ~] = fmincon(f1, params_pre, [],[],[],[],lb,ub,[],opts1);
        toc;
        %         if i == 1
        %             [params1, ~] = fmincon(f1, params_pre, [],[],[],[],lb,ub,[],opts1);
        %         else
        %             [params1, ~] = fmincon(f1, reshape(param_wHwMwGL_mat(i-1,:,rep),1,3*provNum), ...
        %                 [],[],[],[],lb,ub,[],opts1);
        %         end
        
        param_wHwMwGL_mat_group2(i, :, rep) = params1';
        
        R_train = zeros(1, provNum);
        params1_test = params1;
        R_test  = zeros(1, provNum);
        params1_valid = params1;
        R_valid = zeros(1, provNum);
        y_reg_train = SEIR_data_gen_test_pred_inf_determ_2(params1, ini, dt, trainDays_end-startDay+1, 1, R_train);
        y_reg_valid = SEIR_data_gen_test_pred_inf_determ_2(params1_valid, ini, dt, validDays_end-startDay+1, 1, R_valid);
        y_reg_valid = y_reg_valid(validDays_start:validDays_end,:);
        y_reg_test  = SEIR_data_gen_test_pred_inf_determ_2(params1_test, ini, dt, testDays_end-startDay+1, 1, R_test);
        y_reg_test  = y_reg_test(testDays_start:testDays_end,:);
        %-------------------------
        trainErr_1 = abs(y_reg_train(:, 1: provNum) - y_true(1: trainDays_end-startDay+1, 1: provNum));
        trainErr_vec_1 = (sum(trainErr_1))' ./ (y_true(trainDays_end-startDay+1, (3*provNum+1):(4*provNum))  ...
            -   y_true(1, (3*provNum+1):(4*provNum))   )' ;
        validErr_1  = abs(y_reg_valid(:, 1: provNum) ...
            - y_true_valid(validDays_start: validDays_end, 1: provNum));
        validErr_vec_1 = (sum(validErr_1))' ./ (y_true_valid(validDays_end, (3*provNum+1):(4*provNum))  ...
            -   y_true_valid(validDays_start, (3*provNum+1):(4*provNum))   )' ;
        testErr_1  = abs(y_reg_test(:, 1: provNum) ...
            - y_true_test(testDays_start: testDays_end, 1: provNum));
        testErr_vec_1 = (sum(testErr_1))' ./ (y_true_test(testDays_end, (3*provNum+1):(4*provNum))  ...
            -   y_true_test(testDays_start, (3*provNum+1):(4*provNum))   )' ;
        trainErr_wHwMwGL_mat1_group2(i,rep) = sum(trainErr_vec_1);
        validErr_wHwMwGL_mat1_group2(i,rep)  = sum(validErr_vec_1);
        testErr_wHwMwGL_mat1_group2(i,rep)  = sum(testErr_vec_1);
        %--------------------------
        trainErr_2  = abs(y_reg_train(:, 1: provNum) ...
            - y_true(1: trainDays_end-startDay+1, 1: provNum)) ...
            ./ y_true(1: trainDays_end-startDay+1, 1: provNum);
        trainErr_vec_2 = (sum(trainErr_2))' / trainDays_end ;
        validErr_2  = abs(y_reg_valid(:, 1: provNum) ...
            - y_true_valid(validDays_start: validDays_end, 1: provNum)) ...
            ./ y_true_valid(validDays_start: validDays_end, 1: provNum);
        validErr_vec_2 = (sum(validErr_2))' / (validDays_end - validDays_start + 1) ;
        testErr_2  = abs(y_reg_test(:, 1: provNum) ...
            - y_true_test(testDays_start: testDays_end, 1: provNum)) ...
            ./ y_true_test(testDays_start: testDays_end, 1: provNum);
        testErr_vec_2 = (sum(testErr_2))' / (testDays_end - testDays_start + 1) ;
        trainErr_wHwMwGL_mat2_group2(i,rep) = sum(trainErr_vec_2);
        validErr_wHwMwGL_mat2_group2(i,rep) = sum(validErr_vec_2);
        testErr_wHwMwGL_mat2_group2(i,rep)  = sum(testErr_vec_2);
        %--------------------------
        trainErr_3  = ( abs(y_reg_train(:, 1: provNum) ...
            - y_true(1: trainDays_end-startDay+1, 1: provNum)) ...
            ./ y_true(1: trainDays_end-startDay+1, 1: provNum) ).^2;
        trainErr_vec_3 = sqrt((sum(trainErr_3 .* y_true(1: trainDays_end-startDay+1, 1: provNum)...
            ./ (ones(trainDays_end,1)*(y_true(trainDays_end-startDay+1, (3*provNum+1):(4*provNum))  ...
            -   y_true(1, (3*provNum+1):(4*provNum)))) ))' ) ;
        validErr_3  = ( abs(y_reg_valid(:, 1: provNum) ...
            - y_true_valid(validDays_start: validDays_end, 1: provNum)) ...
            ./ y_true_valid(validDays_start: validDays_end, 1: provNum) ).^2;
        validErr_vec_3 = sqrt((sum(validErr_3 .* y_true_valid(validDays_start: validDays_end, 1: provNum)...
            ./ (ones(validDays_end-validDays_start+1,1)*(y_true_valid(validDays_end, (3*provNum+1):(4*provNum))  ...
            -   y_true_valid(validDays_start, (3*provNum+1):(4*provNum)))) ))' ) ;
        testErr_3  = ( abs(y_reg_test(:, 1: provNum) ...
            - y_true_test(testDays_start: testDays_end, 1: provNum)) ...
            ./ y_true_test(testDays_start: testDays_end, 1: provNum) ).^2;
        testErr_vec_3 = sqrt((sum(testErr_3 .* y_true_test(testDays_start: testDays_end, 1: provNum)...
            ./ (ones(testDays_end-testDays_start+1,1)*(y_true_test(testDays_end, (3*provNum+1):(4*provNum))  ...
            -   y_true_test(testDays_start, (3*provNum+1):(4*provNum)))) ))' ) ;
        trainErr_wHwMwGL_mat3_group2(i,rep) = sum(trainErr_vec_3);
        validErr_wHwMwGL_mat3_group2(i,rep) = sum(validErr_vec_3);
        testErr_wHwMwGL_mat3_group2(i,rep)  = sum(testErr_vec_3);
        %--------------------------
        trainErr_vec_4 = sqrt( (sum(trainErr_3))' / (trainDays_end-startDay+1) );
        validErr_vec_4 = sqrt( (sum(validErr_3))' / (validDays_end - validDays_start+1) );
        testErr_vec_4 = sqrt( (sum(testErr_3))' / (testDays_end-testDays_start+1) );
        trainErr_wHwMwGL_mat4_group2(i,rep) = sum(trainErr_vec_4);
        validErr_wHwMwGL_mat4_group2(i,rep) = sum(validErr_vec_4);
        testErr_wHwMwGL_mat4_group2(i,rep) = sum(testErr_vec_4);
        
    end
    
    fprintf('\n')
    
    
end


ini.part   = group_cell;

%% import data (specific code for the files 30prov_valid_group_1_25.mat etc)


rep_data_1_25 = load("30prov_valid_group_1_25.mat", "trainErr_wHwMwGL_mat1_group2", "trainErr_wHwMwGL_mat2_group2", "trainErr_wHwMwGL_mat3_group2", "trainErr_wHwMwGL_mat4_group2",...
    "validErr_wHwMwGL_mat1_group2", "validErr_wHwMwGL_mat2_group2", "validErr_wHwMwGL_mat3_group2", "validErr_wHwMwGL_mat4_group2",...
    "testErr_wHwMwGL_mat1_group2", "testErr_wHwMwGL_mat2_group2", "testErr_wHwMwGL_mat3_group2", "testErr_wHwMwGL_mat4_group2",...
    "param_wHwMwGL_mat_group2");
rep_data_26_50 = load("30prov_valid_group_26_50.mat", "trainErr_wHwMwGL_mat1_group2", "trainErr_wHwMwGL_mat2_group2", "trainErr_wHwMwGL_mat3_group2", "trainErr_wHwMwGL_mat4_group2",...
    "validErr_wHwMwGL_mat1_group2", "validErr_wHwMwGL_mat2_group2", "validErr_wHwMwGL_mat3_group2", "validErr_wHwMwGL_mat4_group2",...
    "testErr_wHwMwGL_mat1_group2", "testErr_wHwMwGL_mat2_group2", "testErr_wHwMwGL_mat3_group2", "testErr_wHwMwGL_mat4_group2",...
    "param_wHwMwGL_mat_group2");
rep_data_51_75 = load("30prov_valid_group_51_75.mat", "trainErr_wHwMwGL_mat1_group2", "trainErr_wHwMwGL_mat2_group2", "trainErr_wHwMwGL_mat3_group2", "trainErr_wHwMwGL_mat4_group2",...
    "validErr_wHwMwGL_mat1_group2", "validErr_wHwMwGL_mat2_group2", "validErr_wHwMwGL_mat3_group2", "validErr_wHwMwGL_mat4_group2",...
    "testErr_wHwMwGL_mat1_group2", "testErr_wHwMwGL_mat2_group2", "testErr_wHwMwGL_mat3_group2", "testErr_wHwMwGL_mat4_group2",...
    "param_wHwMwGL_mat_group2");
rep_data_76_100 = load("30prov_valid_group_76_100.mat", "trainErr_wHwMwGL_mat1_group2", "trainErr_wHwMwGL_mat2_group2", "trainErr_wHwMwGL_mat3_group2", "trainErr_wHwMwGL_mat4_group2",...
    "validErr_wHwMwGL_mat1_group2", "validErr_wHwMwGL_mat2_group2", "validErr_wHwMwGL_mat3_group2", "validErr_wHwMwGL_mat4_group2",...
    "testErr_wHwMwGL_mat1_group2", "testErr_wHwMwGL_mat2_group2", "testErr_wHwMwGL_mat3_group2", "testErr_wHwMwGL_mat4_group2",...
    "param_wHwMwGL_mat_group2");

rep_data = cell(4,1);
rep_data{1} = rep_data_1_25;
rep_data{2} = rep_data_26_50;
rep_data{3} = rep_data_51_75;
rep_data{4} = rep_data_76_100;

trainErr_wHwMwGL_mat1_group2 = zeros(length(muvec), numrep);
validErr_wHwMwGL_mat1_group2 = zeros(length(muvec), numrep);
testErr_wHwMwGL_mat1_group2  = zeros(length(muvec), numrep);

trainErr_wHwMwGL_mat2_group2 = zeros(length(muvec), numrep);
validErr_wHwMwGL_mat2_group2 = zeros(length(muvec), numrep);
testErr_wHwMwGL_mat2_group2  = zeros(length(muvec), numrep);

trainErr_wHwMwGL_mat3_group2 = zeros(length(muvec), numrep);
validErr_wHwMwGL_mat3_group2 = zeros(length(muvec), numrep);
testErr_wHwMwGL_mat3_group2  = zeros(length(muvec), numrep);

trainErr_wHwMwGL_mat4_group2 = zeros(length(muvec), numrep);
validErr_wHwMwGL_mat4_group2 = zeros(length(muvec), numrep);
testErr_wHwMwGL_mat4_group2  = zeros(length(muvec), numrep);

param_wHwMwGL_mat_group2 = zeros(length(muvec), 3*provNum, numrep);

for i = 0:3
    trainErr_wHwMwGL_mat1_group2(:,(25*i+1):25*(i+1)) = rep_data{i+1}.trainErr_wHwMwGL_mat1_group2(:,(25*i+1):25*(i+1));
    trainErr_wHwMwGL_mat2_group2(:,(25*i+1):25*(i+1)) = rep_data{i+1}.trainErr_wHwMwGL_mat2_group2(:,(25*i+1):25*(i+1));
    trainErr_wHwMwGL_mat3_group2(:,(25*i+1):25*(i+1)) = rep_data{i+1}.trainErr_wHwMwGL_mat3_group2(:,(25*i+1):25*(i+1));
    trainErr_wHwMwGL_mat4_group2(:,(25*i+1):25*(i+1)) = rep_data{i+1}.trainErr_wHwMwGL_mat4_group2(:,(25*i+1):25*(i+1));
    
    validErr_wHwMwGL_mat1_group2(:,(25*i+1):25*(i+1)) = rep_data{i+1}.validErr_wHwMwGL_mat1_group2(:,(25*i+1):25*(i+1));
    validErr_wHwMwGL_mat2_group2(:,(25*i+1):25*(i+1)) = rep_data{i+1}.validErr_wHwMwGL_mat2_group2(:,(25*i+1):25*(i+1));
    validErr_wHwMwGL_mat3_group2(:,(25*i+1):25*(i+1)) = rep_data{i+1}.validErr_wHwMwGL_mat3_group2(:,(25*i+1):25*(i+1));
    validErr_wHwMwGL_mat4_group2(:,(25*i+1):25*(i+1)) = rep_data{i+1}.validErr_wHwMwGL_mat4_group2(:,(25*i+1):25*(i+1));
    
    testErr_wHwMwGL_mat1_group2(:,(25*i+1):25*(i+1))  = rep_data{i+1}.testErr_wHwMwGL_mat1_group2(:,(25*i+1):25*(i+1));
    testErr_wHwMwGL_mat2_group2(:,(25*i+1):25*(i+1))  = rep_data{i+1}.testErr_wHwMwGL_mat2_group2(:,(25*i+1):25*(i+1));
    testErr_wHwMwGL_mat3_group2(:,(25*i+1):25*(i+1))  = rep_data{i+1}.testErr_wHwMwGL_mat3_group2(:,(25*i+1):25*(i+1));
    testErr_wHwMwGL_mat4_group2(:,(25*i+1):25*(i+1))  = rep_data{i+1}.testErr_wHwMwGL_mat4_group2(:,(25*i+1):25*(i+1));
    
    param_wHwMwGL_mat_group2(:,:,(25*i+1):25*(i+1)) = rep_data{i+1}.param_wHwMwGL_mat_group2(:,:,(25*i+1):25*(i+1));
    
end


clear rep_data rep_data_1_25 rep_data_26_50 rep_data_51_75 rep_data_76_100


%%

trainDays_end
validDays_end
totDays

fprintf('relative error with weighted average \n\n');

fprintf('\n');
for i = 1: length(muvec)
    fprintf('w/ H, w/ M, w/ GL, mu = 10^%.1f, trainErr: mean %.4f, std %6.4f, testErr: mean %6.4f, std %.4f \n', log10(muvec(i)), mean(trainErr_wHwMwGL_mat1_group2(i,:)/provNum),...
        std(trainErr_wHwMwGL_mat1_group2(i,:)/provNum),  mean(testErr_wHwMwGL_mat1_group2(i,:)/provNum), std(testErr_wHwMwGL_mat1_group2(i,:)/provNum));
end

%%

fprintf('relative error with equal weights \n\n');

fprintf('\n');
for i = 1: length(muvec)
    %fprintf('mu0 = %6.4f: \n\n', mu0vec(i));
    %     for j = 1: length(mu1vec)
    %         fprintf('w/ H, w/ M, w/ GL, mu1 = %6.4f, trainErr: mean %6.4f, std % 6.4f, testErr: mean %6.4f, std % 6.4f \n', mu1vec(j), mean(trainErr_wHwMwGL_mat(i,j,:)),...
    %             std(trainErr_wHwMwGL_mat(i,j,:)),  mean(testErr_wHwMwGL_mat(i,j,:)), std(testErr_wHwMwGL_mat(i,j,:)));
    %     end
    
    fprintf('w/ H, w/ M, w/ GL, mu1 = %6.4f, trainErr: mean %6.4f, std % 6.4f, testErr: mean %6.4f, std % 6.4f \n', muvec(i), mean(trainErr_wHwMwGL_mat2_group2(i,:)/provNum),...
        std(trainErr_wHwMwGL_mat2_group2(i,:)/provNum),  mean(testErr_wHwMwGL_mat2_group2(i,:)/provNum), std(testErr_wHwMwGL_mat2_group2(i,:)/provNum));
    
end
%%
fprintf('\nrelative MSE with weighted average\n\n');

fprintf('\n');
for i = 1: length(muvec)
    %fprintf('mu0 = %6.4f: \n\n', mu0vec(i));
    %     for j = 1: length(mu1vec)
    %         fprintf('w/ H, w/ M, w/ GL, mu1 = %6.4f, trainErr: mean %6.4f, std % 6.4f, testErr: mean %6.4f, std % 6.4f \n', mu1vec(j), mean(trainErr_wHwMwGL_mat(i,j,:)),...
    %             std(trainErr_wHwMwGL_mat(i,j,:)),  mean(testErr_wHwMwGL_mat(i,j,:)), std(testErr_wHwMwGL_mat(i,j,:)));
    %     end
    
    fprintf('w/ H, w/ M, w/ GL, mu1 = 10^%.1f, trainErr: mean %6.4f, std % 6.4f, testErr: mean %6.4f, std % 6.4f \n', log10(muvec(i)), mean(trainErr_wHwMwGL_mat3_group2(i,:)/provNum),...
        std(trainErr_wHwMwGL_mat3_group2(i,:)/provNum),  mean(testErr_wHwMwGL_mat3_group2(i,:)/provNum), std(testErr_wHwMwGL_mat3_group2(i,:)/provNum));
    
end
%%

fprintf('\nrelative MSE with simple average\n\n');

fprintf('\n');
for i = 1: length(muvec)
    %fprintf('mu0 = %6.4f: \n\n', mu0vec(i));
    %     for j = 1: length(mu1vec)
    %         fprintf('w/ H, w/ M, w/ GL, mu1 = %6.4f, trainErr: mean %6.4f, std % 6.4f, testErr: mean %6.4f, std % 6.4f \n', mu1vec(j), mean(trainErr_wHwMwGL_mat(i,j,:)),...
    %             std(trainErr_wHwMwGL_mat(i,j,:)),  mean(testErr_wHwMwGL_mat(i,j,:)), std(testErr_wHwMwGL_mat(i,j,:)));
    %     end
    
    fprintf('w/ H, w/ M, w/ GL, mu1 = %6.4f, trainErr: mean %6.4f, std % 6.4f, testErr: mean %6.4f, std % 6.4f \n', muvec(i), mean(trainErr_wHwMwGL_mat4_group2(i,:)/provNum),...
        std(trainErr_wHwMwGL_mat4_group2(i,:)/provNum),  mean(testErr_wHwMwGL_mat4_group2(i,:)/provNum), std(testErr_wHwMwGL_mat4_group2(i,:)/provNum));
    
end


return


%%


fprintf('estimated lambda \n\n');

prov_group1_ind = 4;
prov_group2_ind = 1;
prov_group3_ind = 3;
%-----------------------------------------------

fprintf('w/ H, w/ M, wGL: \n\n');
for i = 1: length(muvec)
    
    fprintf('%2d-th mu=10^%.1f, lambda_k in prov %d in Group 1: mean %6.4f, std % 6.4f\n', ...
        i, log10(muvec(i)),prov_group1_ind,mean(param_wHwMwGL_mat_group2(i,2*provNum+prov_group1_ind,:)), std(param_wHwMwGL_mat_group2(i,2*provNum+prov_group1_ind,:)));
    fprintf('                 lambda_k in prov %d in Group 2: mean %6.4f, std % 6.4f\n', ...
        prov_group2_ind,mean(param_wHwMwGL_mat_group2(i,2*provNum+prov_group2_ind,:)), std(param_wHwMwGL_mat_group2(i,2*provNum+prov_group2_ind,:)));
    fprintf('                 lambda_k in prov %d in Group 3: mean %6.4f, std % 6.4f\n\n', ...
        prov_group3_ind,mean(param_wHwMwGL_mat_group2(i,2*provNum+prov_group3_ind,:)), std(param_wHwMwGL_mat_group2(i,2*provNum+prov_group3_ind,:)));
end



%%









%           mismatched graph information in another way








%% a third group division for the same data

% mismatch of group
groupIndex3 = groupIndex;
group_ind_temp = 1:provNum;
group_ind_temp_1 = group_ind_temp(groupIndex==1);
group_ind_temp_2 = group_ind_temp(groupIndex==2);
group_ind_temp_3 = group_ind_temp(groupIndex==3);
groupIndex3(group_ind_temp_1(2)) = 2;
groupIndex3(group_ind_temp_2(2)) = 3;
groupIndex3(group_ind_temp_3(2)) = 1;
% groupIndex2(group_ind_temp_1(1)) = 2;
% groupIndex2(group_ind_temp_2(1)) = 3;
% groupIndex2(group_ind_temp_3(1)) = 1;

group_cell3 = cell(group_num,1);
for i = 1:group_num
    group_cell3{i} = group_ind_temp(groupIndex3 == i);
end

ini.part   = group_cell3;

group_cell3
%

trainErr_wHwMwGL_mat1_group3 = zeros(length(muvec), numrep);
validErr_wHwMwGL_mat1_group3 = zeros(length(muvec), numrep);
testErr_wHwMwGL_mat1_group3  = zeros(length(muvec), numrep);

trainErr_wHwMwGL_mat2_group3 = zeros(length(muvec), numrep);
validErr_wHwMwGL_mat2_group3 = zeros(length(muvec), numrep);
testErr_wHwMwGL_mat2_group3  = zeros(length(muvec), numrep);

trainErr_wHwMwGL_mat3_group3 = zeros(length(muvec), numrep);
validErr_wHwMwGL_mat3_group3 = zeros(length(muvec), numrep);
testErr_wHwMwGL_mat3_group3  = zeros(length(muvec), numrep);

trainErr_wHwMwGL_mat4_group3 = zeros(length(muvec), numrep);
validErr_wHwMwGL_mat4_group3 = zeros(length(muvec), numrep);
testErr_wHwMwGL_mat4_group3  = zeros(length(muvec), numrep);

param_wHwMwGL_mat_group3 = zeros(length(muvec), 3*provNum, numrep);

% for rep = 1: numrep
for rep = 1: 100
    fprintf('the %d th replica\n', rep);
    
    
    y_true       = reshape(traj_obs_mat(:, :, rep), totDays, 5*provNum);
    y_true_test  = reshape(traj_obs_test_mat(:, :, rep), totDays, 5*provNum);
    y_true_valid = reshape(traj_obs_valid_mat(:, :, rep), totDays, 5*provNum);
    
    %%
    
    %----------------------------------------------------------------------
    % w/ Hetero, w/ Migrat, w/ GL
    fprintf('w/ Hetero, w/ Migrat, w/ GL\n');
    ini.traMat = traMat;
    params_pre = [I_pre  * ones(1, length(provPop)), E_pre * ones(1, length(provPop)), lambda_pre*ones(1,length(provPop))];
%     for i = 1: length(muvec)
    for i = 19
        %for j = 1: length(mu1vec)
        fprintf('the %d-th mu1 = 10^%.1f\n', i, log10(muvec(i)));
        
        %         mu = muvec(i); %mu1vec(j);
        %         f1 = @(x)SEIR_fminunc_test_pred_inf_multi_3(x, y_obs, ini, 'pois', mu, beta, dt, validDays, startDay, sigma );
        
        mu1 = muvec(i)/5000; %mu1vec(j);
        mu0 = mu1*beta;%mu1/10;
        f1 = @(x)SEIR_fminunc_test_pred_inf_multi_2(x, y_true, ini, 'pois', mu0, mu1, dt, trainDays_end, startDay, sigma );
        
        opts1 = optimoptions('fmincon');
        opts1.MaxIterations =  2000;
        opts1.MaxFunctionEvaluations = 10^6;
        opts1.Display = 'off';%'off';%'iter-detailed';
        opts1.OptimalityTolerance = 1e-10;
        lb = [ones(1,2*provNum), zeros(1,provNum)];
        ub = [ones(1,2*provNum) * 100, ones(1,provNum) * 1];
        tic
        [params1, ~] = fmincon(f1, params_pre, [],[],[],[],lb,ub,[],opts1);
        toc;
        %         if i == 1
        %             [params1, ~] = fmincon(f1, params_pre, [],[],[],[],lb,ub,[],opts1);
        %         else
        %             [params1, ~] = fmincon(f1, reshape(param_wHwMwGL_mat(i-1,:,rep),1,3*provNum), ...
        %                 [],[],[],[],lb,ub,[],opts1);
        %         end
        
        param_wHwMwGL_mat_group3(i, :, rep) = params1';
        
        R_train = zeros(1, provNum);
        params1_test = params1;
        R_test  = zeros(1, provNum);
        params1_valid = params1;
        R_valid = zeros(1, provNum);
        y_reg_train = SEIR_data_gen_test_pred_inf_determ_2(params1, ini, dt, trainDays_end-startDay+1, 1, R_train);
        y_reg_valid = SEIR_data_gen_test_pred_inf_determ_2(params1_valid, ini, dt, validDays_end-startDay+1, 1, R_valid);
        y_reg_valid = y_reg_valid(validDays_start:validDays_end,:);
        y_reg_test  = SEIR_data_gen_test_pred_inf_determ_2(params1_test, ini, dt, testDays_end-startDay+1, 1, R_test);
        y_reg_test  = y_reg_test(testDays_start:testDays_end,:);
        %-------------------------
        trainErr_1 = abs(y_reg_train(:, 1: provNum) - y_true(1: trainDays_end-startDay+1, 1: provNum));
        trainErr_vec_1 = (sum(trainErr_1))' ./ (y_true(trainDays_end-startDay+1, (3*provNum+1):(4*provNum))  ...
            -   y_true(1, (3*provNum+1):(4*provNum))   )' ;
        validErr_1  = abs(y_reg_valid(:, 1: provNum) ...
            - y_true_valid(validDays_start: validDays_end, 1: provNum));
        validErr_vec_1 = (sum(validErr_1))' ./ (y_true_valid(validDays_end, (3*provNum+1):(4*provNum))  ...
            -   y_true_valid(validDays_start, (3*provNum+1):(4*provNum))   )' ;
        testErr_1  = abs(y_reg_test(:, 1: provNum) ...
            - y_true_test(testDays_start: testDays_end, 1: provNum));
        testErr_vec_1 = (sum(testErr_1))' ./ (y_true_test(testDays_end, (3*provNum+1):(4*provNum))  ...
            -   y_true_test(testDays_start, (3*provNum+1):(4*provNum))   )' ;
        trainErr_wHwMwGL_mat1_group3(i,rep) = sum(trainErr_vec_1);
        validErr_wHwMwGL_mat1_group3(i,rep)  = sum(validErr_vec_1);
        testErr_wHwMwGL_mat1_group3(i,rep)  = sum(testErr_vec_1);
        %--------------------------
        trainErr_2  = abs(y_reg_train(:, 1: provNum) ...
            - y_true(1: trainDays_end-startDay+1, 1: provNum)) ...
            ./ y_true(1: trainDays_end-startDay+1, 1: provNum);
        trainErr_vec_2 = (sum(trainErr_2))' / trainDays_end ;
        validErr_2  = abs(y_reg_valid(:, 1: provNum) ...
            - y_true_valid(validDays_start: validDays_end, 1: provNum)) ...
            ./ y_true_valid(validDays_start: validDays_end, 1: provNum);
        validErr_vec_2 = (sum(validErr_2))' / (validDays_end - validDays_start + 1) ;
        testErr_2  = abs(y_reg_test(:, 1: provNum) ...
            - y_true_test(testDays_start: testDays_end, 1: provNum)) ...
            ./ y_true_test(testDays_start: testDays_end, 1: provNum);
        testErr_vec_2 = (sum(testErr_2))' / (testDays_end - testDays_start + 1) ;
        trainErr_wHwMwGL_mat2_group3(i,rep) = sum(trainErr_vec_2);
        validErr_wHwMwGL_mat2_group3(i,rep) = sum(validErr_vec_2);
        testErr_wHwMwGL_mat2_group3(i,rep)  = sum(testErr_vec_2);
        %--------------------------
        trainErr_3  = ( abs(y_reg_train(:, 1: provNum) ...
            - y_true(1: trainDays_end-startDay+1, 1: provNum)) ...
            ./ y_true(1: trainDays_end-startDay+1, 1: provNum) ).^2;
        trainErr_vec_3 = sqrt((sum(trainErr_3 .* y_true(1: trainDays_end-startDay+1, 1: provNum)...
            ./ (ones(trainDays_end,1)*(y_true(trainDays_end-startDay+1, (3*provNum+1):(4*provNum))  ...
            -   y_true(1, (3*provNum+1):(4*provNum)))) ))' ) ;
        validErr_3  = ( abs(y_reg_valid(:, 1: provNum) ...
            - y_true_valid(validDays_start: validDays_end, 1: provNum)) ...
            ./ y_true_valid(validDays_start: validDays_end, 1: provNum) ).^2;
        validErr_vec_3 = sqrt((sum(validErr_3 .* y_true_valid(validDays_start: validDays_end, 1: provNum)...
            ./ (ones(validDays_end-validDays_start+1,1)*(y_true_valid(validDays_end, (3*provNum+1):(4*provNum))  ...
            -   y_true_valid(validDays_start, (3*provNum+1):(4*provNum)))) ))' ) ;
        testErr_3  = ( abs(y_reg_test(:, 1: provNum) ...
            - y_true_test(testDays_start: testDays_end, 1: provNum)) ...
            ./ y_true_test(testDays_start: testDays_end, 1: provNum) ).^2;
        testErr_vec_3 = sqrt((sum(testErr_3 .* y_true_test(testDays_start: testDays_end, 1: provNum)...
            ./ (ones(testDays_end-testDays_start+1,1)*(y_true_test(testDays_end, (3*provNum+1):(4*provNum))  ...
            -   y_true_test(testDays_start, (3*provNum+1):(4*provNum)))) ))' ) ;
        trainErr_wHwMwGL_mat3_group3(i,rep) = sum(trainErr_vec_3);
        validErr_wHwMwGL_mat3_group3(i,rep) = sum(validErr_vec_3);
        testErr_wHwMwGL_mat3_group3(i,rep)  = sum(testErr_vec_3);
        %--------------------------
        trainErr_vec_4 = sqrt( (sum(trainErr_3))' / (trainDays_end-startDay+1) );
        validErr_vec_4 = sqrt( (sum(validErr_3))' / (validDays_end - validDays_start+1) );
        testErr_vec_4 = sqrt( (sum(testErr_3))' / (testDays_end-testDays_start+1) );
        trainErr_wHwMwGL_mat4_group3(i,rep) = sum(trainErr_vec_4);
        validErr_wHwMwGL_mat4_group3(i,rep) = sum(validErr_vec_4);
        testErr_wHwMwGL_mat4_group3(i,rep) = sum(testErr_vec_4);
        
    end
    
    fprintf('\n')
    
    
end


ini.part   = group_cell;


%%

trainDays_end
validDays_end
totDays

fprintf('relative error with weighted average \n\n');

fprintf('\n');
for i = 1: length(muvec)
    fprintf('w/ H, w/ M, w/ GL, mu = 10^%.1f, trainErr: mean %.4f, std %6.4f, testErr: mean %6.4f, std %.4f \n', log10(muvec(i)), mean(trainErr_wHwMwGL_mat1_group3(i,:)/provNum),...
        std(trainErr_wHwMwGL_mat1_group3(i,:)/provNum),  mean(testErr_wHwMwGL_mat1_group3(i,:)/provNum), std(testErr_wHwMwGL_mat1_group3(i,:)/provNum));
end

%%

fprintf('relative error with equal weights \n\n');

fprintf('\n');
for i = 1: length(muvec)    
    fprintf('w/ H, w/ M, w/ GL, mu1 = %6.4f, trainErr: mean %6.4f, std % 6.4f, testErr: mean %6.4f, std % 6.4f \n', muvec(i), mean(trainErr_wHwMwGL_mat2_group3(i,:)/provNum),...
        std(trainErr_wHwMwGL_mat2_group3(i,:)/provNum),  mean(testErr_wHwMwGL_mat2_group3(i,:)/provNum), std(testErr_wHwMwGL_mat2_group3(i,:)/provNum));
    
end
%%
fprintf('\nrelative MSE with weighted average\n\n');

fprintf('\n');
for i = 1: length(muvec)
    fprintf('w/ H, w/ M, w/ GL, mu1 = 10^%.1f, trainErr: mean %6.4f, std % 6.4f, testErr: mean %6.4f, std % 6.4f \n', log10(muvec(i)), mean(trainErr_wHwMwGL_mat3_group3(i,:)/provNum),...
        std(trainErr_wHwMwGL_mat3_group3(i,:)/provNum),  mean(testErr_wHwMwGL_mat3_group3(i,:)/provNum), std(testErr_wHwMwGL_mat3_group3(i,:)/provNum));
    
end
%%

fprintf('\nrelative MSE with simple average\n\n');

fprintf('\n');
for i = 1: length(muvec)    
    fprintf('w/ H, w/ M, w/ GL, mu1 = %6.4f, trainErr: mean %6.4f, std % 6.4f, testErr: mean %6.4f, std % 6.4f \n', muvec(i), mean(trainErr_wHwMwGL_mat4_group3(i,:)/provNum),...
        std(trainErr_wHwMwGL_mat4_group3(i,:)/provNum),  mean(testErr_wHwMwGL_mat4_group3(i,:)/provNum), std(testErr_wHwMwGL_mat4_group3(i,:)/provNum));
    
end


return


%%


fprintf('estimated lambda \n\n');

prov_group1_ind = 4;
prov_group2_ind = 1;
prov_group3_ind = 3;
%-----------------------------------------------

fprintf('w/ H, w/ M, wGL: \n\n');
for i = 1: length(muvec)
    
    fprintf('%2d-th mu=10^%.1f, lambda_k in prov %d in Group 1: mean %6.4f, std % 6.4f\n', ...
        i, log10(muvec(i)),prov_group1_ind,mean(param_wHwMwGL_mat_group3(i,2*provNum+prov_group1_ind,:)), std(param_wHwMwGL_mat_group3(i,2*provNum+prov_group1_ind,:)));
    fprintf('                 lambda_k in prov %d in Group 2: mean %6.4f, std % 6.4f\n', ...
        prov_group2_ind,mean(param_wHwMwGL_mat_group3(i,2*provNum+prov_group2_ind,:)), std(param_wHwMwGL_mat_group3(i,2*provNum+prov_group2_ind,:)));
    fprintf('                 lambda_k in prov %d in Group 3: mean %6.4f, std % 6.4f\n\n', ...
        prov_group3_ind,mean(param_wHwMwGL_mat_group3(i,2*provNum+prov_group3_ind,:)), std(param_wHwMwGL_mat_group3(i,2*provNum+prov_group3_ind,:)));
end


%%



















%% beta = 0.01

beta2 = 0.01

trainErr_wHwMwGL_mat1_beta2 = zeros(length(muvec), numrep);
validErr_wHwMwGL_mat1_beta2 = zeros(length(muvec), numrep);
testErr_wHwMwGL_mat1_beta2  = zeros(length(muvec), numrep);

trainErr_wHwMwGL_mat2_beta2 = zeros(length(muvec), numrep);
validErr_wHwMwGL_mat2_beta2 = zeros(length(muvec), numrep);
testErr_wHwMwGL_mat2_beta2  = zeros(length(muvec), numrep);

trainErr_wHwMwGL_mat3_beta2 = zeros(length(muvec), numrep);
validErr_wHwMwGL_mat3_beta2 = zeros(length(muvec), numrep);
testErr_wHwMwGL_mat3_beta2  = zeros(length(muvec), numrep);

trainErr_wHwMwGL_mat4_beta2 = zeros(length(muvec), numrep);
validErr_wHwMwGL_mat4_beta2 = zeros(length(muvec), numrep);
testErr_wHwMwGL_mat4_beta2  = zeros(length(muvec), numrep);

param_wHwMwGL_mat_beta2 = zeros(length(muvec), 3*provNum, numrep);

% for rep = 1: numrep
for rep = 1: 100
    fprintf('the %d th replica\n', rep);
    
    
    y_true       = reshape(traj_obs_mat(:, :, rep), totDays, 5*provNum);
    y_true_test  = reshape(traj_obs_test_mat(:, :, rep), totDays, 5*provNum);
    y_true_valid = reshape(traj_obs_valid_mat(:, :, rep), totDays, 5*provNum);
    
    %%
    
    %----------------------------------------------------------------------
    % w/ Hetero, w/ Migrat, w/ GL
    fprintf('w/ Hetero, w/ Migrat, w/ GL\n');
    ini.traMat = traMat;
    params_pre = [I_pre  * ones(1, length(provPop)), E_pre * ones(1, length(provPop)), lambda_pre*ones(1,length(provPop))];
    for i = 1: length(muvec)
%     for i = 19
        %for j = 1: length(mu1vec)
        fprintf('the %d-th mu1 = 10^%.1f\n', i, log10(muvec(i)));
        
        %         mu = muvec(i); %mu1vec(j);
        %         f1 = @(x)SEIR_fminunc_test_pred_inf_multi_3(x, y_obs, ini, 'pois', mu, beta, dt, validDays, startDay, sigma );
        
        mu1 = muvec(i)/5000; %mu1vec(j);
        mu0 = mu1*beta2;%mu1/10;
        f1 = @(x)SEIR_fminunc_test_pred_inf_multi_2(x, y_true, ini, 'pois', mu0, mu1, dt, trainDays_end, startDay, sigma );
        
        opts1 = optimoptions('fmincon');
        opts1.MaxIterations =  2000;
        opts1.MaxFunctionEvaluations = 10^6;
        opts1.Display = 'off';%'off';%'iter-detailed';
        opts1.OptimalityTolerance = 1e-10;
        lb = [ones(1,2*provNum), zeros(1,provNum)];
        ub = [ones(1,2*provNum) * 100, ones(1,provNum) * 1];
        tic
        [params1, ~] = fmincon(f1, params_pre, [],[],[],[],lb,ub,[],opts1);
        toc;
        %         if i == 1
        %             [params1, ~] = fmincon(f1, params_pre, [],[],[],[],lb,ub,[],opts1);
        %         else
        %             [params1, ~] = fmincon(f1, reshape(param_wHwMwGL_mat(i-1,:,rep),1,3*provNum), ...
        %                 [],[],[],[],lb,ub,[],opts1);
        %         end
        
        param_wHwMwGL_mat_beta2(i, :, rep) = params1';
        
        R_train = zeros(1, provNum);
        params1_test = params1;
        R_test  = zeros(1, provNum);
        params1_valid = params1;
        R_valid = zeros(1, provNum);
        y_reg_train = SEIR_data_gen_test_pred_inf_determ_2(params1, ini, dt, trainDays_end-startDay+1, 1, R_train);
        y_reg_valid = SEIR_data_gen_test_pred_inf_determ_2(params1_valid, ini, dt, validDays_end-startDay+1, 1, R_valid);
        y_reg_valid = y_reg_valid(validDays_start:validDays_end,:);
        y_reg_test  = SEIR_data_gen_test_pred_inf_determ_2(params1_test, ini, dt, testDays_end-startDay+1, 1, R_test);
        y_reg_test  = y_reg_test(testDays_start:testDays_end,:);
        %-------------------------
        trainErr_1 = abs(y_reg_train(:, 1: provNum) - y_true(1: trainDays_end-startDay+1, 1: provNum));
        trainErr_vec_1 = (sum(trainErr_1))' ./ (y_true(trainDays_end-startDay+1, (3*provNum+1):(4*provNum))  ...
            -   y_true(1, (3*provNum+1):(4*provNum))   )' ;
        validErr_1  = abs(y_reg_valid(:, 1: provNum) ...
            - y_true_valid(validDays_start: validDays_end, 1: provNum));
        validErr_vec_1 = (sum(validErr_1))' ./ (y_true_valid(validDays_end, (3*provNum+1):(4*provNum))  ...
            -   y_true_valid(validDays_start, (3*provNum+1):(4*provNum))   )' ;
        testErr_1  = abs(y_reg_test(:, 1: provNum) ...
            - y_true_test(testDays_start: testDays_end, 1: provNum));
        testErr_vec_1 = (sum(testErr_1))' ./ (y_true_test(testDays_end, (3*provNum+1):(4*provNum))  ...
            -   y_true_test(testDays_start, (3*provNum+1):(4*provNum))   )' ;
        trainErr_wHwMwGL_mat1_beta2(i,rep) = sum(trainErr_vec_1);
        validErr_wHwMwGL_mat1_beta2(i,rep)  = sum(validErr_vec_1);
        testErr_wHwMwGL_mat1_beta2(i,rep)  = sum(testErr_vec_1);
        %--------------------------
        trainErr_2  = abs(y_reg_train(:, 1: provNum) ...
            - y_true(1: trainDays_end-startDay+1, 1: provNum)) ...
            ./ y_true(1: trainDays_end-startDay+1, 1: provNum);
        trainErr_vec_2 = (sum(trainErr_2))' / trainDays_end ;
        validErr_2  = abs(y_reg_valid(:, 1: provNum) ...
            - y_true_valid(validDays_start: validDays_end, 1: provNum)) ...
            ./ y_true_valid(validDays_start: validDays_end, 1: provNum);
        validErr_vec_2 = (sum(validErr_2))' / (validDays_end - validDays_start + 1) ;
        testErr_2  = abs(y_reg_test(:, 1: provNum) ...
            - y_true_test(testDays_start: testDays_end, 1: provNum)) ...
            ./ y_true_test(testDays_start: testDays_end, 1: provNum);
        testErr_vec_2 = (sum(testErr_2))' / (testDays_end - testDays_start + 1) ;
        trainErr_wHwMwGL_mat2_beta2(i,rep) = sum(trainErr_vec_2);
        validErr_wHwMwGL_mat2_beta2(i,rep) = sum(validErr_vec_2);
        testErr_wHwMwGL_mat2_beta2(i,rep)  = sum(testErr_vec_2);
        %--------------------------
        trainErr_3  = ( abs(y_reg_train(:, 1: provNum) ...
            - y_true(1: trainDays_end-startDay+1, 1: provNum)) ...
            ./ y_true(1: trainDays_end-startDay+1, 1: provNum) ).^2;
        trainErr_vec_3 = sqrt((sum(trainErr_3 .* y_true(1: trainDays_end-startDay+1, 1: provNum)...
            ./ (ones(trainDays_end,1)*(y_true(trainDays_end-startDay+1, (3*provNum+1):(4*provNum))  ...
            -   y_true(1, (3*provNum+1):(4*provNum)))) ))' ) ;
        validErr_3  = ( abs(y_reg_valid(:, 1: provNum) ...
            - y_true_valid(validDays_start: validDays_end, 1: provNum)) ...
            ./ y_true_valid(validDays_start: validDays_end, 1: provNum) ).^2;
        validErr_vec_3 = sqrt((sum(validErr_3 .* y_true_valid(validDays_start: validDays_end, 1: provNum)...
            ./ (ones(validDays_end-validDays_start+1,1)*(y_true_valid(validDays_end, (3*provNum+1):(4*provNum))  ...
            -   y_true_valid(validDays_start, (3*provNum+1):(4*provNum)))) ))' ) ;
        testErr_3  = ( abs(y_reg_test(:, 1: provNum) ...
            - y_true_test(testDays_start: testDays_end, 1: provNum)) ...
            ./ y_true_test(testDays_start: testDays_end, 1: provNum) ).^2;
        testErr_vec_3 = sqrt((sum(testErr_3 .* y_true_test(testDays_start: testDays_end, 1: provNum)...
            ./ (ones(testDays_end-testDays_start+1,1)*(y_true_test(testDays_end, (3*provNum+1):(4*provNum))  ...
            -   y_true_test(testDays_start, (3*provNum+1):(4*provNum)))) ))' ) ;
        trainErr_wHwMwGL_mat3_beta2(i,rep) = sum(trainErr_vec_3);
        validErr_wHwMwGL_mat3_beta2(i,rep) = sum(validErr_vec_3);
        testErr_wHwMwGL_mat3_beta2(i,rep)  = sum(testErr_vec_3);
        %--------------------------
        trainErr_vec_4 = sqrt( (sum(trainErr_3))' / (trainDays_end-startDay+1) );
        validErr_vec_4 = sqrt( (sum(validErr_3))' / (validDays_end - validDays_start+1) );
        testErr_vec_4 = sqrt( (sum(testErr_3))' / (testDays_end-testDays_start+1) );
        trainErr_wHwMwGL_mat4_beta2(i,rep) = sum(trainErr_vec_4);
        validErr_wHwMwGL_mat4_beta2(i,rep) = sum(validErr_vec_4);
        testErr_wHwMwGL_mat4_beta2(i,rep) = sum(testErr_vec_4);
        
    end
    
    fprintf('\n')
    
    
end


%%















%% beta = 0.2

beta3 = 0.2

trainErr_wHwMwGL_mat1_beta3 = zeros(length(muvec), numrep);
validErr_wHwMwGL_mat1_beta3 = zeros(length(muvec), numrep);
testErr_wHwMwGL_mat1_beta3  = zeros(length(muvec), numrep);

trainErr_wHwMwGL_mat2_beta3 = zeros(length(muvec), numrep);
validErr_wHwMwGL_mat2_beta3 = zeros(length(muvec), numrep);
testErr_wHwMwGL_mat2_beta3  = zeros(length(muvec), numrep);

trainErr_wHwMwGL_mat3_beta3 = zeros(length(muvec), numrep);
validErr_wHwMwGL_mat3_beta3 = zeros(length(muvec), numrep);
testErr_wHwMwGL_mat3_beta3  = zeros(length(muvec), numrep);

trainErr_wHwMwGL_mat4_beta3 = zeros(length(muvec), numrep);
validErr_wHwMwGL_mat4_beta3 = zeros(length(muvec), numrep);
testErr_wHwMwGL_mat4_beta3  = zeros(length(muvec), numrep);

param_wHwMwGL_mat_beta3 = zeros(length(muvec), 3*provNum, numrep);

% for rep = 1: numrep
for rep = 1: 100
    fprintf('the %d th replica\n', rep);
    
    
    y_true       = reshape(traj_obs_mat(:, :, rep), totDays, 5*provNum);
    y_true_test  = reshape(traj_obs_test_mat(:, :, rep), totDays, 5*provNum);
    y_true_valid = reshape(traj_obs_valid_mat(:, :, rep), totDays, 5*provNum);
    
    %%
    
    %----------------------------------------------------------------------
    % w/ Hetero, w/ Migrat, w/ GL
    fprintf('w/ Hetero, w/ Migrat, w/ GL\n');
    ini.traMat = traMat;
    params_pre = [I_pre  * ones(1, length(provPop)), E_pre * ones(1, length(provPop)), lambda_pre*ones(1,length(provPop))];
    for i = 1: length(muvec)
%     for i = 19
        %for j = 1: length(mu1vec)
        fprintf('the %d-th mu1 = 10^%.1f\n', i, log10(muvec(i)));
        
        %         mu = muvec(i); %mu1vec(j);
        %         f1 = @(x)SEIR_fminunc_test_pred_inf_multi_3(x, y_obs, ini, 'pois', mu, beta, dt, validDays, startDay, sigma );
        
        mu1 = muvec(i)/5000; %mu1vec(j);
        mu0 = mu1*beta3;%mu1/10;
        f1 = @(x)SEIR_fminunc_test_pred_inf_multi_2(x, y_true, ini, 'pois', mu0, mu1, dt, trainDays_end, startDay, sigma );
        
        opts1 = optimoptions('fmincon');
        opts1.MaxIterations =  2000;
        opts1.MaxFunctionEvaluations = 10^6;
        opts1.Display = 'off';%'off';%'iter-detailed';
        opts1.OptimalityTolerance = 1e-10;
        lb = [ones(1,2*provNum), zeros(1,provNum)];
        ub = [ones(1,2*provNum) * 100, ones(1,provNum) * 1];
        tic
        [params1, ~] = fmincon(f1, params_pre, [],[],[],[],lb,ub,[],opts1);
        toc;
        %         if i == 1
        %             [params1, ~] = fmincon(f1, params_pre, [],[],[],[],lb,ub,[],opts1);
        %         else
        %             [params1, ~] = fmincon(f1, reshape(param_wHwMwGL_mat(i-1,:,rep),1,3*provNum), ...
        %                 [],[],[],[],lb,ub,[],opts1);
        %         end
        
        param_wHwMwGL_mat_beta3(i, :, rep) = params1';
        
        R_train = zeros(1, provNum);
        params1_test = params1;
        R_test  = zeros(1, provNum);
        params1_valid = params1;
        R_valid = zeros(1, provNum);
        y_reg_train = SEIR_data_gen_test_pred_inf_determ_2(params1, ini, dt, trainDays_end-startDay+1, 1, R_train);
        y_reg_valid = SEIR_data_gen_test_pred_inf_determ_2(params1_valid, ini, dt, validDays_end-startDay+1, 1, R_valid);
        y_reg_valid = y_reg_valid(validDays_start:validDays_end,:);
        y_reg_test  = SEIR_data_gen_test_pred_inf_determ_2(params1_test, ini, dt, testDays_end-startDay+1, 1, R_test);
        y_reg_test  = y_reg_test(testDays_start:testDays_end,:);
        %-------------------------
        trainErr_1 = abs(y_reg_train(:, 1: provNum) - y_true(1: trainDays_end-startDay+1, 1: provNum));
        trainErr_vec_1 = (sum(trainErr_1))' ./ (y_true(trainDays_end-startDay+1, (3*provNum+1):(4*provNum))  ...
            -   y_true(1, (3*provNum+1):(4*provNum))   )' ;
        validErr_1  = abs(y_reg_valid(:, 1: provNum) ...
            - y_true_valid(validDays_start: validDays_end, 1: provNum));
        validErr_vec_1 = (sum(validErr_1))' ./ (y_true_valid(validDays_end, (3*provNum+1):(4*provNum))  ...
            -   y_true_valid(validDays_start, (3*provNum+1):(4*provNum))   )' ;
        testErr_1  = abs(y_reg_test(:, 1: provNum) ...
            - y_true_test(testDays_start: testDays_end, 1: provNum));
        testErr_vec_1 = (sum(testErr_1))' ./ (y_true_test(testDays_end, (3*provNum+1):(4*provNum))  ...
            -   y_true_test(testDays_start, (3*provNum+1):(4*provNum))   )' ;
        trainErr_wHwMwGL_mat1_beta3(i,rep) = sum(trainErr_vec_1);
        validErr_wHwMwGL_mat1_beta3(i,rep)  = sum(validErr_vec_1);
        testErr_wHwMwGL_mat1_beta3(i,rep)  = sum(testErr_vec_1);
        %--------------------------
        trainErr_2  = abs(y_reg_train(:, 1: provNum) ...
            - y_true(1: trainDays_end-startDay+1, 1: provNum)) ...
            ./ y_true(1: trainDays_end-startDay+1, 1: provNum);
        trainErr_vec_2 = (sum(trainErr_2))' / trainDays_end ;
        validErr_2  = abs(y_reg_valid(:, 1: provNum) ...
            - y_true_valid(validDays_start: validDays_end, 1: provNum)) ...
            ./ y_true_valid(validDays_start: validDays_end, 1: provNum);
        validErr_vec_2 = (sum(validErr_2))' / (validDays_end - validDays_start + 1) ;
        testErr_2  = abs(y_reg_test(:, 1: provNum) ...
            - y_true_test(testDays_start: testDays_end, 1: provNum)) ...
            ./ y_true_test(testDays_start: testDays_end, 1: provNum);
        testErr_vec_2 = (sum(testErr_2))' / (testDays_end - testDays_start + 1) ;
        trainErr_wHwMwGL_mat2_beta3(i,rep) = sum(trainErr_vec_2);
        validErr_wHwMwGL_mat2_beta3(i,rep) = sum(validErr_vec_2);
        testErr_wHwMwGL_mat2_beta3(i,rep)  = sum(testErr_vec_2);
        %--------------------------
        trainErr_3  = ( abs(y_reg_train(:, 1: provNum) ...
            - y_true(1: trainDays_end-startDay+1, 1: provNum)) ...
            ./ y_true(1: trainDays_end-startDay+1, 1: provNum) ).^2;
        trainErr_vec_3 = sqrt((sum(trainErr_3 .* y_true(1: trainDays_end-startDay+1, 1: provNum)...
            ./ (ones(trainDays_end,1)*(y_true(trainDays_end-startDay+1, (3*provNum+1):(4*provNum))  ...
            -   y_true(1, (3*provNum+1):(4*provNum)))) ))' ) ;
        validErr_3  = ( abs(y_reg_valid(:, 1: provNum) ...
            - y_true_valid(validDays_start: validDays_end, 1: provNum)) ...
            ./ y_true_valid(validDays_start: validDays_end, 1: provNum) ).^2;
        validErr_vec_3 = sqrt((sum(validErr_3 .* y_true_valid(validDays_start: validDays_end, 1: provNum)...
            ./ (ones(validDays_end-validDays_start+1,1)*(y_true_valid(validDays_end, (3*provNum+1):(4*provNum))  ...
            -   y_true_valid(validDays_start, (3*provNum+1):(4*provNum)))) ))' ) ;
        testErr_3  = ( abs(y_reg_test(:, 1: provNum) ...
            - y_true_test(testDays_start: testDays_end, 1: provNum)) ...
            ./ y_true_test(testDays_start: testDays_end, 1: provNum) ).^2;
        testErr_vec_3 = sqrt((sum(testErr_3 .* y_true_test(testDays_start: testDays_end, 1: provNum)...
            ./ (ones(testDays_end-testDays_start+1,1)*(y_true_test(testDays_end, (3*provNum+1):(4*provNum))  ...
            -   y_true_test(testDays_start, (3*provNum+1):(4*provNum)))) ))' ) ;
        trainErr_wHwMwGL_mat3_beta3(i,rep) = sum(trainErr_vec_3);
        validErr_wHwMwGL_mat3_beta3(i,rep) = sum(validErr_vec_3);
        testErr_wHwMwGL_mat3_beta3(i,rep)  = sum(testErr_vec_3);
        %--------------------------
        trainErr_vec_4 = sqrt( (sum(trainErr_3))' / (trainDays_end-startDay+1) );
        validErr_vec_4 = sqrt( (sum(validErr_3))' / (validDays_end - validDays_start+1) );
        testErr_vec_4 = sqrt( (sum(testErr_3))' / (testDays_end-testDays_start+1) );
        trainErr_wHwMwGL_mat4_beta3(i,rep) = sum(trainErr_vec_4);
        validErr_wHwMwGL_mat4_beta3(i,rep) = sum(validErr_vec_4);
        testErr_wHwMwGL_mat4_beta3(i,rep) = sum(testErr_vec_4);
        
    end
    
    fprintf('\n')
    
    
end















%%

map_ind_to_ind = containers.Map([group_cell{1}, group_cell{2}, group_cell{3}], (1:30));

ind_to_ind = [group_cell{1}, group_cell{2}, group_cell{3}];

%%



rep = 6;


FontSize = 36;
FontSize_legend = 36;
% rep = 100;


% group_cell

prov_group1_ind = 4;
prov_group2_ind = 1;
prov_group3_ind = 3;


muind = 19;

param_woHwoM_mat(2*provNum+1,rep)
param_woHwM_mat(2*provNum+1,rep)
[param_wHwoM_mat(2*provNum+prov_group1_ind,rep),param_wHwoM_mat(2*provNum+prov_group2_ind,rep),param_wHwoM_mat(2*provNum+prov_group3_ind,rep)]
[param_wHwM_mat(2*provNum+prov_group1_ind,rep),param_wHwM_mat(2*provNum+prov_group2_ind,rep),param_wHwM_mat(2*provNum+prov_group3_ind,rep)]
[param_wHwMwGL_mat(muind,2*provNum+prov_group1_ind,rep),param_wHwMwGL_mat(muind,2*provNum+prov_group2_ind,rep),param_wHwMwGL_mat(muind,2*provNum+prov_group3_ind,rep)]


round([param_woHwoM_mat(prov_group1_ind,rep),param_woHwoM_mat(prov_group2_ind,rep),param_woHwoM_mat(prov_group3_ind,rep)])
round([param_wHwM_mat(prov_group1_ind,rep),param_wHwM_mat(prov_group2_ind,rep),param_wHwM_mat(prov_group3_ind,rep)])


% [reshape(param_wHwMwGL_mat(muind,[2*provNum+prov_group1_ind,2*provNum+prov_group2_ind,2*provNum+prov_group3_ind],:),3,100)',...
%     param_wHwM_mat(2*provNum+prov_group1_ind,:)',param_wHwM_mat(2*provNum+prov_group2_ind,:)',param_wHwM_mat(2*provNum+prov_group3_ind,:)']

y_true_test = reshape(traj_obs_test_mat(:,:,rep),totDays,5*provNum);
params1 = param_wHwMwGL_mat(muind, :, rep);
% params1 = param_true;
R_train = zeros(1, provNum);
params1_test = params1;
% params1_test(11) = 0.4;
% params1_test(1:4) = 10;
% params1_test = params1;
params1_test = param;
% params1_test(1: provNum) = y_true_test(testDays_start-1, (3*provNum+1):(4*provNum)); % change I
% params1_test((provNum+1):(2*provNum)) = y_true_test(testDays_start-1, (2*provNum+1):(3*provNum)); % change E
R_test = zeros(1, provNum);%y_true_test(testDays_start-1, (4*provNum+1):(5*provNum)); % change R
params1_valid = params1;
% params1_valid(1: provNum) = y_true_valid(validDays_start-1, (3*provNum+1):(4*provNum)); % change I
% params1_valid((provNum+1):(2*provNum)) = y_true_valid(validDays_start-1, (2*provNum+1):(3*provNum)); % change E
R_valid = zeros(1, provNum);%y_true_test(validDays_start-1, (4*provNum+1):(5*provNum)); % change R

traj_train_wHwMwGL_temp = SEIR_data_gen_test_pred_inf_determ_2(params1, ini, dt, trainDays_end-startDay+1, 1, R_train);
traj_test_wHwMwGL_temp  = SEIR_data_gen_test_pred_inf_determ_2(params1_test, ini, dt, testDays_end-startDay+1, 1, R_test);
traj_test_wHwMwGL_temp  = traj_test_wHwMwGL_temp(testDays_start:testDays_end,:);

%------------------------------------------------------------

% params1_group2 = param_wHwMwGL_mat_group2(muind, :, rep);
% traj_train_wHwMwGL_temp_group2 = SEIR_data_gen_test_pred_inf_determ_2(params1_group2, ini, dt, trainDays_end-startDay+1, 1, R_train);
% traj_test_wHwMwGL_temp_group2  = SEIR_data_gen_test_pred_inf_determ_2(params1_group2, ini, dt, testDays_end-startDay+1, 1, R_test);
% traj_test_wHwMwGL_temp_group2  = traj_test_wHwMwGL_temp_group2(testDays_start:testDays_end,:);

%------------------------------------------------------------

% params1_group3 = param_wHwMwGL_mat_group3(muind, :, rep);
% traj_train_wHwMwGL_temp_group3 = SEIR_data_gen_test_pred_inf_determ_2(params1_group3, ini, dt, trainDays_end-startDay+1, 1, R_train);
% traj_test_wHwMwGL_temp_group3  = SEIR_data_gen_test_pred_inf_determ_2(params1_group3, ini, dt, testDays_end-startDay+1, 1, R_test);
% traj_test_wHwMwGL_temp_group3  = traj_test_wHwMwGL_temp_group3(testDays_start:testDays_end,:);



% traj_valid_wHwMwGL_temp = SEIR_data_gen_test_pred_inf_determ_2(params1_valid, ini, dt, validDays_end-startDay+1, 1, R_valid);
% traj_valid_wHwMwGL_temp = traj_valid_wHwMwGL_temp(validDays_start:validDays_end,:);


% param_true = param;
traj_ground_truth_determ = SEIR_data_gen_test_pred_inf_determ_2(param_true, ini, dt, testDays_end-startDay+1, 1, R_train);


% ini.traMat = zeros(provNum);
% tempxx2 = SEIR_data_gen_test_pred_inf_determ_wo_h_2(param_woHwoM_mat(:, rep)', ini, dt, testDays_end-startDay+1, 1, R_train);
% ini.traMat = traMat;


figure(3), clf;
% subplot(2,2,1);
subplot(1,3,1);
idx = prov_group1_ind;
temp1 = [traj_train_woHwoM_mat(:, idx, rep)', traj_test_woHwoM_mat(:, idx, rep)'];
temp2 = [traj_train_woHwM_mat(:, idx, rep)', traj_test_woHwM_mat(:, idx, rep)'];
temp3 = [traj_train_wHwoM_mat(:, idx, rep)', traj_test_wHwoM_mat(:, idx, rep)'];
temp4 = [traj_train_wHwM_mat(:, idx, rep)', traj_test_wHwM_mat(:, idx, rep)'];
temp5 = [traj_train_wHwMwGL_temp(:, idx)', traj_test_wHwMwGL_temp(:, idx)'];
% temp5_group2 = [traj_train_wHwMwGL_temp_group2(:, idx)', traj_test_wHwMwGL_temp_group2(:, idx)'];
% temp5_group3 = [traj_train_wHwMwGL_temp_group3(:, idx)', traj_test_wHwMwGL_temp_group3(:, idx)'];
temp6 = [traj_obs_mat(startDay:trainDays_end, idx, rep)', traj_obs_test_mat(testDays_start:testDays_end, idx, rep)'];
hold on;
plot(1:testDays_end, temp1, 'x-', 'LineWidth', 2, 'MarkerSize', 8);
plot(1:testDays_end, temp2, 'd-', 'LineWidth', 2, 'color', '#cc99ff', 'MarkerSize', 8);
plot(1:testDays_end, temp3, 'p-', 'LineWidth', 2, 'MarkerSize', 8);
plot(1:testDays_end, temp4, 'h-', 'LineWidth', 2, 'MarkerSize', 8);
plot(1:testDays_end, temp5, 's-', 'LineWidth', 2, 'MarkerSize', 8);
% plot(1:testDays_end, temp5_group2, 'x-', 'LineWidth', 2, 'color', '#99ff33');
% plot(1:testDays_end, temp5_group3, 'x-', 'LineWidth', 2, 'color', '#336600');
plot(1:testDays_end, temp6, 'o-', 'LineWidth', 2, 'color', '#4DBEEE');
% plot(1:testDays_end, traj_ground_truth_determ(:, idx)', 'o-', 'LineWidth', 2, 'color', 'black');
set(gca, 'FontSize', FontSize);
title(sprintf('Province %d in Group 1', map_ind_to_ind(prov_group1_ind)), 'FontSize', 36, 'interpreter', 'latex');
xlabel('day', 'interpreter', 'latex'); ylabel('cases', 'interpreter', 'latex');
xline(validDays_end, 'LineWidth', 2);
xline(trainDays_end, 'LineWidth', 2, 'color', '#A2142F');
grid on;
% yticks([0 1e4 2e4 3e4 4e4 5e4])
% yticklabels({'0', '10^4','2\times10^4','3\times10^4','4\times10^4','5\times10^4'})
% yticklabels({'0', '10,000','20,000','30,000','40,000','50,000'})
% legend('woHwoM', 'woHwM', 'wHwoM', 'wHwM', 'wHwMwGL','True','Location','best');
% legend('Model 1', 'Model 2', 'Model 3', 'Model 4', 'Model 5 (P)', 'Model 5 (P'')', 'Model 5 (P'''')','True','Location','best', 'FontSize', 24);
legend('Model 1', 'Model 2', 'Model 3', 'Model 4', 'Model 5','True','Location','best', 'FontSize', FontSize_legend, 'interpreter', 'latex');


% subplot(2,2,2);
subplot(1,3,2);
idx = prov_group2_ind;
temp1 = [traj_train_woHwoM_mat(:, idx, rep)', traj_test_woHwoM_mat(:, idx, rep)'];
temp2 = [traj_train_woHwM_mat(:, idx, rep)', traj_test_woHwM_mat(:, idx, rep)'];
temp3 = [traj_train_wHwoM_mat(:, idx, rep)', traj_test_wHwoM_mat(:, idx, rep)'];
temp4 = [traj_train_wHwM_mat(:, idx, rep)', traj_test_wHwM_mat(:, idx, rep)'];
temp5 = [traj_train_wHwMwGL_temp(:, idx)', traj_test_wHwMwGL_temp(:, idx)'];
% temp5_group2 = [traj_train_wHwMwGL_temp_group2(:, idx)', traj_test_wHwMwGL_temp_group2(:, idx)'];
% temp5_group3 = [traj_train_wHwMwGL_temp_group3(:, idx)', traj_test_wHwMwGL_temp_group3(:, idx)'];
temp6 = [traj_obs_mat(startDay:trainDays_end, idx, rep)', traj_obs_test_mat(testDays_start:testDays_end, idx, rep)'];
hold on;
plot(1:testDays_end, temp1, 'x-', 'LineWidth', 2, 'MarkerSize', 8);
plot(1:testDays_end, temp2, 'd-', 'LineWidth', 2, 'color', '#cc99ff', 'MarkerSize', 8);
plot(1:testDays_end, temp3, 'p-', 'LineWidth', 2, 'MarkerSize', 8);
plot(1:testDays_end, temp4, 'h-', 'LineWidth', 2, 'MarkerSize', 8);
plot(1:testDays_end, temp5, 's-', 'LineWidth', 2, 'MarkerSize', 8);
% plot(1:testDays_end, temp5_group2, 'x-', 'LineWidth', 2, 'color', '#99ff33');
% plot(1:testDays_end, temp5_group3, 'x-', 'LineWidth', 2, 'color', '#336600');
plot(1:testDays_end, temp6, 'o-', 'LineWidth', 2, 'color', '#4DBEEE');
% plot(1:testDays_end, traj_ground_truth_determ(:, idx)', 'o-', 'LineWidth', 2, 'color', 'black');
set(gca, 'FontSize', FontSize);
title(sprintf('Province %d in Group 2', map_ind_to_ind(prov_group2_ind)), 'FontSize', 36, 'interpreter', 'latex');
xlabel('day', 'interpreter', 'latex'); ylabel('cases', 'interpreter', 'latex');
xline(validDays_end, 'LineWidth', 2);
xline(trainDays_end, 'LineWidth', 2, 'color', '#A2142F');
grid on;
% yticks([0 1e4 2e4 3e4 4e4 5e4])
% yticklabels({'0', '10^4','2\times10^4','3\times10^4','4\times10^4','5\times10^4'})
% yticklabels({'0', '10,000','20,000','30,000','40,000','50,000'})
% legend('woHwoM', 'woHwM', 'wHwoM', 'wHwM', 'wHwMwGL','True','Location','best');
% legend('Model 1', 'Model 2', 'Model 3', 'Model 4', 'Model 5 (P)', 'Model 5 (P'')', 'Model 5 (P'''')','True','Location','best', 'FontSize', 24);
legend('Model 1', 'Model 2', 'Model 3', 'Model 4', 'Model 5','True','Location','best', 'FontSize', FontSize_legend, 'interpreter', 'latex');


% subplot(2,2,2);
subplot(1,3,3);
idx = prov_group3_ind;
temp1 = [traj_train_woHwoM_mat(:, idx, rep)', traj_test_woHwoM_mat(:, idx, rep)'];
temp2 = [traj_train_woHwM_mat(:, idx, rep)', traj_test_woHwM_mat(:, idx, rep)'];
temp3 = [traj_train_wHwoM_mat(:, idx, rep)', traj_test_wHwoM_mat(:, idx, rep)'];
temp4 = [traj_train_wHwM_mat(:, idx, rep)', traj_test_wHwM_mat(:, idx, rep)'];
temp5 = [traj_train_wHwMwGL_temp(:, idx)', traj_test_wHwMwGL_temp(:, idx)'];
% temp5_group2 = [traj_train_wHwMwGL_temp_group2(:, idx)', traj_test_wHwMwGL_temp_group2(:, idx)'];
% temp5_group3 = [traj_train_wHwMwGL_temp_group3(:, idx)', traj_test_wHwMwGL_temp_group3(:, idx)'];
temp6 = [traj_obs_mat(startDay:trainDays_end, idx, rep)', traj_obs_test_mat(testDays_start:testDays_end, idx, rep)'];
hold on;
plot(1:testDays_end, temp1, 'x-', 'LineWidth', 2, 'MarkerSize', 8);
plot(1:testDays_end, temp2, 'd-', 'LineWidth', 2, 'color', '#cc99ff', 'MarkerSize', 8);
plot(1:testDays_end, temp3, 'p-', 'LineWidth', 2, 'MarkerSize', 8);
plot(1:testDays_end, temp4, 'h-', 'LineWidth', 2, 'MarkerSize', 8);
plot(1:testDays_end, temp5, 's-', 'LineWidth', 2, 'MarkerSize', 8);
% plot(1:testDays_end, temp5_group2, 'x-', 'LineWidth', 2, 'color', '#99ff33');
% plot(1:testDays_end, temp5_group3, 'x-', 'LineWidth', 2, 'color', '#336600');
plot(1:testDays_end, temp6, 'o-', 'LineWidth', 2, 'color', '#4DBEEE');
set(gca, 'FontSize', FontSize);
title(sprintf('Province %d in Group 3', map_ind_to_ind(prov_group3_ind)), 'FontSize', 36, 'interpreter', 'latex');
xlabel('day', 'interpreter', 'latex'); ylabel('cases', 'interpreter', 'latex');
xline(validDays_end, 'LineWidth', 2);
xline(trainDays_end, 'LineWidth', 2, 'color', '#A2142F');
grid on;
yticks([0 1e4 2e4 3e4])
yticklabels({'0', '10000','20000','30000'})
% yticklabels({'0', '10,000','20,000','30,000','40,000','50,000'})
% legend('woHwoM', 'woHwM', 'wHwoM', 'wHwM', 'wHwMwGL','True','Location','best');
% legend('Model 1', 'Model 2', 'Model 3', 'Model 4', 'Model 5 (P)', 'Model 5 (P'')', 'Model 5 (P'''')','True','Location','best', 'FontSize', 24);
legend('Model 1', 'Model 2', 'Model 3', 'Model 4', 'Model 5','True','Location','best', 'FontSize', FontSize_legend, 'interpreter', 'latex');


% subplot(2,2,4);
% idx = 4;
% temp1 = [traj_train_woHwoM_mat(:, idx, rep)', traj_test_woHwoM_mat(:, idx, rep)'];
% temp2 = [traj_train_woHwM_mat(:, idx, rep)', traj_test_woHwM_mat(:, idx, rep)'];
% temp3 = [traj_train_wHwoM_mat(:, idx, rep)', traj_test_wHwoM_mat(:, idx, rep)'];
% temp4 = [traj_train_wHwM_mat(:, idx, rep)', traj_test_wHwM_mat(:, idx, rep)'];
% temp5 = [traj_train_wHwMwGL_temp(:, idx)', traj_test_wHwMwGL_temp(:, idx)'];
% temp6 = [traj_obs_mat(startDay:trainDays_end, idx, rep)', traj_obs_test_mat(testDays_start:testDays_end, idx, rep)'];
% hold on;
% plot(1:testDays_end, temp1, 'x-', 'LineWidth', 2);
% plot(1:testDays_end, temp2, 'x-', 'LineWidth', 2);
% plot(1:testDays_end, temp3, 'x-', 'LineWidth', 2);
% plot(1:testDays_end, temp4, 'x-', 'LineWidth', 2);
% plot(1:testDays_end, temp5, 'x-', 'LineWidth', 2);
% plot(1:testDays_end, temp6, 'o-', 'LineWidth', 2);
% set(gca, 'FontSize', 24);
% title('Province 4', 'FontSize', 30);
% xlabel('day'); ylabel('cases');
% xline(validDays_end, 'LineWidth', 2);
% xline(trainDays_end, 'LineWidth', 2, 'color', '#A2142F');
% grid on;
% % yticks([0 1e4 2e4])
% % yticklabels({'0', '10^4','2\times10^4'})
% % yticklabels({'0', '10,000','20,000'});
% % legend('woHwoM', 'woHwM', 'wHwoM', 'wHwM', 'wHwMwGL','True','Location','best');
% legend('Model 1', 'Model 2', 'Model 3', 'Model 4', 'Model 5','True','Location','best', 'FontSize', 24);

sgtitle({"True trajectory \& Fitted trajectories", ''}, 'FontSize', 42, 'interpreter', 'latex');

%

%rep = 83;
figure(4), clf;
% subplot(2,2,2);
subplot(1,3,1);
idx = prov_group1_ind;
temp1 = [abs(traj_train_woHwoM_mat(:, idx, rep)-traj_obs_mat(startDay:trainDays_end, idx, rep))', abs(traj_test_woHwoM_mat(:, idx, rep)-traj_obs_test_mat(testDays_start:testDays_end, idx, rep))'];
temp2 = [abs(traj_train_woHwM_mat(:, idx, rep)-traj_obs_mat(startDay:trainDays_end, idx, rep))', abs(traj_test_woHwM_mat(:, idx, rep)-traj_obs_test_mat(testDays_start:testDays_end, idx, rep))'];
temp3 = [abs(traj_train_wHwoM_mat(:, idx, rep)-traj_obs_mat(startDay:trainDays_end, idx, rep))', abs(traj_test_wHwoM_mat(:, idx, rep)-traj_obs_test_mat(testDays_start:testDays_end, idx, rep))'];
temp4 = [abs(traj_train_wHwM_mat(:, idx, rep)-traj_obs_mat(startDay:trainDays_end, idx, rep))', abs(traj_test_wHwM_mat(:, idx, rep)-traj_obs_test_mat(testDays_start:testDays_end, idx, rep))'];
% temp5 = [abs(traj_train_wHwMwGL_temp(:, idx)-traj_obs_mat(startDay:trainDays_end, idx, rep))', abs(traj_test_wHwMwGL_temp(:, idx)-traj_obs_test_mat(testDays_start:testDays_end, idx, rep))'];
temp5 = [abs(traj_train_wHwMwGL_temp(:, idx)-traj_obs_mat(startDay:trainDays_end, idx, rep))', abs(traj_test_wHwMwGL_temp(:, idx)-traj_obs_test_mat(testDays_start:testDays_end, idx, rep))'];
% temp5_group2 = [abs(traj_train_wHwMwGL_temp_group2(:, idx)-traj_obs_mat(startDay:trainDays_end, idx, rep))', abs(traj_test_wHwMwGL_temp_group2(:, idx)-traj_obs_test_mat(testDays_start:testDays_end, idx, rep))'];
% temp5_group3 = [abs(traj_train_wHwMwGL_temp_group3(:, idx)-traj_obs_mat(startDay:trainDays_end, idx, rep))', abs(traj_test_wHwMwGL_temp_group3(:, idx)-traj_obs_test_mat(testDays_start:testDays_end, idx, rep))'];
hold on;
plot(1:testDays_end, temp1, 'x-', 'LineWidth', 2, 'MarkerSize', 8);
plot(1:testDays_end, temp2, 'd-', 'LineWidth', 2, 'color', '#cc99ff', 'MarkerSize', 8);
plot(1:testDays_end, temp3, 'p-', 'LineWidth', 2, 'MarkerSize', 8);
plot(1:testDays_end, temp4, 'h-', 'LineWidth', 2, 'MarkerSize', 8);
plot(1:testDays_end, temp5, 's-', 'LineWidth', 2, 'MarkerSize', 8);
% plot(1:testDays_end, temp5_group2, 'x-', 'LineWidth', 2, 'color', '#99ff33');
% plot(1:testDays_end, temp5_group3, 'x-', 'LineWidth', 2, 'color', '#336600');

set(gca, 'FontSize', FontSize);
title(sprintf('Province %d in Group 1', map_ind_to_ind(prov_group1_ind)), 'FontSize', 36, 'interpreter', 'latex');
xlabel('day', 'interpreter', 'latex');
ylabel('Absolute Errors', 'interpreter', 'latex');
xline(validDays_end, 'LineWidth', 2);
xline(trainDays_end, 'LineWidth', 2, 'color', '#A2142F');
grid on;
% yticks([0 5e3 1e4 1.5e4])
% yticklabels({'0', '5,000','10,000','15,000'})
% legend('woHwoM', 'woHwM', 'wHwoM', 'wHwM', 'wHwMwGL','Location','best');
% legend('Model 1', 'Model 2', 'Model 3', 'Model 4', 'Model 5 (P)', 'Model 5 (P'')', 'Model 5 (P'''')','Location','best', 'FontSize', 24);
legend('Model 1', 'Model 2', 'Model 3', 'Model 4', 'Model 5','Location','best', 'FontSize', FontSize_legend, 'interpreter', 'latex');


% subplot(2,2,2);
subplot(1,3,2);
idx = prov_group2_ind;
temp1 = [abs(traj_train_woHwoM_mat(:, idx, rep)-traj_obs_mat(startDay:trainDays_end, idx, rep))', abs(traj_test_woHwoM_mat(:, idx, rep)-traj_obs_test_mat(testDays_start:testDays_end, idx, rep))'];
temp2 = [abs(traj_train_woHwM_mat(:, idx, rep)-traj_obs_mat(startDay:trainDays_end, idx, rep))', abs(traj_test_woHwM_mat(:, idx, rep)-traj_obs_test_mat(testDays_start:testDays_end, idx, rep))'];
temp3 = [abs(traj_train_wHwoM_mat(:, idx, rep)-traj_obs_mat(startDay:trainDays_end, idx, rep))', abs(traj_test_wHwoM_mat(:, idx, rep)-traj_obs_test_mat(testDays_start:testDays_end, idx, rep))'];
temp4 = [abs(traj_train_wHwM_mat(:, idx, rep)-traj_obs_mat(startDay:trainDays_end, idx, rep))', abs(traj_test_wHwM_mat(:, idx, rep)-traj_obs_test_mat(testDays_start:testDays_end, idx, rep))'];
% temp5 = [abs(traj_train_wHwMwGL_temp(:, idx)-traj_obs_mat(startDay:trainDays_end, idx, rep))', abs(traj_test_wHwMwGL_temp(:, idx)-traj_obs_test_mat(testDays_start:testDays_end, idx, rep))'];
temp5 = [abs(traj_train_wHwMwGL_temp(:, idx)-traj_obs_mat(startDay:trainDays_end, idx, rep))', abs(traj_test_wHwMwGL_temp(:, idx)-traj_obs_test_mat(testDays_start:testDays_end, idx, rep))'];
% temp5_group2 = [abs(traj_train_wHwMwGL_temp_group2(:, idx)-traj_obs_mat(startDay:trainDays_end, idx, rep))', abs(traj_test_wHwMwGL_temp_group2(:, idx)-traj_obs_test_mat(testDays_start:testDays_end, idx, rep))'];
% temp5_group3 = [abs(traj_train_wHwMwGL_temp_group3(:, idx)-traj_obs_mat(startDay:trainDays_end, idx, rep))', abs(traj_test_wHwMwGL_temp_group3(:, idx)-traj_obs_test_mat(testDays_start:testDays_end, idx, rep))'];
hold on;
plot(1:testDays_end, temp1, 'x-', 'LineWidth', 2, 'MarkerSize', 8);
plot(1:testDays_end, temp2, 'd-', 'LineWidth', 2, 'color', '#cc99ff', 'MarkerSize', 8);
plot(1:testDays_end, temp3, 'p-', 'LineWidth', 2, 'MarkerSize', 8);
plot(1:testDays_end, temp4, 'h-', 'LineWidth', 2, 'MarkerSize', 8);
plot(1:testDays_end, temp5, 's-', 'LineWidth', 2, 'MarkerSize', 8);
% plot(1:testDays_end, temp5_group2, 'x-', 'LineWidth', 2, 'color', '#99ff33');
% plot(1:testDays_end, temp5_group3, 'x-', 'LineWidth', 2, 'color', '#336600');

set(gca, 'FontSize', FontSize);
title(sprintf('Province %d in Group 2', map_ind_to_ind(prov_group2_ind)), 'FontSize', 36, 'interpreter', 'latex');
xlabel('day', 'interpreter', 'latex');
ylabel('Absolute Errors', 'interpreter', 'latex');
xline(validDays_end, 'LineWidth', 2);
xline(trainDays_end, 'LineWidth', 2, 'color', '#A2142F');
grid on;
% yticks([0 5e3 1e4 1.5e4])
% yticklabels({'0', '5,000','10,000','15,000'})
% legend('woHwoM', 'woHwM', 'wHwoM', 'wHwM', 'wHwMwGL','Location','best');
% legend('Model 1', 'Model 2', 'Model 3', 'Model 4', 'Model 5 (P)', 'Model 5 (P'')', 'Model 5 (P'''')','Location','best', 'FontSize', 24);
legend('Model 1', 'Model 2', 'Model 3', 'Model 4', 'Model 5','Location','best', 'FontSize', FontSize_legend, 'interpreter', 'latex');


% subplot(2,2,2);
subplot(1,3,3);
idx = prov_group3_ind;
temp1 = [abs(traj_train_woHwoM_mat(:, idx, rep)-traj_obs_mat(startDay:trainDays_end, idx, rep))', abs(traj_test_woHwoM_mat(:, idx, rep)-traj_obs_test_mat(testDays_start:testDays_end, idx, rep))'];
temp2 = [abs(traj_train_woHwM_mat(:, idx, rep)-traj_obs_mat(startDay:trainDays_end, idx, rep))', abs(traj_test_woHwM_mat(:, idx, rep)-traj_obs_test_mat(testDays_start:testDays_end, idx, rep))'];
temp3 = [abs(traj_train_wHwoM_mat(:, idx, rep)-traj_obs_mat(startDay:trainDays_end, idx, rep))', abs(traj_test_wHwoM_mat(:, idx, rep)-traj_obs_test_mat(testDays_start:testDays_end, idx, rep))'];
temp4 = [abs(traj_train_wHwM_mat(:, idx, rep)-traj_obs_mat(startDay:trainDays_end, idx, rep))', abs(traj_test_wHwM_mat(:, idx, rep)-traj_obs_test_mat(testDays_start:testDays_end, idx, rep))'];
% temp5 = [abs(traj_train_wHwMwGL_temp(:, idx)-traj_obs_mat(startDay:trainDays_end, idx, rep))', abs(traj_test_wHwMwGL_temp(:, idx)-traj_obs_test_mat(testDays_start:testDays_end, idx, rep))'];
temp5 = [abs(traj_train_wHwMwGL_temp(:, idx)-traj_obs_mat(startDay:trainDays_end, idx, rep))', abs(traj_test_wHwMwGL_temp(:, idx)-traj_obs_test_mat(testDays_start:testDays_end, idx, rep))'];
% temp5_group2 = [abs(traj_train_wHwMwGL_temp_group2(:, idx)-traj_obs_mat(startDay:trainDays_end, idx, rep))', abs(traj_test_wHwMwGL_temp_group2(:, idx)-traj_obs_test_mat(testDays_start:testDays_end, idx, rep))'];
% temp5_group3 = [abs(traj_train_wHwMwGL_temp_group3(:, idx)-traj_obs_mat(startDay:trainDays_end, idx, rep))', abs(traj_test_wHwMwGL_temp_group3(:, idx)-traj_obs_test_mat(testDays_start:testDays_end, idx, rep))'];
hold on;
plot(1:testDays_end, temp1, 'x-', 'LineWidth', 2, 'MarkerSize', 8);
plot(1:testDays_end, temp2, 'd-', 'LineWidth', 2, 'color', '#cc99ff', 'MarkerSize', 8);
plot(1:testDays_end, temp3, 'p-', 'LineWidth', 2, 'MarkerSize', 8);
plot(1:testDays_end, temp4, 'h-', 'LineWidth', 2, 'MarkerSize', 8);
plot(1:testDays_end, temp5, 's-', 'LineWidth', 2, 'MarkerSize', 8);
% plot(1:testDays_end, temp5_group2, 'x-', 'LineWidth', 2, 'color', '#99ff33');
% plot(1:testDays_end, temp5_group3, 'x-', 'LineWidth', 2, 'color', '#336600');

set(gca, 'FontSize', FontSize);
% set(gcf,'renderer','Painters')
title(sprintf('Province %d in Group 3', map_ind_to_ind(prov_group3_ind)), 'FontSize', 36, 'interpreter', 'latex');
xlabel('day', 'interpreter', 'latex');
ylabel('Absolute Errors', 'interpreter', 'latex');
xline(validDays_end, 'LineWidth', 2);
xline(trainDays_end, 'LineWidth', 2, 'color', '#A2142F');
grid on;
% legend('woHwoM', 'woHwM', 'wHwoM', 'wHwM', 'wHwMwGL','Location','best');
% legend('Model 1', 'Model 2', 'Model 3', 'Model 4', 'Model 5 (P)', 'Model 5 (P'')', 'Model 5 (P'''')','Location','best', 'FontSize', 24);
legend('Model 1', 'Model 2', 'Model 3', 'Model 4', 'Model 5','Location','best', 'FontSize', FontSize_legend, 'interpreter', 'latex');


% subplot(2,2,4);
% idx = 4;
% temp1 = [abs(traj_train_woHwoM_mat(:, idx, rep)-traj_obs_mat(startDay:trainDays_end, idx, rep))', abs(traj_test_woHwoM_mat(:, idx, rep)-traj_obs_test_mat(testDays_start:testDays_end, idx, rep))'];
% temp2 = [abs(traj_train_woHwM_mat(:, idx, rep)-traj_obs_mat(startDay:trainDays_end, idx, rep))', abs(traj_test_woHwM_mat(:, idx, rep)-traj_obs_test_mat(testDays_start:testDays_end, idx, rep))'];
% temp3 = [abs(traj_train_wHwoM_mat(:, idx, rep)-traj_obs_mat(startDay:trainDays_end, idx, rep))', abs(traj_test_wHwoM_mat(:, idx, rep)-traj_obs_test_mat(testDays_start:testDays_end, idx, rep))'];
% temp4 = [abs(traj_train_wHwM_mat(:, idx, rep)-traj_obs_mat(startDay:trainDays_end, idx, rep))', abs(traj_test_wHwM_mat(:, idx, rep)-traj_obs_test_mat(testDays_start:testDays_end, idx, rep))'];
% temp5 = [abs(traj_train_wHwMwGL_temp(:, idx)-traj_obs_mat(startDay:trainDays_end, idx, rep))', abs(traj_test_wHwMwGL_temp(:, idx)-traj_obs_test_mat(testDays_start:testDays_end, idx, rep))'];
% hold on;
% plot(1:testDays_end, temp1, 'x-', 'LineWidth', 2);
% plot(1:testDays_end, temp2, 'x-', 'LineWidth', 2);
% plot(1:testDays_end, temp3, 'x-', 'LineWidth', 2);
% plot(1:testDays_end, temp4, 'x-', 'LineWidth', 2);
% plot(1:testDays_end, temp5, 'x-', 'LineWidth', 2);
% set(gca, 'FontSize', 24);
% title('Province 4', 'FontSize', 30);
% xlabel('day');
% ylabel('Absolute Errors');
% xline(validDays_end, 'LineWidth', 2);
% xline(trainDays_end, 'LineWidth', 2, 'color', '#A2142F');
% grid on;
% 
% % yticks([0 2e3 4e3 6e3 8e3 1e4])
% % yticklabels({'0', '2,000','4,000','6,000', '8,000', '10,000'})
% % legend('woHwoM', 'woHwM', 'wHwoM', 'wHwM', 'wHwMwGL','Location','best');
% legend('Model 1', 'Model 2', 'Model 3', 'Model 4', 'Model 5','Location','best', 'FontSize', 24);

sgtitle({"Absolute errors of predicted trajectories for all models", ''}, 'FontSize', 42, 'interpreter', 'latex');




%%
% A1 = ones(length(group_cell{1}),length(group_cell{1}));
% A2 = ones(length(group_cell{2}),length(group_cell{2}));
% A3 = ones(length(group_cell{3}),length(group_cell{3}));
% 
% 
% W = blkdiag(A1, A2, A3);
% L = diag(sum(W,2)) - W;

% [V,D] = eig(L, 'vector');

coord = zeros(provNum, 2);
size = 0.1;
% coord(1,:) = [-1*size,3*size];
% coord(2,:) = [0*size,3*size];
% coord(3,:) = [1*size,3*size];
% coord(4,:) = [-1*size,2*size];
% coord(5,:) = [0*size,2*size];
% coord(6,:) = [1*size,2*size];
% coord(7,:) = [-1*size,1*size];
% coord(8,:) = [0*size,1*size];
% coord(9,:) = [1*size,1*size];
% 
% coord(10,:) = [-3.9*size,-1*size];
% coord(11,:) = [-3*size,-1*size];
% coord(12,:) = [-2.1*size,-1*size];
% coord(13,:) = [-1.2*size,-1*size];
% coord(14,:) = [-3.9*size,-2*size];
% coord(15,:) = [-3*size,-2*size];
% coord(16,:) = [-2.1*size,-2*size];
% coord(17,:) = [-1.2*size,-2*size];
% coord(18,:) = [-3.9*size,-3*size];
% coord(19,:) = [-3*size,-3*size];
% coord(20,:) = [-2.1*size,-3*size];
% coord(21,:) = [-1.2*size,-3*size];

% coord(22,:) = [1.2*size,-1*size];
% coord(23,:) = [2.1*size,-1*size];
% coord(24,:) = [3*size,-1*size];
% coord(25,:) = [1.2*size,-2*size];
% coord(26,:) = [2.1*size,-2*size];
% coord(27,:) = [3*size,-2*size];
% coord(28,:) = [1.2*size,-3*size];
% coord(29,:) = [2.1*size,-3*size];
% coord(30,:) = [3*size,-3*size];


coord(1,:) = [-2*size,5.5*size];
coord(2,:) = [-0.5*size,5.5*size];
coord(3,:) = [1*size,5.5*size];
coord(4,:) = [-2*size,4*size];
coord(5,:) = [-0.5*size,4*size];
coord(6,:) = [1*size,4*size];
coord(7,:) = [-2*size,2.5*size];
coord(8,:) = [-0.5*size,2.5*size];
coord(9,:) = [1*size,2.5*size];

coord(10,:) = [-5.5*size,0*size];
coord(11,:) = [-4*size,0*size];
coord(12,:) = [-2.5*size,0*size];
coord(13,:) = [-1*size,0*size];
coord(14,:) = [-5.5*size,-1.5*size];
coord(15,:) = [-4*size,-1.5*size];
coord(16,:) = [-2.5*size,-1.5*size];
coord(17,:) = [-1*size,-1.5*size];
coord(18,:) = [-5.5*size,-3*size];
coord(19,:) = [-4*size,-3*size];
coord(20,:) = [-2.5*size,-3*size];
coord(21,:) = [-1*size,-3*size];


coord(22,:) = [1.5*size,-5.5*size];
coord(23,:) = [0*size,-5.5*size];
coord(24,:) = [-1.5*size,-5.5*size];
coord(25,:) = [1.5*size,-7*size];
coord(26,:) = [0*size,-7*size];
coord(27,:) = [-1.5*size,-7*size];
coord(28,:) = [1.5*size,-8.5*size];
coord(29,:) = [0*size,-8.5*size];
coord(30,:) = [-1.5*size,-8.5*size];


% lambdavec_new = lambdavec(ind_to_ind)';


% lambdavec_temp = mean(param_wHwMwGL_mat(19,61:90,:),3);
% lambdavec_new = lambdavec_temp(ind_to_ind)';

lambdavec_temp = mean(param_woHwoM_mat(61,:));
lambdavec_new = lambdavec_temp * ones(30,1);


mark = 'ooooooooossssssssssssddddddddd';


figure(200), clf;
%--------------------------------------------------------------------------
subplot(1,3,1)
hold on
lambdavec_new = lambdavec(ind_to_ind)';
rectangle('Position',[-2.8*size 1.6*size 4.8*size 4.8*size],'Curvature',0.2, ...
    'EdgeColor', 'none', 'FaceColor', [0.8 1 1])

rectangle('Position',[-6.3*size -3.8*size 6.4*size 4.8*size],'Curvature',0.2, ...
    'EdgeColor', 'none', 'FaceColor', [0.9 1 0.9])

rectangle('Position',[-2.3*size -9.5*size 4.8*size 4.8*size],'Curvature',0.2, ...
    'EdgeColor', 'none', 'FaceColor', [1 1 0.87])

for i = 1:30
    scatter(coord(i,1), coord(i,2), 1000, lambdavec_new(i), 'filled', 'marker', mark(i));
end
% grid on
axis equal
hcb = colorbar;
hcb.FontSize = 28;
% hcb.Title.String = "$\lambda_k$";
hcb.Title.FontSize = 40;
hcb.Title.Interpreter = 'latex';
hcb.Title.Position = [92 290.7525 0];
% hcb.Limits = [0.2929 0.5144];
colormap parula
axis tight
xlim([-6.6*size 3.7*size])
ylim([-10.7*size 7.7*size])
text(coord(:,1)+0.01, coord(:,2)+0.07, cellstr(num2str((1:30)')), 'Fontsize', 30, 'FontWeight', 'bold');%, 'Interpreter', 'latex');
text(-6.3*size, 4*size, 'Group 1', 'Fontsize', 32, 'FontWeight', 'bold', 'color', '#0000cc');%, 'Interpreter', 'latex');
text(.2*size, -1.5*size, 'Group 2', 'Fontsize', 32, 'FontWeight', 'bold', 'color', '#006600');%, 'Interpreter', 'latex');
text(-5.9*size, -7*size, 'Group 3', 'Fontsize', 32, 'FontWeight', 'bold', 'color', '#cc9900');%s, 'Interpreter', 'latex');

set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])

hcb.Title.String = "$\lambda_k$";
%--------------------------------------------------------------------------

subplot(1,3,2)
hold on
lambdavec_temp = mean(param_wHwMwGL_mat(19,61:90,:),3);
lambdavec_new = lambdavec_temp(ind_to_ind)';

rectangle('Position',[-2.8*size 1.6*size 4.8*size 4.8*size],'Curvature',0.2, ...
    'EdgeColor', 'none', 'FaceColor', [0.8 1 1])

rectangle('Position',[-6.3*size -3.8*size 6.4*size 4.8*size],'Curvature',0.2, ...
    'EdgeColor', 'none', 'FaceColor', [0.9 1 0.9])

rectangle('Position',[-2.3*size -9.5*size 4.8*size 4.8*size],'Curvature',0.2, ...
    'EdgeColor', 'none', 'FaceColor', [1 1 0.87])

for i = 1:30
    scatter(coord(i,1), coord(i,2), 1000, lambdavec_new(i), 'filled', 'marker', mark(i));
end
% grid on
axis equal
hcb = colorbar;
hcb.FontSize = 28;
% hcb.Title.String = "$\lambda_k$";
hcb.Title.FontSize = 40;
hcb.Title.Interpreter = 'latex';
hcb.Title.Position = [103 290.7525 0];
% hcb.Limits = [0.2929 0.5144];
colormap parula
axis tight
xlim([-6.6*size 3.7*size])
ylim([-10.7*size 7.7*size])
text(coord(:,1)-0.0, coord(:,2)+0.06, cellstr(num2str((1:30)')), 'Fontsize', 30, 'FontWeight', 'bold');%, 'Interpreter', 'latex');
text(-6.3*size, 4*size, 'Group 1', 'Fontsize', 32, 'FontWeight', 'bold', 'color', '#0000cc');%, 'Interpreter', 'latex');
text(.2*size, -1.5*size, 'Group 2', 'Fontsize', 32, 'FontWeight', 'bold', 'color', '#006600');%, 'Interpreter', 'latex');
text(-5.9*size, -7*size, 'Group 3', 'Fontsize', 32, 'FontWeight', 'bold', 'color', '#cc9900');%s, 'Interpreter', 'latex');

set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])

hcb.Title.String = "$\widehat{\lambda}_k^{[5]}$";

%--------------------------------------------------------------------------

subplot(1,3,3)
hold on
lambdavec_temp = mean(param_woHwoM_mat(61,:));
lambdavec_new = lambdavec_temp * ones(30,1);

rectangle('Position',[-2.8*size 1.6*size 4.8*size 4.8*size],'Curvature',0.2, ...
    'EdgeColor', 'none', 'FaceColor', [0.8 1 1])

rectangle('Position',[-6.3*size -3.8*size 6.4*size 4.8*size],'Curvature',0.2, ...
    'EdgeColor', 'none', 'FaceColor', [0.9 1 0.9])

rectangle('Position',[-2.3*size -9.5*size 4.8*size 4.8*size],'Curvature',0.2, ...
    'EdgeColor', 'none', 'FaceColor', [1 1 0.87])

for i = 1:30
    scatter(coord(i,1), coord(i,2), 1000, lambdavec_new(i), 'filled', 'marker', mark(i));
end
% grid on
axis equal
hcb = colorbar;
hcb.FontSize = 28;
% hcb.Title.String = "$\lambda_k$";
hcb.Title.FontSize = 40;
hcb.Title.Interpreter = 'latex';
hcb.Title.Position = [92 290.7525 0];
% hcb.Limits = [0.2929 0.5144];
colormap parula
axis tight
xlim([-6.6*size 3.7*size])
ylim([-10.7*size 7.7*size])
text(coord(:,1)-0.0, coord(:,2)+0.06, cellstr(num2str((1:30)')), 'Fontsize', 30, 'FontWeight', 'bold');%, 'Interpreter', 'latex');
text(-6.3*size, 4*size, 'Group 1', 'Fontsize', 32, 'FontWeight', 'bold', 'color', '#0000cc');%, 'Interpreter', 'latex');
text(.2*size, -1.5*size, 'Group 2', 'Fontsize', 32, 'FontWeight', 'bold', 'color', '#006600');%, 'Interpreter', 'latex');
text(-5.9*size, -7*size, 'Group 3', 'Fontsize', 32, 'FontWeight', 'bold', 'color', '#cc9900');%s, 'Interpreter', 'latex');

set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'ytick',[])
set(gca,'yticklabel',[])

hcb.Title.String = "$\widehat{\lambda}_k^{[1]}$";










