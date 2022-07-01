% four provinces

% validation set last for the same time period as the testing set, for example:
% training   set: 1-10days
% validation set: 11-30days
% testing    set: 11-30 days
% trajectories from which these sets are from are all independent

% predicted testing trajectories are computed in a "hidden way":
% using the fitted [I0, E0, lambda], compute the trajectory from the 1st
% day instead of using the ground truth information from the 10th day
% namely, the only known facts are newly confirmed cases

clear all;

% rng(621);
rng(100000);

%%

provNum = 4;
provPop = 1e6 * ones(1, provNum);

traVal  = 5e3; % prov 1 to prov 2 and prov 2 to prov 1
traMat  = traVal * ones(provNum);
traMat  = traMat - diag(diag(traMat)); % traval from prov i to prov i = 0

delta     = 0.14;
gamma     = 0.14;
lambdavec = [0.5, 0.47, 0.4, 0.37] % prov 1 & 2; prov 3 & 4
% lambdavec = [0.4, 0.37, 0.35, 0.32]

E0      = 30 * ones(1, provNum);
I0      = 10 * ones(1, provNum);
param   = [I0, E0, lambdavec];

totDays  = 20;
startDay = 1;
trainDays_end   = 10;
validDays_start = trainDays_end+1;
validDays_end   = totDays;
testDays_start  = trainDays_end+1;
testDays_end    = totDays;

ini.totPop = provPop;
%ini.traVal = traVal;
ini.part   = cell(2, 1);
% ini.part{1} = [1;2];
% ini.part{2} = [3;4];
ini.part{1} = [1;3];
ini.part{2} = [2;4];


ini.delta  = delta;
ini.gamma  = gamma;
%ini.lambda = lambda;

dt      = 0.1;
%day_end = 25;

%%
numrep = 1e2;
sigma  = 10^(-6);
% sigma  = 10^(-3);
% sigma  = 1;

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
muvec = 10.^ (1.5:0.1:3.5)
% muvec = 10.^(-2);% * 5000;
beta  = 0.1; % mu0 = beta * mu1
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

trainErr_woHwoM_mat2  = zeros(numrep,1); % relative error with weighted average
trainErr_woHwM_mat2   = zeros(numrep,1);
trainErr_wHwoM_mat2   = zeros(numrep,1);
trainErr_wHwM_mat2    = zeros(numrep,1);
trainErr_wHwMwGL_mat2 = zeros(length(muvec), numrep);

validErr_woHwoM_mat2  = zeros(numrep,1); 
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
    
    y_true_test  = max(simu_data_generate_random_approx(param, ini, totDays, startDay, dt),1);
    y_true_valid = max(simu_data_generate_random_approx(param, ini, totDays, startDay, dt),1);%y_true;%max(SEIR_data_gen_test_pred_inf_random_dt_2(param, ini, totDays, startDay, dt),1);
    y_true       = max(simu_data_generate_random_approx(param, ini, totDays, startDay, dt),1);
    
    y_obs   = y_true(:, 1: (5*provNum)); % deltaC, deltaCR,  E, I, R
    traj_obs_mat(:, :, rep) = y_obs;
    traj_obs_test_mat(:, :, rep) = y_true_test(:, 1: (5*provNum));
    traj_obs_valid_mat(:, :, rep) = y_true_valid(:, 1: (5*provNum));
    
    %%
    %-----------------------------------------------------------------------------------------------------
    % wo Hetero, wo Migrat
    fprintf('wo Hetero, wo Migrat\n');
    ini.traMat = zeros(provNum);
    
    f1 = @(x)fmin_simu_data_without_hetero(x, y_obs, ini, 'pois', dt, trainDays_end, startDay );
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
    y_reg_train = simu_data_generate_determ_without_hetero(params1, ini, dt, trainDays_end-startDay+1, 1, R_train);
    y_reg_valid = simu_data_generate_determ_without_hetero(params1_valid, ini, dt, validDays_end-startDay+1, 1, R_valid);
    y_reg_valid = y_reg_valid(validDays_start:validDays_end,:);
    y_reg_test  = simu_data_generate_determ_without_hetero(params1_test, ini, dt, testDays_end-startDay+1, 1, R_test);
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
    
    f1 = @(x)fmin_simu_data_without_hetero(x, y_obs, ini, 'pois', dt, trainDays_end, startDay );
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
    y_reg_train = simu_data_generate_determ_without_hetero(params1, ini, dt, trainDays_end-startDay+1, 1, R_train);
    y_reg_valid = simu_data_generate_determ_without_hetero(params1_valid, ini, dt, validDays_end-startDay+1, 1, R_valid);
    y_reg_valid = y_reg_valid(validDays_start:validDays_end,:);
    y_reg_test  = simu_data_generate_determ_without_hetero(params1_test, ini, dt, testDays_end-startDay+1, 1, R_test);
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
    ini.traMat = zeros(provNum); 
    
    f1 = @(x)fmin_simu_data(x, y_obs, ini, 'pois', 0, beta, dt, trainDays_end, startDay, 0);
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
    y_reg_train = simu_data_generate_determ(params1, ini, dt, trainDays_end-startDay+1, 1, R_train);
    y_reg_valid = simu_data_generate_determ(params1_valid, ini, dt, validDays_end-startDay+1, 1, R_valid);
    y_reg_valid = y_reg_valid(validDays_start:validDays_end,:);
    y_reg_test  = simu_data_generate_determ(params1_test, ini, dt, testDays_end-startDay+1, 1, R_test);
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
    ini.traMat = traMat; 
    
    f1 = @(x)fmin_simu_data(x, y_obs, ini, 'pois', 0, beta, dt, trainDays_end, startDay, 0);
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
    y_reg_train = simu_data_generate_determ(params1, ini, dt, trainDays_end-startDay+1, 1, R_train);
    y_reg_valid = simu_data_generate_determ(params1_valid, ini, dt, validDays_end-startDay+1, 1, R_valid);
    y_reg_valid = y_reg_valid(validDays_start:validDays_end,:);
    y_reg_test  = simu_data_generate_determ(params1_test, ini, dt, testDays_end-startDay+1, 1, R_test);
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
        
        mu = muvec(i); %mu1vec(j);
        f1 = @(x)fmin_simu_data(x, y_obs, ini, 'pois', mu, beta, dt, trainDays_end, startDay, sigma );
                
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
        y_reg_train = simu_data_generate_determ(params1, ini, dt, trainDays_end-startDay+1, 1, R_train);
        y_reg_valid = simu_data_generate_determ(params1_valid, ini, dt, validDays_end-startDay+1, 1, R_valid);
        y_reg_valid = y_reg_valid(validDays_start:validDays_end,:);
        y_reg_test  = simu_data_generate_determ(params1_test, ini, dt, testDays_end-startDay+1, 1, R_test);
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

muind_plot = 13;

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

sgtitle({'Mean of total weighted relative validation errors over 100 replicas', '(simulated data with 4 provinces)'}, 'Fontsize', 42, 'Interpreter', 'latex');

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

sgtitle({'Mean of total weighted relative testing errors over 100 replicas', '(simulated data with 4 provinces)'}, 'Fontsize', 42, 'Interpreter', 'latex');



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
ylabel('testing error', 'Interpreter', 'latex');
set(gca, 'FontSize', 36);
title('$\textnormal{MSE}^{[\rm Val]}_{(s)}$', 'FontSize', 36, 'Interpreter', 'latex');
grid on;

sgtitle({'Mean of total simple averaged relative validation errors over 100 replicas', '(simulated data with 4 provinces)'}, 'Fontsize', 42, 'Interpreter', 'latex');

figure(8), clf;
subplot(1,2,1);
hold on;
plot(log10(muvec), mean(testErr_wHwMwGL_mat2,2), 'x-', 'LineWidth', 2);
plot(log10(muvec(muind_plot)), mean(testErr_wHwMwGL_mat2(muind_plot,:),2), 'p', 'MarkerSize', 24, 'MarkerFaceColor', '#EDB120', 'MarkerEdgeColor', '#EDB120');
yline(mean(testErr_wHwM_mat2), '-', '$\mu=0$', 'LineWidth', 2, 'FontSize', 30, 'LabelHorizontalAlignment', 'left', 'Interpreter', 'latex');
xlabel('$\log_{10} \mu$', 'Interpreter', 'latex');
ylabel('testing error', 'Interpreter', 'latex');
set(gca, 'FontSize', 36);
% ylim([0
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

sgtitle({'Mean of total simply averaged relative testing errors over 100 replicas', '(simulated data with 4 provinces)'}, 'Fontsize', 42, 'Interpreter', 'latex');


% return;


%%

trainDays_end
validDays_end
totDays

fprintf('relative error with weighted average \n\n');

fprintf('wo H, wo M, trainErr: mean %6.4f, std % 6.4f, testErr: mean %6.4f, std % 6.4f \n', mean(trainErr_woHwoM_mat1),...
    std(trainErr_woHwoM_mat1), mean(testErr_woHwoM_mat1), std(testErr_woHwoM_mat1));
fprintf('wo H, w/ M, trainErr: mean %6.4f, std % 6.4f, testErr: mean %6.4f, std % 6.4f \n', mean(trainErr_woHwM_mat1),...
    std(trainErr_woHwM_mat1), mean(testErr_woHwM_mat1), std(testErr_woHwM_mat1));
fprintf('w/ H, wo M, trainErr: mean %6.4f, std % 6.4f, testErr: mean %6.4f, std % 6.4f \n', mean(trainErr_wHwoM_mat1),...
    std(trainErr_wHwoM_mat1), mean(testErr_wHwoM_mat1), std(testErr_wHwoM_mat1));
fprintf('w/ H, w/ M, trainErr: mean %6.4f, std % 6.4f, testErr: mean %6.4f, std % 6.4f \n', mean(trainErr_wHwM_mat1),...
    std(trainErr_wHwM_mat1), mean(testErr_wHwM_mat1), std(testErr_wHwM_mat1));
fprintf('\n');
for i = 1: length(muvec)
    fprintf('w/ H, w/ M, w/ GL, mu = %6.1f, trainErr: mean %.4f, std %6.4f, testErr: mean %6.4f, std %.4f \n', muvec(i), mean(trainErr_wHwMwGL_mat1(i,:)),...
        std(trainErr_wHwMwGL_mat1(i,:)),  mean(testErr_wHwMwGL_mat1(i,:)), std(testErr_wHwMwGL_mat1(i,:)));
end

%%

fprintf('relative error with equal weights \n\n');

fprintf('wo H, wo M, trainErr: mean %6.4f, std % 6.4f, testErr: mean %6.4f, std % 6.4f \n', mean(trainErr_woHwoM_mat2),...
    std(trainErr_woHwoM_mat2), mean(testErr_woHwoM_mat2), std(testErr_woHwoM_mat2));
fprintf('wo H, w/ M, trainErr: mean %6.4f, std % 6.4f, testErr: mean %6.4f, std % 6.4f \n', mean(trainErr_woHwM_mat2),...
    std(trainErr_woHwM_mat2), mean(testErr_woHwM_mat2), std(testErr_woHwM_mat2));
fprintf('w/ H, wo M, trainErr: mean %6.4f, std % 6.4f, testErr: mean %6.4f, std % 6.4f \n', mean(trainErr_wHwoM_mat2),...
    std(trainErr_wHwoM_mat2), mean(testErr_wHwoM_mat2), std(testErr_wHwoM_mat2));
fprintf('w/ H, w/ M, trainErr: mean %6.4f, std % 6.4f, testErr: mean %6.4f, std % 6.4f \n', mean(trainErr_wHwM_mat2),...
    std(trainErr_wHwM_mat2), mean(testErr_wHwM_mat2), std(testErr_wHwM_mat2));
fprintf('\n');
for i = 1: length(muvec)
    %fprintf('mu0 = %6.4f: \n\n', mu0vec(i));
    %     for j = 1: length(mu1vec)
    %         fprintf('w/ H, w/ M, w/ GL, mu1 = %6.4f, trainErr: mean %6.4f, std % 6.4f, testErr: mean %6.4f, std % 6.4f \n', mu1vec(j), mean(trainErr_wHwMwGL_mat(i,j,:)),...
    %             std(trainErr_wHwMwGL_mat(i,j,:)),  mean(testErr_wHwMwGL_mat(i,j,:)), std(testErr_wHwMwGL_mat(i,j,:)));
    %     end
    
    fprintf('w/ H, w/ M, w/ GL, mu1 = %6.4f, trainErr: mean %6.4f, std % 6.4f, testErr: mean %6.4f, std % 6.4f \n', muvec(i), mean(trainErr_wHwMwGL_mat2(i,:)),...
        std(trainErr_wHwMwGL_mat2(i,:)),  mean(testErr_wHwMwGL_mat2(i,:)), std(testErr_wHwMwGL_mat2(i,:)));
    
end

%%

fprintf('\nrelative MSE with weighted average\n\n');

fprintf('wo H, wo M, trainErr: mean %6.4f, std % 6.4f, testErr: mean %6.4f, std % 6.4f \n', mean(trainErr_woHwoM_mat3),...
    std(trainErr_woHwoM_mat3), mean(testErr_woHwoM_mat3), std(testErr_woHwoM_mat3));
fprintf('wo H, w/ M, trainErr: mean %6.4f, std % 6.4f, testErr: mean %6.4f, std % 6.4f \n', mean(trainErr_woHwM_mat3),...
    std(trainErr_woHwM_mat3), mean(testErr_woHwM_mat3), std(testErr_woHwM_mat3));
fprintf('w/ H, wo M, trainErr: mean %6.4f, std % 6.4f, testErr: mean %6.4f, std % 6.4f \n', mean(trainErr_wHwoM_mat3),...
    std(trainErr_wHwoM_mat3), mean(testErr_wHwoM_mat3), std(testErr_wHwoM_mat3));
fprintf('w/ H, w/ M, trainErr: mean %6.4f, std % 6.4f, testErr: mean %6.4f, std % 6.4f \n', mean(trainErr_wHwM_mat3),...
    std(trainErr_wHwM_mat3), mean(testErr_wHwM_mat3), std(testErr_wHwM_mat3));
fprintf('\n');
for i = 1: length(muvec)
    %fprintf('mu0 = %6.4f: \n\n', mu0vec(i));
    %     for j = 1: length(mu1vec)
    %         fprintf('w/ H, w/ M, w/ GL, mu1 = %6.4f, trainErr: mean %6.4f, std % 6.4f, testErr: mean %6.4f, std % 6.4f \n', mu1vec(j), mean(trainErr_wHwMwGL_mat(i,j,:)),...
    %             std(trainErr_wHwMwGL_mat(i,j,:)),  mean(testErr_wHwMwGL_mat(i,j,:)), std(testErr_wHwMwGL_mat(i,j,:)));
    %     end
    
    fprintf('w/ H, w/ M, w/ GL, mu1 = %6.4f, trainErr: mean %6.4f, std % 6.4f, testErr: mean %6.4f, std % 6.4f \n', muvec(i), mean(trainErr_wHwMwGL_mat3(i,:)),...
        std(trainErr_wHwMwGL_mat3(i,:)),  mean(testErr_wHwMwGL_mat3(i,:)), std(testErr_wHwMwGL_mat3(i,:)));
    
end

%%
fprintf('\nrelative MSE with simple average\n\n');

fprintf('wo H, wo M, trainErr: mean %6.4f, std % 6.4f, testErr: mean %6.4f, std % 6.4f \n', mean(trainErr_woHwoM_mat4),...
    std(trainErr_woHwoM_mat4), mean(testErr_woHwoM_mat4), std(testErr_woHwoM_mat4));
fprintf('wo H, w/ M, trainErr: mean %6.4f, std % 6.4f, testErr: mean %6.4f, std % 6.4f \n', mean(trainErr_woHwM_mat4),...
    std(trainErr_woHwM_mat4), mean(testErr_woHwM_mat4), std(testErr_woHwM_mat4));
fprintf('w/ H, wo M, trainErr: mean %6.4f, std % 6.4f, testErr: mean %6.4f, std % 6.4f \n', mean(trainErr_wHwoM_mat4),...
    std(trainErr_wHwoM_mat4), mean(testErr_wHwoM_mat4), std(testErr_wHwoM_mat4));
fprintf('w/ H, w/ M, trainErr: mean %6.4f, std % 6.4f, testErr: mean %6.4f, std % 6.4f \n', mean(trainErr_wHwM_mat4),...
    std(trainErr_wHwM_mat4), mean(testErr_wHwM_mat4), std(testErr_wHwM_mat4));
fprintf('\n');
for i = 1: length(muvec)
    %fprintf('mu0 = %6.4f: \n\n', mu0vec(i));
    %     for j = 1: length(mu1vec)
    %         fprintf('w/ H, w/ M, w/ GL, mu1 = %6.4f, trainErr: mean %6.4f, std % 6.4f, testErr: mean %6.4f, std % 6.4f \n', mu1vec(j), mean(trainErr_wHwMwGL_mat(i,j,:)),...
    %             std(trainErr_wHwMwGL_mat(i,j,:)),  mean(testErr_wHwMwGL_mat(i,j,:)), std(testErr_wHwMwGL_mat(i,j,:)));
    %     end
    
    fprintf('w/ H, w/ M, w/ GL, mu1 = %6.4f, trainErr: mean %6.4f, std % 6.4f, testErr: mean %6.4f, std % 6.4f \n', muvec(i), mean(trainErr_wHwMwGL_mat4(i,:)),...
        std(trainErr_wHwMwGL_mat4(i,:)),  mean(testErr_wHwMwGL_mat4(i,:)), std(testErr_wHwMwGL_mat4(i,:)));
    
end


return


%%


fprintf('estimated lambda \n\n');

fprintf('wo H, wo M, lambda: mean %6.4f, std % 6.4f \n', mean(param_woHwoM_mat(9,:)), std(param_woHwoM_mat(9,:)));
fprintf('wo H, w/ M, lambda: mean %6.4f, std % 6.4f \n', mean(param_woHwM_mat(9,:)), std(param_woHwM_mat(9,:)));
fprintf('w/ H, wo M, lambda1: mean %6.4f, std % 6.4f, lambda2: mean %6.4f, std % 6.4f \n', ...
    mean(param_wHwoM_mat(9,:)), std(param_wHwoM_mat(9,:)), mean(param_wHwoM_mat(10,:)), std(param_wHwoM_mat(10,:)));
fprintf('            lambda3: mean %6.4f, std % 6.4f, lambda4: mean %6.4f, std % 6.4f \n', ...
    mean(param_wHwoM_mat(11,:)), std(param_wHwoM_mat(11,:)), mean(param_wHwoM_mat(12,:)), std(param_wHwoM_mat(12,:)));
fprintf('w/ H, w/ M, lambda1: mean %6.4f, std % 6.4f, lambda2: mean %6.4f, std % 6.4f \n', ...
    mean(param_wHwM_mat(9,:)), std(param_wHwM_mat(9,:)), mean(param_wHwM_mat(10,:)), std(param_wHwM_mat(10,:)));
fprintf('            lambda3: mean %6.4f, std % 6.4f, lambda4: mean %6.4f, std % 6.4f \n', ...
    mean(param_wHwM_mat(11,:)), std(param_wHwM_mat(11,:)), mean(param_wHwM_mat(12,:)), std(param_wHwM_mat(12,:)));
fprintf('\n');
for i = 1: length(muvec)
    %fprintf('mu0 = %6.4f: \n\n', mu0vec(i));
    %     for j = 1: length(mu1vec)
    %         fprintf('w/ H, w/ M, w/ GL, mu1 = %6.4f, trainErr: mean %6.4f, std % 6.4f, testErr: mean %6.4f, std % 6.4f \n', mu1vec(j), mean(trainErr_wHwMwGL_mat(i,j,:)),...
    %             std(trainErr_wHwMwGL_mat(i,j,:)),  mean(testErr_wHwMwGL_mat(i,j,:)), std(testErr_wHwMwGL_mat(i,j,:)));
    %     end
    
    fprintf('w/ H, w/ M, wGL, mu1 %6.4f, lambda1: mean %6.4f, std % 6.4f, lambda2: mean %6.4f, std % 6.4f \n', muvec(i), ...
        mean(reshape(param_wHwMwGL_mat(i,9,:), numrep,1)), std(reshape(param_wHwMwGL_mat(i,9,:), numrep,1)), ...
        mean(reshape(param_wHwMwGL_mat(i,10,:), numrep,1)), std(reshape(param_wHwMwGL_mat(i,10,:), numrep,1)));
    fprintf('                            lambda3: mean %6.4f, std % 6.4f, lambda4: mean %6.4f, std % 6.4f \n', ...
        mean(reshape(param_wHwMwGL_mat(i,11,:), numrep,1)), std(reshape(param_wHwMwGL_mat(i,11,:), numrep,1)), ...
        mean(reshape(param_wHwMwGL_mat(i,12,:), numrep,1)), std(reshape(param_wHwMwGL_mat(i,12,:), numrep,1)));
    
    
    
end


%%


rep = 50;


FontSize = 36;
FontSize_legend = 36;


mu1ind = 13;
y_true_test = reshape(traj_obs_test_mat(:,:,rep),totDays,5*provNum);
params1 = param_wHwMwGL_mat(mu1ind, :, rep);
R_train = zeros(1, provNum);
params1_test = params1;
params1_test = param;
R_test = zeros(1, provNum);%y_true_test(testDays_start-1, (4*provNum+1):(5*provNum)); % change R
params1_valid = params1;
R_valid = zeros(1, provNum);%y_true_test(validDays_start-1, (4*provNum+1):(5*provNum)); % change R

traj_train_wHwMwGL_temp = SEIR_data_gen_test_pred_inf_determ_2(params1, ini, dt, trainDays_end-startDay+1, 1, R_train);
traj_test_wHwMwGL_temp  = SEIR_data_gen_test_pred_inf_determ_2(params1_test, ini, dt, testDays_end-startDay+1, 1, R_test);
traj_test_wHwMwGL_temp  = traj_test_wHwMwGL_temp(testDays_start:testDays_end,:);
traj_valid_wHwMwGL_temp = SEIR_data_gen_test_pred_inf_determ_2(params1_valid, ini, dt, validDays_end-startDay+1, 1, R_valid);
traj_valid_wHwMwGL_temp = traj_valid_wHwMwGL_temp(validDays_start:validDays_end,:);


figure(3), clf;
subplot(2,2,1);
idx = 1;
temp1 = [traj_train_woHwoM_mat(:, idx, rep)', traj_test_woHwoM_mat(:, idx, rep)'];
temp2 = [traj_train_woHwM_mat(:, idx, rep)', traj_test_woHwM_mat(:, idx, rep)'];
temp3 = [traj_train_wHwoM_mat(:, idx, rep)', traj_test_wHwoM_mat(:, idx, rep)'];
temp4 = [traj_train_wHwM_mat(:, idx, rep)', traj_test_wHwM_mat(:, idx, rep)'];
temp5 = [traj_train_wHwMwGL_temp(:, idx)', traj_test_wHwMwGL_temp(:, idx)'];
temp6 = [traj_obs_mat(startDay:trainDays_end, idx, rep)', traj_obs_test_mat(testDays_start:testDays_end, idx, rep)'];
hold on;
plot(1:testDays_end, temp1, 'x-', 'LineWidth', 2, 'MarkerSize', 8);
plot(1:testDays_end, temp2, 'd-', 'LineWidth', 2, 'color', '#cc99ff', 'MarkerSize', 8);
plot(1:testDays_end, temp3, 'p-', 'LineWidth', 2, 'MarkerSize', 8);
plot(1:testDays_end, temp4, 'h-', 'LineWidth', 2, 'MarkerSize', 8);
plot(1:testDays_end, temp5, 's-', 'LineWidth', 2, 'MarkerSize', 8);
plot(1:testDays_end, temp6, 'o-', 'LineWidth', 2, 'color', '#4DBEEE', 'MarkerSize', 8);
set(gca, 'FontSize', FontSize);
title('Province 1', 'FontSize', 36, 'interpreter', 'latex');
xlabel('day', 'interpreter', 'latex'); ylabel('cases', 'interpreter', 'latex');
xline(validDays_end, 'LineWidth', 2);
xline(trainDays_end, 'LineWidth', 2, 'color', '#A2142F');
grid on;
legend('Model 1', 'Model 2', 'Model 3', 'Model 4', 'Model 5','True','Location','best', 'FontSize', FontSize_legend, 'interpreter', 'latex');
%

subplot(2,2,2);
idx = 2;
temp1 = [traj_train_woHwoM_mat(:, idx, rep)', traj_test_woHwoM_mat(:, idx, rep)'];
temp2 = [traj_train_woHwM_mat(:, idx, rep)', traj_test_woHwM_mat(:, idx, rep)'];
temp3 = [traj_train_wHwoM_mat(:, idx, rep)', traj_test_wHwoM_mat(:, idx, rep)'];
temp4 = [traj_train_wHwM_mat(:, idx, rep)', traj_test_wHwM_mat(:, idx, rep)'];
temp5 = [traj_train_wHwMwGL_temp(:, idx)', traj_test_wHwMwGL_temp(:, idx)'];
temp6 = [traj_obs_mat(startDay:trainDays_end, idx, rep)', traj_obs_test_mat(testDays_start:testDays_end, idx, rep)'];
hold on;
plot(1:testDays_end, temp1, 'x-', 'LineWidth', 2, 'MarkerSize', 8);
plot(1:testDays_end, temp2, 'd-', 'LineWidth', 2, 'color', '#cc99ff', 'MarkerSize', 8);
plot(1:testDays_end, temp3, 'p-', 'LineWidth', 2, 'MarkerSize', 8);
plot(1:testDays_end, temp4, 'h-', 'LineWidth', 2, 'MarkerSize', 8);
plot(1:testDays_end, temp5, 's-', 'LineWidth', 2, 'MarkerSize', 8);
plot(1:testDays_end, temp6, 'o-', 'LineWidth', 2, 'color', '#4DBEEE', 'MarkerSize', 8);
set(gca, 'FontSize', FontSize);
title('Province 2', 'FontSize', 36, 'interpreter', 'latex');
xlabel('day', 'interpreter', 'latex'); ylabel('cases', 'interpreter', 'latex');
xline(validDays_end, 'LineWidth', 2);
xline(trainDays_end, 'LineWidth', 2, 'color', '#A2142F');
grid on;
legend('Model 1', 'Model 2', 'Model 3', 'Model 4', 'Model 5','True','Location','best', 'FontSize', FontSize_legend, 'interpreter', 'latex');

subplot(2,2,3);
idx = 3;
temp1 = [traj_train_woHwoM_mat(:, idx, rep)', traj_test_woHwoM_mat(:, idx, rep)'];
temp2 = [traj_train_woHwM_mat(:, idx, rep)', traj_test_woHwM_mat(:, idx, rep)'];
temp3 = [traj_train_wHwoM_mat(:, idx, rep)', traj_test_wHwoM_mat(:, idx, rep)'];
temp4 = [traj_train_wHwM_mat(:, idx, rep)', traj_test_wHwM_mat(:, idx, rep)'];
temp5 = [traj_train_wHwMwGL_temp(:, idx)', traj_test_wHwMwGL_temp(:, idx)'];
temp6 = [traj_obs_mat(startDay:trainDays_end, idx, rep)', traj_obs_test_mat(testDays_start:testDays_end, idx, rep)'];
hold on;
plot(1:testDays_end, temp1, 'x-', 'LineWidth', 2, 'MarkerSize', 8);
plot(1:testDays_end, temp2, 'd-', 'LineWidth', 2, 'color', '#cc99ff', 'MarkerSize', 8);
plot(1:testDays_end, temp3, 'p-', 'LineWidth', 2, 'MarkerSize', 8);
plot(1:testDays_end, temp4, 'h-', 'LineWidth', 2, 'MarkerSize', 8);
plot(1:testDays_end, temp5, 's-', 'LineWidth', 2, 'MarkerSize', 8);
plot(1:testDays_end, temp6, 'o-', 'LineWidth', 2, 'color', '#4DBEEE', 'MarkerSize', 8);
set(gca, 'FontSize', FontSize);
title('Province 3', 'FontSize', 36, 'interpreter', 'latex');
xlabel('day', 'interpreter', 'latex'); ylabel('cases', 'interpreter', 'latex');
xline(validDays_end, 'LineWidth', 2);
xline(trainDays_end, 'LineWidth', 2, 'color', '#A2142F');
grid on;
legend('Model 1', 'Model 2', 'Model 3', 'Model 4', 'Model 5','True','Location','best', 'FontSize', FontSize_legend, 'interpreter', 'latex');

subplot(2,2,4);
idx = 4;
temp1 = [traj_train_woHwoM_mat(:, idx, rep)', traj_test_woHwoM_mat(:, idx, rep)'];
temp2 = [traj_train_woHwM_mat(:, idx, rep)', traj_test_woHwM_mat(:, idx, rep)'];
temp3 = [traj_train_wHwoM_mat(:, idx, rep)', traj_test_wHwoM_mat(:, idx, rep)'];
temp4 = [traj_train_wHwM_mat(:, idx, rep)', traj_test_wHwM_mat(:, idx, rep)'];
temp5 = [traj_train_wHwMwGL_temp(:, idx)', traj_test_wHwMwGL_temp(:, idx)'];
temp6 = [traj_obs_mat(startDay:trainDays_end, idx, rep)', traj_obs_test_mat(testDays_start:testDays_end, idx, rep)'];
hold on;
plot(1:testDays_end, temp1, 'x-', 'LineWidth', 2, 'MarkerSize', 8);
plot(1:testDays_end, temp2, 'd-', 'LineWidth', 2, 'color', '#cc99ff', 'MarkerSize', 8);
plot(1:testDays_end, temp3, 'p-', 'LineWidth', 2, 'MarkerSize', 8);
plot(1:testDays_end, temp4, 'h-', 'LineWidth', 2, 'MarkerSize', 8);
plot(1:testDays_end, temp5, 's-', 'LineWidth', 2, 'MarkerSize', 8);
plot(1:testDays_end, temp6, 'o-', 'LineWidth', 2, 'color', '#4DBEEE', 'MarkerSize', 8);
set(gca, 'FontSize', FontSize);
title('Province 4', 'FontSize', 36, 'interpreter', 'latex');
xlabel('day', 'interpreter', 'latex'); ylabel('cases', 'interpreter', 'latex');
xline(validDays_end, 'LineWidth', 2);
xline(trainDays_end, 'LineWidth', 2, 'color', '#A2142F');
grid on;
legend('Model 1', 'Model 2', 'Model 3', 'Model 4', 'Model 5','True','Location','best', 'FontSize', FontSize_legend, 'interpreter', 'latex');

sgtitle("True trajectory \& Fitted trajectories", 'FontSize', 42, 'interpreter', 'latex');

%

%rep = 83;
figure(4), clf;
subplot(2,2,1);
idx = 1;
temp1 = [abs(traj_train_woHwoM_mat(:, idx, rep)-traj_obs_mat(startDay:trainDays_end, idx, rep))', abs(traj_test_woHwoM_mat(:, idx, rep)-traj_obs_test_mat(testDays_start:testDays_end, idx, rep))'];
temp2 = [abs(traj_train_woHwM_mat(:, idx, rep)-traj_obs_mat(startDay:trainDays_end, idx, rep))', abs(traj_test_woHwM_mat(:, idx, rep)-traj_obs_test_mat(testDays_start:testDays_end, idx, rep))'];
temp3 = [abs(traj_train_wHwoM_mat(:, idx, rep)-traj_obs_mat(startDay:trainDays_end, idx, rep))', abs(traj_test_wHwoM_mat(:, idx, rep)-traj_obs_test_mat(testDays_start:testDays_end, idx, rep))'];
temp4 = [abs(traj_train_wHwM_mat(:, idx, rep)-traj_obs_mat(startDay:trainDays_end, idx, rep))', abs(traj_test_wHwM_mat(:, idx, rep)-traj_obs_test_mat(testDays_start:testDays_end, idx, rep))'];
temp5 = [abs(traj_train_wHwMwGL_temp(:, idx)-traj_obs_mat(startDay:trainDays_end, idx, rep))', abs(traj_test_wHwMwGL_temp(:, idx)-traj_obs_test_mat(testDays_start:testDays_end, idx, rep))'];
hold on;
plot(1:testDays_end, temp1, 'x-', 'LineWidth', 2, 'MarkerSize', 8);
plot(1:testDays_end, temp2, 'd-', 'LineWidth', 2, 'color', '#cc99ff', 'MarkerSize', 8);
plot(1:testDays_end, temp3, 'p-', 'LineWidth', 2, 'MarkerSize', 8);
plot(1:testDays_end, temp4, 'h-', 'LineWidth', 2, 'MarkerSize', 8);
plot(1:testDays_end, temp5, 's-', 'LineWidth', 2, 'MarkerSize', 8);
set(gca, 'FontSize', FontSize);
title('Province 1', 'FontSize', 36, 'interpreter', 'latex');
xlabel('day', 'interpreter', 'latex');
ylabel('Absolute Errors', 'interpreter', 'latex');
xline(validDays_end, 'LineWidth', 2);
xline(trainDays_end, 'LineWidth', 2, 'color', '#A2142F');
grid on;
legend('Model 1', 'Model 2', 'Model 3', 'Model 4', 'Model 5','Location','best', 'FontSize', FontSize_legend, 'interpreter', 'latex');

subplot(2,2,2);
idx = 2;
temp1 = [abs(traj_train_woHwoM_mat(:, idx, rep)-traj_obs_mat(startDay:trainDays_end, idx, rep))', abs(traj_test_woHwoM_mat(:, idx, rep)-traj_obs_test_mat(testDays_start:testDays_end, idx, rep))'];
temp2 = [abs(traj_train_woHwM_mat(:, idx, rep)-traj_obs_mat(startDay:trainDays_end, idx, rep))', abs(traj_test_woHwM_mat(:, idx, rep)-traj_obs_test_mat(testDays_start:testDays_end, idx, rep))'];
temp3 = [abs(traj_train_wHwoM_mat(:, idx, rep)-traj_obs_mat(startDay:trainDays_end, idx, rep))', abs(traj_test_wHwoM_mat(:, idx, rep)-traj_obs_test_mat(testDays_start:testDays_end, idx, rep))'];
temp4 = [abs(traj_train_wHwM_mat(:, idx, rep)-traj_obs_mat(startDay:trainDays_end, idx, rep))', abs(traj_test_wHwM_mat(:, idx, rep)-traj_obs_test_mat(testDays_start:testDays_end, idx, rep))'];
temp5 = [abs(traj_train_wHwMwGL_temp(:, idx)-traj_obs_mat(startDay:trainDays_end, idx, rep))', abs(traj_test_wHwMwGL_temp(:, idx)-traj_obs_test_mat(testDays_start:testDays_end, idx, rep))'];
hold on;
plot(1:testDays_end, temp1, 'x-', 'LineWidth', 2, 'MarkerSize', 8);
plot(1:testDays_end, temp2, 'd-', 'LineWidth', 2, 'color', '#cc99ff', 'MarkerSize', 8);
plot(1:testDays_end, temp3, 'p-', 'LineWidth', 2, 'MarkerSize', 8);
plot(1:testDays_end, temp4, 'h-', 'LineWidth', 2, 'MarkerSize', 8);
plot(1:testDays_end, temp5, 's-', 'LineWidth', 2, 'MarkerSize', 8);
set(gca, 'FontSize', FontSize);
title('Province 2', 'FontSize', 36, 'interpreter', 'latex');
xlabel('day', 'interpreter', 'latex');
ylabel('Absolute Errors', 'interpreter', 'latex');
xline(validDays_end, 'LineWidth', 2);
xline(trainDays_end, 'LineWidth', 2, 'color', '#A2142F');
grid on;
legend('Model 1', 'Model 2', 'Model 3', 'Model 4', 'Model 5','Location','best', 'FontSize', FontSize_legend, 'interpreter', 'latex');

subplot(2,2,3);
idx = 3;
temp1 = [abs(traj_train_woHwoM_mat(:, idx, rep)-traj_obs_mat(startDay:trainDays_end, idx, rep))', abs(traj_test_woHwoM_mat(:, idx, rep)-traj_obs_test_mat(testDays_start:testDays_end, idx, rep))'];
temp2 = [abs(traj_train_woHwM_mat(:, idx, rep)-traj_obs_mat(startDay:trainDays_end, idx, rep))', abs(traj_test_woHwM_mat(:, idx, rep)-traj_obs_test_mat(testDays_start:testDays_end, idx, rep))'];
temp3 = [abs(traj_train_wHwoM_mat(:, idx, rep)-traj_obs_mat(startDay:trainDays_end, idx, rep))', abs(traj_test_wHwoM_mat(:, idx, rep)-traj_obs_test_mat(testDays_start:testDays_end, idx, rep))'];
temp4 = [abs(traj_train_wHwM_mat(:, idx, rep)-traj_obs_mat(startDay:trainDays_end, idx, rep))', abs(traj_test_wHwM_mat(:, idx, rep)-traj_obs_test_mat(testDays_start:testDays_end, idx, rep))'];
temp5 = [abs(traj_train_wHwMwGL_temp(:, idx)-traj_obs_mat(startDay:trainDays_end, idx, rep))', abs(traj_test_wHwMwGL_temp(:, idx)-traj_obs_test_mat(testDays_start:testDays_end, idx, rep))'];
hold on;
plot(1:testDays_end, temp1, 'x-', 'LineWidth', 2, 'MarkerSize', 8);
plot(1:testDays_end, temp2, 'd-', 'LineWidth', 2, 'color', '#cc99ff', 'MarkerSize', 8);
plot(1:testDays_end, temp3, 'p-', 'LineWidth', 2, 'MarkerSize', 8);
plot(1:testDays_end, temp4, 'h-', 'LineWidth', 2, 'MarkerSize', 8);
plot(1:testDays_end, temp5, 's-', 'LineWidth', 2, 'MarkerSize', 8);
set(gca, 'FontSize', FontSize);
title('Province 3', 'FontSize', 36, 'interpreter', 'latex');
xlabel('day', 'interpreter', 'latex');
ylabel('Absolute Errors', 'interpreter', 'latex');
xline(validDays_end, 'LineWidth', 2);
xline(trainDays_end, 'LineWidth', 2, 'color', '#A2142F');
grid on;
legend('Model 1', 'Model 2', 'Model 3', 'Model 4', 'Model 5','Location','best', 'FontSize', FontSize_legend, 'interpreter', 'latex');

subplot(2,2,4);
idx = 4;
temp1 = [abs(traj_train_woHwoM_mat(:, idx, rep)-traj_obs_mat(startDay:trainDays_end, idx, rep))', abs(traj_test_woHwoM_mat(:, idx, rep)-traj_obs_test_mat(testDays_start:testDays_end, idx, rep))'];
temp2 = [abs(traj_train_woHwM_mat(:, idx, rep)-traj_obs_mat(startDay:trainDays_end, idx, rep))', abs(traj_test_woHwM_mat(:, idx, rep)-traj_obs_test_mat(testDays_start:testDays_end, idx, rep))'];
temp3 = [abs(traj_train_wHwoM_mat(:, idx, rep)-traj_obs_mat(startDay:trainDays_end, idx, rep))', abs(traj_test_wHwoM_mat(:, idx, rep)-traj_obs_test_mat(testDays_start:testDays_end, idx, rep))'];
temp4 = [abs(traj_train_wHwM_mat(:, idx, rep)-traj_obs_mat(startDay:trainDays_end, idx, rep))', abs(traj_test_wHwM_mat(:, idx, rep)-traj_obs_test_mat(testDays_start:testDays_end, idx, rep))'];
temp5 = [abs(traj_train_wHwMwGL_temp(:, idx)-traj_obs_mat(startDay:trainDays_end, idx, rep))', abs(traj_test_wHwMwGL_temp(:, idx)-traj_obs_test_mat(testDays_start:testDays_end, idx, rep))'];
hold on;
plot(1:testDays_end, temp1, 'x-', 'LineWidth', 2, 'MarkerSize', 8);
plot(1:testDays_end, temp2, 'd-', 'LineWidth', 2, 'color', '#cc99ff', 'MarkerSize', 8);
plot(1:testDays_end, temp3, 'p-', 'LineWidth', 2, 'MarkerSize', 8);
plot(1:testDays_end, temp4, 'h-', 'LineWidth', 2, 'MarkerSize', 8);
plot(1:testDays_end, temp5, 's-', 'LineWidth', 2, 'MarkerSize', 8);
set(gca, 'FontSize', FontSize);
title('Province 4', 'FontSize', 36, 'interpreter', 'latex');
xlabel('day', 'interpreter', 'latex');
ylabel('Absolute Errors', 'interpreter', 'latex');
xline(validDays_end, 'LineWidth', 2);
xline(trainDays_end, 'LineWidth', 2, 'color', '#A2142F');
grid on;
legend('Model 1', 'Model 2', 'Model 3', 'Model 4', 'Model 5','Location','best', 'FontSize', FontSize_legend, 'interpreter', 'latex');

sgtitle("Absolute errors of predicted trajectories for all models", 'FontSize', 42, 'interpreter', 'latex');


