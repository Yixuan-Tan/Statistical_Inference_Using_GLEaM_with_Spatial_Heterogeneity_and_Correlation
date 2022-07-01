clear all;
load("data_preprocessed.mat");

tot_days = length(dates_string);

map_country_to_ind = containers.Map(country_names, (1:11));

newly_confirmed_data = diff(accu_confirmed_data);
newly_recovered_data = diff(accu_recovered_data);


%% taking average

% generate a circulant matrix
% averaged newly confirmed / recovered cases

c      = zeros(tot_days-7,1);
c(1)   = 1;
r      = zeros(tot_days-1,1);
r(1:7) = ones(7,1);
T      = toeplitz(c,r);
T      = T/7;

newly_confirmed_data_ave = T * newly_confirmed_data;
newly_recovered_data_ave = T * newly_recovered_data;

% generate a circulant matrix
% averaged accumulated confirmed / recovered cases

c      = zeros(tot_days-6,1);
c(1)   = 1;
r      = zeros(tot_days,1);
r(1:7) = ones(7,1);
T      = toeplitz(c,r);
T      = T/7;

accu_confirmed_data_ave = T * accu_confirmed_data;
accu_recovered_data_ave = T * accu_recovered_data;

% change date

dates_string_ave = dates_string;
dates_string_ave(1:3) = [];
dates_string_ave(end-2:end) = [];

%%

figure(1), clf;

for i = 1:11
    subplot(3,4,i);
    country_name_temp = country_names(i);
    country_ind_temp  = map_country_to_ind(country_name_temp);
    plot(1:(tot_days-7), newly_confirmed_data_ave(:,country_ind_temp), 'x-');
    
    grid on;

    title(sprintf('%s', country_name_temp));
    set(gca, 'FontSize', 14);
end
sgtitle('newly confirmed cases with average', 'FontSize', 20);



figure(2), clf;

for i = 1:11
    subplot(3,4,i);
    country_name_temp = country_names(i);
    country_ind_temp  = map_country_to_ind(country_name_temp);
    plot(1:(tot_days-7), newly_recovered_data_ave(:,country_ind_temp), 'x-');
    
    grid on;

    title(sprintf('%s', country_name_temp));
    set(gca, 'FontSize', 14);
end
sgtitle('newly recovered cases with average', 'FontSize', 20);

%%

clear c r i T country_ind_temp country_name_temp
save('data_averaged.mat');
clear all;



























