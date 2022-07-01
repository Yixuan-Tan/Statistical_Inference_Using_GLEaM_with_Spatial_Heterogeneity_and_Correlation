% confirmed_data = readmatrix("Confirmed.xlsx");
% deaths_data    = readmatrix("Deaths.xlsx");
% recovered_data = readmatrix("Recovered.xlsx");

% confirmed_data(:,1) = [];
% deaths_data(:,1)    = [];
% recovered_data(:,1) = [];

clear all

confirmed_data = readtable("Confirmed.xlsx");
deaths_data    = readtable("Deaths.xlsx");
recovered_data = readtable("Recovered.xlsx");

country_names    = convertCharsToStrings(confirmed_data.Properties.VariableNames);
country_names(1) = [];
country_names    = country_names';

dates_string = string(confirmed_data.Country_Region);

confirmed_data = removevars(confirmed_data,'Country_Region');
deaths_data    = removevars(deaths_data,'Country_Region');
recovered_data = removevars(recovered_data,'Country_Region');

accu_confirmed_data = table2array(confirmed_data);
accu_recovered_data = table2array(recovered_data) + table2array(deaths_data);

clear confirmed_data deaths_data recovered_data

save('data_preprocessed.mat');

clear all




