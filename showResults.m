function []=showResults(SROCCs,PLCCs,RMESs,metric)
mu_CC=mean(PLCCs);
mu_SROCC=mean(SROCCs);
mu_RMSE=mean(RMESs);
med_CC=median(PLCCs);
med_SROCC=median(SROCCs);
med_RMSE=median(RMESs);
result_mean=[metric ': ' '|' num2str(mu_CC) '|' num2str(mu_SROCC) '|'  num2str(mu_RMSE) '|']
result_med=[metric ': ' '|' num2str(med_CC) '|' num2str(med_SROCC) '|'  num2str(med_RMSE) '|']