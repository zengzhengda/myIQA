function [best_C,best_gamma]=SVR_choosing_paremeter_no_norm(IQA_feature,IQA_label,dataset_name)
% svr²ÎÊýÑ¡Ôñ
% [IQA_label, IQA_feature] = libsvmread('traindata_scaled');
training_label_vector1=IQA_label;
training_instance_matrix=IQA_feature;
folds = 5;

[C,gamma] = meshgrid(0:2:8, -8:2:0); 

%# grid search, and cross-validation
cv_acc = zeros(numel(C),1);
for i=1:numel(C)
    cv_acc(i) = svmtrain(training_label_vector1, training_instance_matrix,  sprintf('-s %f -t %f -c %f -g %f -v %d', 3, 2, 2^C(i), 2^gamma(i), folds));
end

%# pair (C,gamma) with best accuracy
[~,idx] = min(cv_acc);

%# contour plot of paramter selection
contour(C, gamma, reshape(cv_acc,size(C))), colorbar
hold on
plot(C(idx), gamma(idx), 'rx')
text(C(idx), gamma(idx), sprintf('Acc = %.2f %%',cv_acc(idx)), 'HorizontalAlign','left', 'VerticalAlign','top')
hold off
xlabel('log_2(C)'), ylabel('log_2(\gamma)'), title('Cross-Validation Accuracy')

%# now you can train you model using best_C and best_gamma
best_C = 2^C(idx),   %6
best_gamma = 2^gamma(idx),   %-6

save best_svr_param_cid2013_robust_no_norm.mat best_C best_gamma;
