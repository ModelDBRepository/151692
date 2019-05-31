
%% Correlation

function [correlation] = correl (n_samples, mean_value, standard_deviation, n_connections)

%a = 0:0.1:0.9;
a=0;
%correlation = zeros(1, size(a,2));

%for j = 1:size(a,2) 
    
    corr = zeros(1, n_samples);    
    for i =1:n_samples
    
        x_unp = normrnd(mean_value, standard_deviation, 1, n_connections) .* (rand(1,n_connections) > a);          %Also: mean_value + (standard_deviation).*randn(1,N);   
        y_unp = normrnd(mean_value, standard_deviation, 1, n_connections) .* (rand(1,n_connections) > a);          %Also: mean_value + (standard_deviation).*randn(1,N);
    
        x = x_unp((x_unp+y_unp)~=0);    %uses logical index to cut away the zero elements from the distributions when they have the same index in x and y (w_{ij}=w_{ji}=0)
        y = y_unp((x_unp+y_unp)~=0);
    
        distr       = [      x+y;       abs(x-y)];
        mean_val    = [mean(x+y), mean(abs(x-y))];
        variance    = [ var(x+y),  var(abs(x-y))];
    
        covariance = sum((distr(1,:)-mean_val(1)) .* (distr(2,:)-mean_val(2))) / n_connections;
        corr(i) = covariance / (sqrt(variance(1))*sqrt(variance(2)));
    
    end
    
    correlation = mean(corr);    

%    correlation(j) = mean(corr);    
%    display ('Iteration done')
%end