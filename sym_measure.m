
%% Numerical simulation for the symmetry measure. Uniform distribution

function [s] = sym_measure (matrix)

%%
%Parameters and variables necessary for running this code without external call. Also a loop over the value of a is needed
%n_samples = 10000;                                            %number of matrixes in the sample
%N = 200;                                                      %number of neurons
%max_w = 1;
%a = 0:0.1:0.9;                                                %pruning values
%number_points = size(a,2);                                    %number of points in the plot
%sample_mean = zeros(number_points);
%sample_variance = zeros(number_points);

%%
%symm = zeros(1,n_samples);

%for iter = 1:n_samples
    
    %sample_matrix = max_w .* rand(N) .* (rand(N) > a);      %generate a random NxN matrix from zero to max_w and introduce pruning a
        
        upper = triu(matrix,1);                      %extract the upper triangle matrix
        lower = tril(matrix,-1)';                    %extract the lower triangle matrix and transpose it
        
        x = upper(:);                                       %convert the matrix into a vector
        y = lower(:);                                       %convert the matrix into a vector
        
        temp = x + y;                                       %sum vector elements==sum the reciprocal elements of the matrix
        
        nonzero_index = find(temp~=0.);                     %create a vector whoose elements are the index of the non zero elements in temp
        
        K = length(nonzero_index);                          %counts how many elements of temp are nonzero==counts the number of pairs connections for which at least one direction is nonzero
        
        if K > 0
            s = 1 - sum ( abs(x(nonzero_index)-y(nonzero_index)) ./ (x(nonzero_index)+y(nonzero_index)) ) / K;
        else
            s = 0;
        end
        
        %sprintf('Point number %d iteration number %d',n,iter)
    
%end

%sample_mean = mean(symm);
%sample_variance = var(symm);


    
    
   