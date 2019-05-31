
%% Mean and variance of symmetry measure for uniform distribution and for:
% 1. Random network (simulation and theory)
% 2. Random symmetric network (simulation)
% 3. Random asymmetric network (simulation)
% 4. Random fully symmetric network (simulation)
% 5. Random fully asymmetric network (simulation)

% This file calls the function sym_measure.m
% This file saves workspace in uniform.mat

close all
clear all

%% Parameters
n_samples = 100000; %number of networks
n_neurons = 10;     %number of neurons
a = 0:0.1:0.9;      %pruning values
n_points = size(a,2);
n_bins = 200;
max_w = 1;          %maximum weights value 

%% Plotting parameters
numericFontSize = 25;
axesFontSize = 30;
lineThickness = 2;
markLine = 1;
markSize = 12;

%% Variables
s_rand = zeros(1,n_samples);
s_full_sym = zeros(1,n_samples);
s_full_asym = zeros(1,n_samples);
s_rand_sym = zeros(1,n_samples);
s_rand_asym = zeros(1,n_samples);

sample_mean = zeros(1,n_points);
sample_variance = zeros(1,n_points);
full_sym_mean = zeros(1,n_points);
full_sym_variance = zeros(1,n_points);
full_asym_mean = zeros(1,n_points);
full_asym_variance = zeros(1,n_points);
rand_sym_mean = zeros(1,n_points);
rand_sym_variance = zeros(1,n_points);
rand_asym_mean = zeros(1,n_points);
rand_asym_variance = zeros(1,n_points);

theoretical_mean = zeros(1,n_points);
symbolic_mean = zeros(1,n_points);
theoretical_variance = zeros(1,n_points);

syms u v; %symbolic variables

%% Sampling
sym_matrix = max_w * (ones(n_neurons) - diag(diag(ones(n_neurons))));   %fully symmetric
asym_matrix = triu(sym_matrix,1);                                       %fully asymmetric

asym_lower = 0.1;
asym_range = max_w - asym_lower;
asym_scale = 0.01 * asym_lower;

for indx = 1:n_points
    
    %% Simulation
    for sample = 1:n_samples
    
        rand_sample = max_w .* rand(n_neurons) .* (rand(n_neurons) > a(indx));
        rand_sample = rand_sample - diag(diag(rand_sample));
        s_rand(sample) = sym_measure (rand_sample);    
        
        full_sym_sample = sym_matrix .* (rand(n_neurons) > a(indx));
        full_sym_sample = full_sym_sample - diag(diag(full_sym_sample));
        s_full_sym(sample) = sym_measure (full_sym_sample);    
        
        full_asym_sample = asym_matrix .* (rand(n_neurons) > a(indx));
        full_asym_sample = full_asym_sample - diag(diag(full_asym_sample));
        s_full_asym(sample) = sym_measure (full_asym_sample);    
    
        rand_sym_sample = triu( max_w * rand(n_neurons), 1 );
        rand_sym_sample = (rand_sym_sample + rand_sym_sample') .* (rand(n_neurons) > a(indx));
        rand_sym_sample = rand_sym_sample - diag(diag(rand_sym_sample));
        s_rand_sym(sample) = sym_measure (rand_sym_sample);  
    
        temp = rand(n_neurons);
        rand_asym_sample =  triu (asym_lower + asym_range * rand(n_neurons), 1);
        rand_asym_sample = rand_asym_sample .* (temp >= (asym_range/2)) + 0.001 * rand_asym_sample .* (temp < (asym_range/2));
        rand_asym_sample = rand_asym_sample + tril( (rand_asym_sample' >= asym_lower) .* asym_scale .* rand(n_neurons) + (rand_asym_sample' < asym_lower) .* (asym_lower + asym_range * rand(n_neurons)), -1);
        rand_asym_sample = rand_asym_sample .* (rand(n_neurons) > a(indx));
        rand_asym_sample = rand_asym_sample - diag(diag(rand_asym_sample));
        s_rand_asym(sample) = sym_measure (rand_asym_sample); 
    
    end
    
    sample_mean(indx) = mean(s_rand);
    sample_variance(indx) = var(s_rand);
    full_sym_mean(indx) = mean(s_full_sym);
    full_sym_variance(indx) = var(s_full_sym);
    full_asym_mean(indx) = mean(s_full_asym);
    full_asym_variance(indx) = var(s_full_asym);
    rand_sym_mean(indx) = mean(s_rand_sym);
    rand_sym_variance(indx) = var(s_rand_sym);
    rand_asym_mean(indx) = mean(s_rand_asym);
    rand_asym_variance(indx) = var(s_rand_asym);   
    
    %quantities for plotting distribution of s
    if indx == 1 
        s_rand_nop = s_rand;    %store no_prune case
    end
    if indx == 5
        s_rand_p04 = s_rand;    %store prune=0.4 case
    end 
    
    %quantities for plotting connectivity matrix
    if indx == 3
        w_rand_p02 = rand_sample;       %store prune=0.2 case only one trial (the last)
        s_rand_p02 = s_rand(n_samples); %store prune=0.2 case only one trial (the last)
    end
       
    %% Theoretical calculation 
     
    theoretical_mean(indx) = 1 - ( ((1-a(indx))/(1+a(indx))) * (2*log(2)-1) ) - 2*a(indx)/(1+a(indx));
    
    symbolic_mean(indx) = 1 - ( ((1-a(indx))/(1+a(indx))) * (1-theoretical_mean(1)) ) - (2*a(indx)/(1+a(indx))) * vpa(int(int( dirac(u-v), v, u-0.000000001, 2-u), u, 0, 1),6);
    %NOTE: because the delta falls right on the border of the integration domain this create some troubles. We need to add a -eps. 
    %Also, we should multiply the delta by (u/v). However, despite the result should not change because the delta sets u=v, it changes in matlab. Try with wolfram alpha to verify that it does not change
    
    theoretical_variance(indx) = 1 + ( (8*log(2)*(a(indx)-1)-7*a(indx)+5) / (1+a(indx)) ) - ( ( (2-2*log(2)+2*a(indx)*log(2)-2*a(indx)) / (1+a(indx)) ) ^2 );
    %---NOTE
    %Simulation: we first apply the measure definition, which means that we compute the mean of Z_{ij} = |w_{ij}-w_{ij}|/(w_{ij}+w_{ij}) over the network and then we compute the mean and variance over the n_samples networks
    %Theoretical calculus: is the opposite. The formula derived is by performing the mean of one conneciton Z_{ij} over the n_samples networks; therefore, after we have to mean over the connections in a network, which
    %by definition means to divide by number of connecitons (which in turn is simply apply the definition of the measure)
    n_connections = (1-a(indx)^2) * 0.5 * n_neurons * (n_neurons-1);
    theoretical_variance(indx) = theoretical_variance(indx) / n_connections;
    
    display('Iteration done');
    
end


%% a=0 check
syms u;

prune_index = 1;    %no pruning

f = exp(-((u-theoretical_mean(prune_index))^2/(2*theoretical_variance(prune_index)))) / sqrt(2*pi*theoretical_variance(prune_index));

%normalization and theory-simulation agreement check
norm = double(int(f,0,1));

sprintf('Normalization: %f.', norm)

figure(1);

subplot(1,2,1)
[counts,bin_location] = hist(s_rand_nop,n_bins);
bin_dim_hist = ( max(bin_location) - min(bin_location) ) / n_bins;
bar(bin_location, counts/(sum(bin_dim_hist*counts)));
hold on
ezplot(f);
set(gca,'fontsize',numericFontSize);
xlabel('s','fontsize',axesFontSize);
ylabel('PDF(s)','fontsize',axesFontSize);
title('');

subplot(1,2,2)
hist(s_rand_p04,100);
set(gca,'fontsize',numericFontSize);
xlabel('s','fontsize',axesFontSize);
ylabel('PDF(s)','fontsize',axesFontSize);


%% P-values
syms u;
p_up = zeros(1, size(a, 2));
p_lo = zeros(1, size(a, 2));
p_value_area_check = zeros(1, size(a, 2));

for prune_index = 1 : 9
    
    p_value_area_th = 0.9500;
    p_value_tol = 0.0001;
    
    f = exp(-((u-theoretical_mean(prune_index))^2/(2*theoretical_variance(prune_index)))) / sqrt(2*pi*theoretical_variance(prune_index));
    
    p_up(prune_index) = (1.386*sqrt(2)*sqrt(theoretical_variance(prune_index))) + theoretical_mean(prune_index); %we are using the properties of translation of a normal gaussian and the quintile function
    p_lo(prune_index) = theoretical_mean(prune_index) - (p_up(prune_index)-theoretical_mean(prune_index));
    p_value_area_check(prune_index) = double(vpa(int(f,p_lo(prune_index),p_up(prune_index)),6));

    if abs(p_value_area_check(prune_index) - p_value_area_th) > p_value_tol
        display('Error in the estimated area corresponding to the p-value');   
    else
        display('P-value has been computed correctly');   
    end
end

%% Data saving and plot

save uniform

figure(2);

L1=errorbar(a, theoretical_mean, sqrt(theoretical_variance));
set(L1, 'LineWidth', lineThickness+2);
set(L1, 'LineStyle', '--');
set(L1, 'Color', [1 1 1] * 0.);
hold on

L2=errorbar(a, sample_mean, sqrt(sample_variance));
set(L2, 'LineWidth', lineThickness);
set(L2, 'LineStyle', '-');
set(L2, 'Color', [1 1 1] * 0.5);
hold on

% L3=errorbar(a, full_sym_mean, sqrt(full_sym_variance));
% set(L3, 'LineWidth', lineThickness);
% hold on

L4=errorbar(a, rand_sym_mean, sqrt(rand_sym_variance));
set(L4, 'LineWidth', lineThickness);
set(L4, 'LineStyle', '-.');
set(L4, 'Color', [1 1 1] * 0.7);
hold on

% L5=errorbar(a, full_asym_mean, sqrt(full_asym_variance));
% set(L5, 'LineWidth', lineThickness);
% hold on

L6=errorbar(a, rand_asym_mean, sqrt(rand_asym_variance));
set(L6, 'LineWidth', lineThickness);
set(L6, 'LineStyle', ':');
set(L6, 'Color', [1 1 1] * 0.7);

set(gca,'fontsize',numericFontSize);
xlabel('a','fontsize',axesFontSize);
ylabel('E[s]','fontsize',axesFontSize);
axis([-0.05 0.95 -0.05 1.05]);

print(gcf, '-depsc2', '-loose', 'Uniform_stat'); % Print the figure in eps (first option) and uncropped (second object)


figure(3);

imagesc(w_rand_p02);
axis square;
colormap(gray);
colorbar;
set(gca,'fontsize',numericFontSize);
xlabel('Neuron post','fontsize',axesFontSize);
ylabel('Neuron pre','fontsize',axesFontSize);

print(gcf, '-depsc2', '-loose', 'Uniform_random_connectivity_example'); % Print the figure in eps (first option) and uncropped (second object)

sprintf('Symmetry measure of the sample network shown in Figure 3 (a=0.2): %f.',s_rand_p02)


