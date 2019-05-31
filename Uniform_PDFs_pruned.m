
%% Symbolic calculation and numeric simulation for uniform distriutions with pruning
                
clear all
close all
         
%% Parameters                
N = 10000000;   %number of pair connections
a = 0.1;        %pruning value. 0.1 is good to see the effect in sum and abs_diff but for the joint we need a lower value, for instance 0.01

%% Plotting parameters
numericFontSize = 25;
axesFontSize = 30;
lineThickness = 2;
markLine = 1;
markSize = 12;

n_bins = 200;                   %bins for the histogram

%% Parameters and variables for symbolic calculation
syms u;

%% Simulation samples
x_unp = rand(N,1) .* (rand(N,1) > a);
y_unp = rand(N,1) .* (rand(N,1) > a);

x = x_unp((x_unp+y_unp)~=0);    %uses logical index to cut away the zero elements from the distributions when they have the same index in x and y (w_{ij}=w_{ji}=0)
y = y_unp((x_unp+y_unp)~=0);    %uses logical index to cut away the zero elements from the distributions when they have the same index in x and y (w_{ij}=w_{ji}=0)

%% Single distributions
sum_var = x + y;
abs_diff_var = abs(x-y);

%% Plots

% Initial Distribution
bin_dim = 1 / n_bins;
[counts,bin_location] = hist(x_unp, n_bins);               %Note: bin_location==sample
figure(1);    
bar(bin_location, counts/(sum(bin_dim*counts)));  
hold on
x_ax = 0:0.001:1;
y_ax = (1-a) * ones(1, size(x_ax,2));
h = plot(x_ax,y_ax);
set(h, 'color', 'k', 'LineWidth', lineThickness);

set(gca,'fontsize',numericFontSize);
xlabel('w','fontsize',axesFontSize);
ylabel('f(w)','fontsize',axesFontSize);
axis([-0.2 1.2 0 22]);
colormap(gray);
caxis([-1 1.4]);
title('');

print(gcf, '-depsc2', '-loose', 'Uniform_initialPDF_prune0_1'); % Print the figure in eps (first option) and uncropped (second object) 


% Sum
bin_dim = 2 / n_bins;
[counts,bin_location] = hist(sum_var, n_bins);         %Note: bin_location==sample
figure(2);    
bar(bin_location, counts/(sum(bin_dim*counts)));  
hold on
h1 = ezplot( ((1-a)/(1+a))*u + (2*a)/(1+a), [0,1] );
set(h1,'color', 'k', 'LineWidth', lineThickness);
hold on
h2 = ezplot( -((1-a)/(1+a))*u + 2*((1-a)/(1+a)), [1,2] );
set(h2,'color','k', 'LineWidth', lineThickness);

set(gca,'fontsize',numericFontSize);
xlabel('Z_2','fontsize',axesFontSize);
ylabel('f(Z_2)','fontsize',axesFontSize);
axis([-0.2 2.2 0 1.5]);
colormap(gray);
caxis([-1 1.4]);
title('');

print(gcf, '-depsc2', '-loose', 'Uniform_sumPDF_prune0_1'); % Print the figure in eps (first option) and uncropped (second object) 


% Difference
bin_dim = 2 / n_bins;
[counts,bin_location] = hist(x-y, n_bins);            %Note: bin_location==sample
figure(3);    
bar(bin_location, counts/(sum(bin_dim*counts)));  
hold on
h1 = ezplot( ((1-a)/(1+a))*u + 1/(1+a), [-1,0]);
set(h1,'color', 'k', 'LineWidth', lineThickness);
hold on
h2 = ezplot( -((1-a)/(1+a))*u + 1/(1+a), [0,1]);
set(h2,'color', 'k', 'LineWidth', lineThickness);

set(gca,'fontsize',numericFontSize);
xlabel('diff','fontsize',axesFontSize);
ylabel('f(diff)','fontsize',axesFontSize);
axis([-1.2 1.2 0 1.5]);
colormap(gray);
caxis([-1 1.4]);
title('');


% Absolute difference
bin_dim = 1 / n_bins;
[counts,bin_location] = hist(abs_diff_var, n_bins);    %Note: bin_location==sample
figure(4);    
bar(bin_location, counts/(sum(bin_dim*counts)));  
hold on
h = ezplot( -2*((1-a)/(1+a))*u + 2/(1+a), [0,1]);
set(h,'color','k', 'LineWidth', lineThickness);

set(gca,'fontsize',numericFontSize);
xlabel('Z_1','fontsize',axesFontSize);
ylabel('f(Z_1)','fontsize',axesFontSize);
axis([-0.2 1.2 0 2.5]);
colormap(gray);
caxis([-1 1.4]);
title('');

print(gcf, '-depsc2', '-loose', 'Uniform_absdiffPDF_prune0_1'); % Print the figure in eps (first option) and uncropped (second object)


% Joint Distribution 
bin_dim_x = 1 / n_bins;
bin_dim_y = 2 / n_bins;
[bidim_counts,Cell_location] = hist3([sum_var, abs_diff_var], [n_bins,n_bins]);   
figure(5);
mesh(Cell_location{1},Cell_location{2},bidim_counts'/(sum(sum(bin_dim_x*bin_dim_y*bidim_counts))))      %simulation
%hold on
%surf(Cell_location{1},Cell_location{2}, ((1-a)/(1+a))*ones(size(Cell_location{1},2),size(Cell_location{2},2))')

set(gca,'fontsize',numericFontSize);
axis on
%xlabel('Z_2','fontsize',axesFontSize);
%ylabel('Z_1','fontsize',axesFontSize);
%zlabel('f(Z_1,Z_2)','fontsize',axesFontSize);
%set(gca, 'ZTick',[0:(1-a)/(1+a):3]);
colormap(gray);
caxis([-2 6]);
view([145 45]);
title('');

print(gcf, '-depsc2', '-loose', 'Uniform_jointPDF_prune0_01'); % Print the figure in eps (first option) and uncropped (second object)

