
%% Symbolic calculation and numeric simulation for uniform distriutions with no pruning
                
clear all
close all
         
%% Parameters                
N = 10000000;   %number of pair connections

%% Plotting parameters
numericFontSize = 25;
axesFontSize = 30;
lineThickness = 2;
markLine = 1;
markSize = 12;

n_bins = 200;   %bins for the histogram

%% Parameters and variables for symbolic calculation
syms u;

%% Simulation samples
x = rand(N,1);
y = rand(N,1);

%% Single distributions
sum_var = x + y;
abs_diff_var = abs(x-y);

%% Plots

% Initial Distribution
bin_dim = 1 / n_bins;
[counts,bin_location] = hist(x, n_bins);               %Note: bin_location==sample
figure(1);    
bar(bin_location, counts/(sum(bin_dim*counts)));  
hold on
x_ax = 0:0.001:1;
y_ax = ones(1, size(x_ax,2));
h = plot(x_ax,y_ax);
set(h, 'color', 'k', 'LineWidth', lineThickness);

set(gca,'fontsize',numericFontSize);
xlabel('w','fontsize',axesFontSize);
ylabel('f(w)','fontsize',axesFontSize);
axis([-0.2 1.2 0 1.2]);
colormap(gray);
caxis([-1 1.4]);
title('');

print(gcf, '-depsc2', '-loose', 'Uniform_initialPDF_noprune'); % Print the figure in eps (first option) and uncropped (second object) 


% Sum
bin_dim = 2 / n_bins;
[counts,bin_location] = hist(sum_var, n_bins);         %Note: bin_location==sample
figure(2);    
bar(bin_location, counts/(sum(bin_dim*counts)));  
hold on
h1 = ezplot(u,[0,1]);
set(h1,'color', 'k', 'LineWidth', lineThickness);
hold on
h2 = ezplot(-u+2,[1,2]);
set(h2,'color','k', 'LineWidth', lineThickness);

set(gca,'fontsize',numericFontSize);
xlabel('Z_2','fontsize',axesFontSize);
ylabel('f(Z_2)','fontsize',axesFontSize);
axis([-0.2 2.2 0 1.5]);
colormap(gray);
caxis([-1 1.4]);
title('');

print(gcf, '-depsc2', '-loose', 'Uniform_sumPDF_noprune'); % Print the figure in eps (first option) and uncropped (second object) 
 
 
% Difference
bin_dim = 2 / n_bins;
[counts,bin_location] = hist(x-y, n_bins);            %Note: bin_location==sample
figure(3);    
bar(bin_location, counts/(sum(bin_dim*counts)));  
hold on
h1 = ezplot(u+1,[-1,0]);
set(h1,'color', 'k', 'LineWidth', lineThickness);
hold on
h2 = ezplot(-u+1,[0,1]);
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
h = ezplot(-2*u+2,[0,1]);
set(h,'color','k', 'LineWidth', lineThickness);

set(gca,'fontsize',numericFontSize);
xlabel('Z_1','fontsize',axesFontSize);
ylabel('f(Z_1)','fontsize',axesFontSize);
axis([-0.2 1.2 0 2.5]);
colormap(gray);
caxis([-1 1.4]);
title('');

print(gcf, '-depsc2', '-loose', 'Uniform_absdiffPDF_noprune'); % Print the figure in eps (first option) and uncropped (second object) 


% Joint Distribution
bin_dim_x = 1 / n_bins;
bin_dim_y = 2 / n_bins;
[bidim_counts,Cell_location] = hist3([sum_var, abs_diff_var], [n_bins,n_bins]);   
figure(5);
mesh(Cell_location{1},Cell_location{2},bidim_counts'/(sum(sum(bin_dim_x*bin_dim_y*bidim_counts))))      %simulation
%hold on
%surf(Cell_location{1},Cell_location{2}, ones(size(Cell_location{1},2),size(Cell_location{2},2))')

set(gca,'fontsize',numericFontSize);
axis on
% xlabel('Z_2','fontsize',axesFontSize);
% ylabel('Z_1','fontsize',axesFontSize);
% zlabel('f(Z_1,Z_2)','fontsize',axesFontSize);
colormap(gray);
caxis([-1.5 4]);
view([145 45]);
title('');

print(gcf, '-depsc2', '-loose', 'Uniform_jointPDF_noprune'); % Print the figure in eps (first option) and uncropped (second object) 

