
%% Symbolic calculation and numeric simulation for gaussian distriutions with no pruning

clear all
close all

%% Parameters
N = 10000000;               %number of pair connections
mean_value = 0.5;           %mean value of the gaussian distribution of connections
standard_deviation = 1/10;  %sqrt of variance. With 1/10 we are considering up to 5*sd. 3*sd is 99.7%

%% Plotting parameters
numericFontSize = 25;
axesFontSize = 30;
lineThickness = 2;
markLine = 1;
markSize = 12;

minv      = [ 0,    0,    -1,    0];
maxv      = [ 1,    2,     1,    1];
n_bins   = 200;      
bin_dim  = (maxv-minv) / n_bins;

%% Parameters and variables for symbolic calculation
syms u v w          %symbolic variables

mu = [mean_value,            2*mean_value,            0,                       0                     ];       %values used to create the functions: single, sum, difference, absolute difference
va = [standard_deviation^2,  2*standard_deviation^2,  2*standard_deviation^2,  2*standard_deviation^2];       %values used to create the functions: single, sum, difference, absolute difference

mu_sym = zeros(1, size(mu,2));      %mean value obtained by symbolic integration
va_sym = zeros(1, size(va,2));      %variance obtained by symbolic integration
mu_bid_sym = zeros(2);              %mean value obtained by symbolic integration
va_bid_sym = zeros(2);              %variance obtained by symbolic integration

tol=0.99;    %tolerance for the integration of the distribution

%% Simulation samples
x = normrnd(mean_value, standard_deviation, 1, N);         
y = normrnd(mean_value, standard_deviation, 1, N);         

distr    = [      x;       x+y;       x-y;       abs(x-y)]; %vector of distributions
mean_val = [mean(x), mean(x+y), mean(x-y), mean(abs(x-y))]; %corresponding mean value
variance = [ var(x),  var(x+y),  var(x-y),  var(abs(x-y))]; %corresponding variance

%% Printing vectors
variable = char('x, y', 'x+y', 'x-y', '|x-y|'); %creates a matrix where each string is spread over a row, every character being an element of the matrix
x_axis_name = char ('w','Z_2','diff','Z_1');
y_axis_name = char ('f(w)','f(Z_2)','f(diff)','f(Z_1)');
figure_name = char ('Gaussian_initialPDF_noprune','Gaussian_sumPDF_noprune','Gaussian_diffPDF_noprune','Gaussian_absdiffPDF_noprune');

%% Single distributions
for i = 1:size(distr,1)
    
    sample = minv(i):bin_dim(i):maxv(i);
    
    %---SYMBOLIC
    f = exp(-((u-mu(i))^2/(2*va(i)))) / sqrt(2*pi*va(i));
    area = int(f,minv(i),maxv(i));
    
    %area check, normalization and sampling
    if double(area) > 1
        display('Error, something is going wrong: the integral of a distribution cannot be larger than 1')
    elseif double(area) < tol
        f = f / double(area);             
        val = normpdf(sample, mu(i), sqrt(va(i))) / double(area);
        area = int(f,minv(i),maxv(i));
    else
        val = normpdf(sample, mu(i), sqrt(va(i)));
    end    
    
    %mean and variance by integration - for the absolute difference they have to be different from the ones used to create the function
    mu_sym(i) = int(f*u,minv(i),maxv(i));
    va_sym(i) = int(f*(u-mu_sym(i))^2,minv(i),maxv(i));
    
    
    %---SIMULATION
    [counts,bin_location] = hist(distr(i,:),n_bins);        %Note: bin_location==sample
    bin_dim_hist = ( max(bin_location) - min(bin_location) ) / n_bins;
    
    %---PRINT & PLOT
    sprintf('Variable: %s \nSimulation with N=%d pair connections: Mean %f. Variance %f. \nSymbolic calculation: Integral of the distribution with mean %f and variance %f: %f', variable(i,:), N, mean_val(i), variance(i), mu_sym(i), va_sym(i), double(area))
            %double(area) converts sym to a variable double, so it calculates the value and it can be displayed. Otherwise I could also type vpa(area,6) directly on a sym object: it evaluates the numerical expression
    
    figure(i);
    bar(bin_location, counts/(sum(bin_dim_hist*counts)));     %simulation --- Note: sum(bin_dim(i)*counts)==bin_dim*N
    hold on
    h = plot(sample,val);                                   %theoretical
    set(h, 'color', 'k', 'LineWidth', lineThickness);
    
    set(gca,'fontsize',numericFontSize);
    xlabel(x_axis_name(i, find(x_axis_name(i,:) ~= ' ')),'fontsize',axesFontSize);
    ylabel(y_axis_name(i, find(y_axis_name(i,:) ~= ' ')),'fontsize',axesFontSize);
    axis([minv(i)-0.5 maxv(i)+0.5 0 ceil(max( counts/(sum(bin_dim_hist*counts)) ))]);    
    
    colormap(gray);
    caxis([-1 1.4]);
            
    if i ~= 3
        print(gcf, '-depsc2', '-loose', figure_name(i, find(figure_name(i,:) ~= ' '))); % Print the figure in eps (first option) and uncropped (second object)    
    end
end


%% Joint distributions-bivariate gaussian

for i = 3:4
    
    samplex = minv(2):bin_dim(2):maxv(2);   %first distribution is the sum
    sampley = minv(i):bin_dim(i):maxv(i);   %second distribution is either the difference or the absolute difference
    
    %covariance and correlation from the simulation: the two variables are not independents
    covariance = sum ((distr(2,:)-mean_val(2)).*(distr(i,:)-mean_val(i))) / N;
    corr = covariance / (sqrt(variance(2))*sqrt(variance(i)));
    sigma = [va(2), corr*sqrt(va(2))*sqrt(va(i)); corr*sqrt(va(2))*sqrt(va(i)), va(i)]; 
    
    %---SYMBOLIC
    f = exp( - (((v-mu(2))^2/(va(2))) + ((w-mu(i))^2/(va(i))) - 2*corr*(v-mu(2))*(w-mu(i))/(sqrt(va(2))*sqrt(va(i))) ) / (2*(1-(corr^2))) )  / ( 2*pi*sqrt (va(2)*va(i)*(1-(corr^2))) );   %bivariate depedent
    %f = exp(-((v-mu(2))^2/(2*va(2)))-((w-mu(i))^2/(2*va(i))))/(2*pi*sqrt(va(2)*va(i)));     %bivariate independent
    area = int((int(f,v,minv(2),maxv(2))),w,minv(i),maxv(i));
    
    %area check, normalization and sampling
    if double(area) > 1
        display('Error, something is going wrong: the integral of a distribution cannot be larger than 1')
    elseif double(area) < tol
        f = f / double(area);
        [X1,X2] = meshgrid(samplex,sampley);                                %Generates the grid from the two vectors                    
        val = mvnpdf([X1(:) X2(:)], [mu(2),mu(i)], sigma) / double(area);   %Note: mvnpdf returns a vector evaluating the PDF at each row of [X1(:) X2(:)]
        val = reshape(val,length(samplex),length(sampley));
        area = int((int(f,v,minv(2),maxv(2))),w,minv(i),maxv(i));
    else
        [X1,X2] = meshgrid(samplex,sampley);                     
        val = mvnpdf([X1(:) X2(:)], [mu(2),mu(i)], sigma);                  %Note: mvnpdf returns a vector evaluating the PDF at each row of [X1(:) X2(:)]
        val = reshape(val,length(samplex),length(sampley));
    end
    
    %mean and variance by integration - for the absolute difference they have to be different from the ones used to create the function
    mu_bid_sym(i-2,1) = int((int(f*v,v,minv(2),maxv(2))),w,minv(i),maxv(i));
    mu_bid_sym(i-2,2) = int((int(f*w,v,minv(2),maxv(2))),w,minv(i),maxv(i));
    va_bid_sym(i-2,1) = int((int(f*(v-mu_bid_sym(i-2,1))^2,v,minv(2),maxv(2))),w,minv(i),maxv(i));
    va_bid_sym(i-2,2) = int((int(f*(w-mu_bid_sym(i-2,2))^2,v,minv(2),maxv(2))),w,minv(i),maxv(i));
          
    
    %---SIMULATION
    grid{1} = samplex;
    grid{2} = sampley;
    [bidim_counts,Cell_location] = hist3([(distr(2,:))', (distr(i,:))'], grid);
    
    %---PRINT AND PLOT
    sprintf('Variables: %s, %s \nSimulation with N=%d pair connections: Mean %f, %f. Variance %f, %f. \nSymbolic calculation: Integral of the distribution with mean %f, %f and variance %f, %f: %f', variable(2,:), variable(i,:), N, mean_val(2), mean_val(i), variance(2), variance(i), mu_bid_sym(i-2,1), mu_bid_sym(i-2,2), va_bid_sym(i-2,1), va_bid_sym(i-2,2), double(area))
           %double(area) converts sym to a variable double, so it calculates the value and it can be displayed. Otherwise I could also type vpa(area,6) directly on a sym object: it evaluates the numerical expression
    
    figure;
    mesh(Cell_location{1},Cell_location{2},bidim_counts'/(sum(sum(bin_dim(2)*bin_dim(i)*bidim_counts))))    %simulation
    hold on
    surf(samplex,sampley,val);                                                                              %theoretical
    
    set(gca,'fontsize',numericFontSize);
    %xlabel('Z_2','fontsize',axesFontSize);
    %ylabel('Z_1','fontsize',axesFontSize);
    %zlabel('f(Z_1,Z_2)','fontsize',axesFontSize);
    axis([0 2 0 1]);
   
    colormap(gray);
    caxis([-1 20]);
    view([130 14]);
       
    if i == 4
        print(gcf, '-depsc2', '-loose', 'Gaussian_jointPDF_noprune'); % Print the figure in eps (first option) and uncropped (second object)    
    end      
end


%% Check the validity of the restriction in [0,1]
%Verify that all the points of the distribution lie in the range [0,1]x[abs(x-y),2-abs(x-y)] and that the integral in this range is 1
%Few points can lie outside, those that in the initial distributions x and y are outside [0,1]

%---SIMULATION
k=0;    %counts the points inside
m=0;    %counts the points outside
for i=1:size(distr,2)
    if and(((distr(2,i))>=(distr(4,i))),((distr(2,i))<=2-(distr(4,i))))
        k=k+1;            
    else
        m=m+1;        
    end
end

%---SYMBOLIC
area_check = int((int(f,v,w,(2-w))),w,0,1);

%---PRINT
sprintf('Checking the joint distribution f(x+y,|x-y|) in the definition interval [0,1]x[|x-y|,2-|x-y|]. \nSimulation. Number of points of the distribution outside: %d out of %f. \nSymbolic calculation. Area of the joint distribution: %f', m, N, double(area_check))

