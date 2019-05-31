
%% Symbolic calculation and numeric simulation for gaussian distriutions with pruning

clear all
close all

%% Parameters
N = 10000000;               %number of pair connections
mean_value = 0.5;           %mean value of the gaussian distribution of connections
standard_deviation = 1/10;  %sqrt of variance. With 1/10 we are considering up to 5*sd. 3*sd is 99.7%

a = 0:0.1:0.9;              %pruning values
number_points = size(a,2);    

%% Plotting parameters
numericFontSize = 25;
axesFontSize = 30;
lineThickness = 2;
markLine = 1;
markSize = 12;

minv      = [ 0,    0,   -1,    0];
maxv      = [ 1,    2,    1,    1];
n_bins   = 200;      
bin_dim  = (maxv-minv) / n_bins;

%% Parameters and variables for symbolic calculation
syms u v w          %symbolic variables

mu = [mean_value,            2*mean_value,            0,                       0                     ];       %values used to create the functions: single, sum, difference, absolute difference
va = [standard_deviation^2,  2*standard_deviation^2,  2*standard_deviation^2,  2*standard_deviation^2];       %values used to create the functions: single, sum, difference, absolute difference

mu_sym = zeros(1, size(mu,2));        %mean value obtained by symbolic integration
va_sym = zeros(1, size(va,2));        %variance obtained by symbolic integration
mu_bid_sym = zeros(1,2);              %mean value obtained by symbolic integration
va_bid_sym = zeros(1,2);              %variance obtained by symbolic integration

tol=0.99;    %tolerance for the integration of the distribution

%% Printing vector
variable = char('x, y', 'x+y', 'x-y', '|x-y|'); %creates a matrix where each string is spread over a row, every character being an element of the matrix
x_axis_name = char ('w','Z_2','diff','Z_1');
y_axis_name = char ('f(w)','f(Z_2)','f(diff)','f(Z_1)');
figure_name = char ('Gaussian_initialPDF_prune0_1','Gaussian_sumPDF_prune0_1','Gaussian_diffPDF_prune0_1','Gaussian_absdiffPDF_prune0_1');

%% Sampling
for n = 1:number_points
    
    %% Simulation samples
    x_unp = normrnd(mean_value, standard_deviation, 1, N) .* (rand(1,N) > a(n));
    y_unp = normrnd(mean_value, standard_deviation, 1, N) .* (rand(1,N) > a(n));
    
    x = x_unp((x_unp+y_unp)~=0);    %uses logical index to cut away the zero elements from the distributions when they have the same index in x and y (w_{ij}=w_{ji}=0)
    y = y_unp((x_unp+y_unp)~=0);    %uses logical index to cut away the zero elements from the distributions when they have the same index in x and y (w_{ij}=w_{ji}=0)
    
    distr       = [      x;       x+y;       x-y;       abs(x-y)];  %vector of distributions
    mean_val    = [mean(x), mean(x+y), mean(x-y), mean(abs(x-y))];  %corresponding mean value
    variance    = [ var(x),  var(x+y),  var(x-y),  var(abs(x-y))];  %corresponding variance
    
    index_zero_x = find(x==0);
    index_zero_y = find(y==0);
    
    x_peak = [x(index_zero_x), x(index_zero_y)];
    y_peak = [y(index_zero_x), y(index_zero_y)];
    mean_peak = mean(x_peak+y_peak);
    var_peak = var(x_peak+y_peak);
 
    %% Single distributions
    for i = 1:size(distr,1)
        
        sample = minv(i):bin_dim(i):maxv(i);

        %---SYMBOLIC
        f1 = exp(-((u-mu(1))^2/(2*va(1)))) / sqrt(2*pi*va(1));
        f2 = exp(-((u-mu(i))^2/(2*va(i)))) / sqrt(2*pi*va(i));
        
        %Single distribution area check
        area_f2 = double(int(f2,minv(i),maxv(i)));
        
        if area_f2 > 1
            display('Error, something is going wrong: the integral of a distribution cannot be larger than 1')
        elseif area_f2 < tol
            f2 = f2 / area_f2;             
        end
        
        %overall distribution
        f = (2*a(n)/(1+a(n)))*f1 + ((1-a(n))/(1+a(n)))*f2;
        area = int(f,minv(i),maxv(i));
    
        %area check, normalization and sampling
        if double(area) > 1
            display('Error, something is going wrong: the integral of a distribution cannot be larger than 1')
        elseif double(area) < tol
            f = f / double(area);             
            val = (((2*a(n)/(1+a(n)))*normpdf(sample,mu(1),sqrt(va(1))))+(((1-a(n))/(1+a(n)))*normpdf(sample,mu(i),sqrt(va(i)))/area_f2)) / double(area);
            area = int(f,minv(i),maxv(i));
        else
            val = (((2*a(n)/(1+a(n)))*normpdf(sample,mu(1),sqrt(va(1))))+(((1-a(n))/(1+a(n)))*normpdf(sample,mu(i),sqrt(va(i)))/area_f2));
        end    
    
        %mean and variance by integration - for the absolute difference they have to be different from the ones used to create the function
        mu_sym(i) = int(f*u,minv(i),maxv(i));
        va_sym(i) = int(f*(u-mu_sym(i))^2,minv(i),maxv(i));
    
           
        %---SIMULATION
        [counts,bin_location] = hist(distr(i,:),n_bins);        %Note: bin_location==sample
        bin_dim_hist = ( max(bin_location) - min(bin_location) ) / n_bins;

    
        %print and plot
        sprintf('Pruning a = %f. \nVariable: %s \nSimulation with N=%d pair connections: Mean %f. Variance %f. \nSymbolic calculation: Integral of the distribution with mean %f and variance %f: %f', a(n), variable(i,:), N, mean_val(i), variance(i), mu_sym(i), va_sym(i), double(area))
            
        if n==1 || n==2 || n==6
            figure;
            bar(bin_location, counts/(sum(bin_dim_hist*counts)));     %simulation --- Note: sum(bin_dim(i)*counts) is no longer equal to bin_dim*N because of the pruning
            hold on
            if i == 1 
                h = plot(sample,(1-a(n))*val);     %theoretical
            else 
                h = plot(sample,val);              %theoretical
            end
            
            set(h, 'color', 'k', 'LineWidth', lineThickness);
            set(gca,'fontsize',numericFontSize);
            xlabel(x_axis_name(i, find(x_axis_name(i,:) ~= ' ')),'fontsize',axesFontSize);
            ylabel(y_axis_name(i, find(y_axis_name(i,:) ~= ' ')),'fontsize',axesFontSize);
            axis([minv(i)-0.5 maxv(i)+0.5 0 ceil(max( counts/(sum(bin_dim_hist*counts)) ))]);
                               
            colormap(gray);
            caxis([-1 1.4]);
        end
        
        if  (n == 2) && (i ~= 3)
            print(gcf, '-depsc2', '-loose', figure_name(i, find(figure_name(i,:) ~= ' '))); % Print the figure in eps (first option) and uncropped (second object)    
        end
        
    end


    %% Joint distributions-bivariate gaussian
    samplex = minv(2):bin_dim(2):maxv(2);   %first distribution is the sum
    sampley = minv(4):bin_dim(4):maxv(4);   %second distribution is or the difference or the absolute difference
    
    if n==1
        
        %covariance and correlation from the simulation
        covariance = sum((distr(2,:)-mean_val(2)).*(distr(4,:)-mean_val(4)))/N;
        corr = covariance/(sqrt(variance(2))*sqrt(variance(4)));
        sigma = [va(2), corr*sqrt(va(2))*sqrt(va(4)); corr*sqrt(va(2))*sqrt(va(4)), va(4)];
    
        %---SYMBOLIC
        g_un = exp( - (((v-mu(2))^2/(va(2))) + ((w-mu(4))^2/(va(4))) - 2*corr*(v-mu(2))*(w-mu(4))/(sqrt(va(2))*sqrt(va(4))) ) / (2*(1-(corr^2))) )  / ( 2*pi*sqrt (va(2)*va(4)*(1-(corr^2))) );
        %g_un = exp(-((v-mu(2))^2/(2*va(2)))-((w-mu(4))^2/(2*va(4))))/(2*pi*sqrt(va(2)*va(4)));
        area_bid = int((int(g_un,v,minv(2),maxv(2))),w,minv(4),maxv(4));
        
        %area check, normalization and sampling
        if double(area_bid)>1
            display('Error, something is going wrong: the integral of a distribution cannot be larger than 1')
        elseif double(area_bid)<tol
            g_un = g_un/double(area_bid);
            [X1,X2] = meshgrid(samplex,sampley);                     
            val = mvnpdf([X1(:) X2(:)], [mu(2),mu(4)], sigma)/double(area_bid);      %Note: mvnpdf returns a vector evaluating the PDF at each row of [X1(:) X2(:)]
            val = reshape(val,length(samplex),length(sampley));
            %area = int((int(g_un,v,minv(2),maxv(2))),w,minv(4),maxv(4));
        else
            [X1,X2] = meshgrid(samplex,sampley);                     
            val = mvnpdf([X1(:) X2(:)], [mu(2),mu(4)], sigma);      %Note: mvnpdf returns a vector evaluating the PDF at each row of [X1(:) X2(:)]
            val = reshape(val,length(samplex),length(sampley));
        end
        
        %mean and variance by integration - for the absolute difference they have to be different from the ones used to create the function
        mu_bid_sym(1,1) = int((int(g_un*v,v,minv(2),maxv(2))),w,minv(4),maxv(4));
        mu_bid_sym(1,2) = int((int(g_un*w,v,minv(2),maxv(2))),w,minv(4),maxv(4));
        va_bid_sym(1,1) = int((int(g_un*(v-mu_bid_sym(1,1))^2,v,minv(2),maxv(2))),w,minv(4),maxv(4));
        va_bid_sym(1,2) = int((int(g_un*(w-mu_bid_sym(1,2))^2,v,minv(2),maxv(2))),w,minv(4),maxv(4));
          
    else
        
        g2 = 0;
        g = ((1-a(n))/(1+a(n)))*g_un + (2*a(n)/(1+a(n)))*g2;            
        [X1,X2] = meshgrid(samplex,sampley);                     
        val = ((1-a(n))/(1+a(n)))*mvnpdf([X1(:) X2(:)], [mu(2),mu(4)], sigma)/double(area_bid);      %Note: mvnpdf returns a vector evaluating the PDF at each row of [X1(:) X2(:)]
        val = reshape(val,length(samplex),length(sampley));        
        %area = int((int(g,v,minv(2),maxv(2))),w,minv(4),maxv(4));
    end
    
     
    %---SIMULATION
    grid{1} = samplex;
    grid{2} = sampley;
    %[bidim_counts,Cell_location] = hist3([(distr(2,:))', (distr(i,:))'], grid);
    [bidim_counts,Cell_location] = hist3([(distr(2,:))', (distr(i,:))'], [n_bins, n_bins]);
    bin_dim_hist_1 = ( max(Cell_location{1}) - min(Cell_location{1}) ) / n_bins;
    bin_dim_hist_2 = ( max(Cell_location{2}) - min(Cell_location{2}) ) / n_bins;

%     %---PEAK DISTRIBUTION (without pruning)
%     %simulation
%     [bidim_counts_peak,Cell_location_peak] = hist3([(x_peak+y_peak)', (abs(x_peak-y_peak))'], grid);
%     surf(Cell_location_peak{1},Cell_location_peak{2},(bidim_counts_peak/max(max(bidim_counts_peak)))'/(2*pi*va(2)));    %(2*pi*va(2)) is the normalization from theory
%     hold on
%     
%     %theory    
%     mu_peak = [0.5 0.5];
%     va_peak = [(va(2)), 0; 0, (va(4))];   %covariance matrix
%     X_peak = mvnrnd(mu_peak,va_peak,10000000);
%     [counts_peak,Cell_peak] = hist3(X_peak, grid);
%     mesh(Cell_peak{1},Cell_peak{2},(counts_peak)'/(sum(sum(bin_dim(2)*bin_dim(4)*counts_peak))));
%     axis([0 1 0 1 0 10]);
    
    %print and plot
        sprintf('Pruning a = %f. \nVariables: %s, %s \nSimulation with N=%d pair connections: Mean %f, %f. Variance %f, %f. \nSymbolic calculation: Integral of the distribution with mean %f, %f and variance %f, %f: %f', a(n), variable(2,:), variable(4,:), N, mean_val(2), mean_val(4), variance(2), variance(4), mu_bid_sym(1,1), mu_bid_sym(1,2), va_bid_sym(1,1), va_bid_sym(1,2), double(area))
    
    if n==1 || n==2 || n==6
        figure;
        mesh(Cell_location{1},Cell_location{2},bidim_counts'/(sum(sum(bin_dim_hist_1*bin_dim_hist_2*bidim_counts))))      %simulation
        hold on
        surf(samplex,sampley,val);                                                                                %theoretical
        
        set(gca,'fontsize',numericFontSize);
        axis([0 2 0 1]);
        %xlabel('Z_2','fontsize',axesFontSize);
        %ylabel('Z_1','fontsize',axesFontSize);
        %zlabel('f(Z_1,Z_2)','fontsize',axesFontSize);
        view([130 14]);
        colormap(gray);
        caxis([0 140]);
    end
        
    if  (n == 2) 
        print(gcf, '-depsc2', '-loose', 'Gaussian_jointPDF_prune0_1'); % Print the figure in eps (first option) and uncropped (second object)    
    end
end


%%%%%-----Checking

%Simulation. Verify that all the points of the distribution lie in the range [0,1]x[abs(x-y),2-abs(x-y)] --- Few points can lie outside, those that in the formers distributions x and y are outside [0,1]
k=0;
m=0;
for i=1:size(distr,2)
    if and(((distr(2,i))>=(distr(4,i))),((distr(2,i))<=2-(distr(4,i))))
        k=k+1;            
    else
        m=m+1;
        %display(i)
    end
end

%Symbolic. Verify that the integral in the range [0,1]x[abs(x-y),2-abs(x-y)] is equal to 1
check = int((int(g_un,v,w,(2-w))),w,0,1);

%Print
sprintf('Checking the joint distribution f(x+y,|x-y|) in the definition interval [0,1]x[|x-y|,2-|x-y|]. \nSimulation. Number of points of the distribution outside: %d out of %f. \nSymbolic calculation. Area of the joint distribution: %f', m, N, double(check))

