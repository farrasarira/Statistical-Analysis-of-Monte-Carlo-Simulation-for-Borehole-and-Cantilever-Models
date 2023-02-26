clc
clear all
close all
%% Monte Carlo Simulation

N = 100000;    % Number of Simulations

%------- Input -------------------------------------------------
rw = zeros(N,1);    % radius of borehole [m]
r = zeros(N,1);     % radius of influence [m]
Tu = zeros(N,1);    % transmissivity of upper aquifer [m2/yr] (uniform distribution)
Hu = zeros(N,1);    % potentiometric head of upper aquifer [m] (uniform distribution)
Tl = zeros(N,1);    % transmissivity of lower aquifer [m2/yr] (uniform distribution)
Hl = zeros(N,1);    % potentiometric head of lower aquifer [m] (uniform distribution)
L = zeros(N,1);     % length of borehole [m] (uniform distribution)
Kw = zeros(N,1);    % hydraulic conductivity of borehole [m/yr] (uniform distribution)

%------- Output ------------------------------------------------
f = zeros(N,1); % Water Flow Rate [m^3/yr] (Random Variables)

i = 1;  

while i<=N
    rw(i) = normrnd(0.10,0.0161812);  
    r(i) = lognrnd(7.71,1.0056);     
    Tu(i) = 63070 + 52530*rand;       
    Hu(i) = 990 + 120*rand;            
    Tl(i) = 63.1 + 52.9*rand;         
    Hl(i) = 700 + 120*rand;            
    L(i) = 1120 + 560*rand;            
    Kw(i) = 9855 + 2190*rand;          

    if rw(i) >= 0.05 && rw(i) <= 0.15 && r(i) >= 100 && r(i) <= 50000       % check whather the random sampling fulfill the constraints
        f(i) = borehole([rw(i), r(i), Tu(i), Hu(i), Tl(i), Hl(i), L(i), Kw(i)]);
        i = i + 1;
    end
end


%% Plot Histogram

figure(1);
binWidth = 5;
lastVal = ceil(max(f));
binEdges = 0:binWidth:lastVal+1;
h = histogram(f,binEdges);
xlabel('Water Flow Rate');
ylabel('Frequency');

%% Plot Probability Distribution Function

pd = fitdist(f,"Gamma");    % Curve Fitting

figure(2);
h = histogram(f,binEdges,"Normalization","pdf");
xlabel('Water Flow Rate []');
ylabel('Probability Density');
hold on;
xgrid = [0:300];
pdfEst = pdf(pd,xgrid);
plot(xgrid,pdfEst,'r-','LineWidth',1.5);
legend('Simulated Probability Distribution','approximated Probability Distribution Function');

%% Calculating Statistical Properties of PDF

mean_f = mean(f);           % mean
fluct_f = f - mean_f;       % fluctuation
var_f = mean(fluct_f.^2);   % variance 
std_f = sqrt(var_f);        % standard deviation
skewness_f = mean(fluct_f.^3)/std_f^3;  % skewness
kurtosis_f = mean(fluct_f.^4)/std_f^4;  % kurtosis

