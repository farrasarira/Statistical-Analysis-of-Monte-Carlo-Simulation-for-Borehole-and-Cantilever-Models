clc
clear all
close all
%% Monte Carlo Simulation

N = 100000;         % Number of Simulations

% ------- INPUT ----------------------------
R = zeros(N,1);     % yield stress
E = zeros(N,1);     % Young's modulus of beam material
X = zeros(N,1);     % horizontal load
Y = zeros(N,1);     % vertical load

D0 = 2.2535;

% ------- OUTPUT ---------------------------
D = zeros(N,1);     % displacement
S = zeros(N,1);     % stress 

i = 1;

while i <= N
    
    R(i) = normrnd(40000,2000);
    E(i) = normrnd(2.9E7,1.45E6);
    X(i) = normrnd(500,100);
    Y(i) = normrnd(1000,100);
    
    [D(i), S(i)] = canti([R(i) E(i) X(i) Y(i)]);
    D(i) = D(i) / D0 - 1;   % Scaling the Displacement
    S(i) = S(i) / R(i) - 1; % Scaling the Stress
    
    i = i + 1;
end

%% Plot Histogram

figure(1);
binWidth = 0.05;
lastVal = ceil(max(D));
firstVal = floor(min(D));
binEdges = firstVal-1:binWidth:lastVal+1;
h = histogram(D,binEdges);
xlabel('Displacement');
ylabel('Frequency');

figure(2);
binWidth = 0.05;
lastVal = ceil(max(S));   
firstVal = floor(min(S));
binEdges = firstVal-1:binWidth:lastVal+1;
h = histogram(S,binEdges);
xlabel('Stress');
ylabel('Frequency');

%% Plot Probability Distribution Function

pd_D= fitdist(D,"Normal");    % Curve Fitting
pd_S= fitdist(S,"Normal");    % Curve Fitting

figure(3);
binWidth = 0.05;
lastVal = ceil(max(D));
firstVal = floor(min(D));
binEdges = firstVal-1:binWidth:lastVal+1;
h = histogram(D,binEdges,'Normalization','pdf');
xlabel('Displacement');
ylabel('Probability Density');
hold on;
xgrid = [firstVal:0.01:lastVal];
pdfEst = pdf(pd_D,xgrid);
plot(xgrid,pdfEst,'r-','LineWidth',1.5);
legend('actual probability density','approximated probability density function');
hold off;

figure(4);
binWidth = 0.05;
lastVal = ceil(max(S));
firstVal = floor(min(S));
binEdges = firstVal-1:binWidth:lastVal+1;
h = histogram(S,binEdges,'Normalization','pdf');
xlabel('Stress');
ylabel('Probability Density');
hold on;
xgrid = [firstVal:0.01:lastVal];
pdfEst = pdf(pd_S,xgrid);
plot(xgrid,pdfEst,'r-','LineWidth',1.5);
legend('actual probability density','approximated probability density function');
hold off;

%% Calculating Statistical Properties of PDF
%----- Displacement -------------------------------
mean_D = mean(D);           % mean
fluct_D = D - mean_D;       % fluctuation
var_D = mean(fluct_D.^2);   % variance 
std_D = sqrt(var_D);        % standard deviation
skewness_D = mean(fluct_D.^3)/std_D^3;  % skewness
kurtosis_D = mean(fluct_D.^4)/std_D^4;  % kurtosis

%----- Stress -------------------------------------
mean_S = mean(S);           % mean
fluct_S = S - mean_S;       % fluctuation
var_S = mean(fluct_S.^2);   % variance 
std_S = sqrt(var_S);        % standard deviation
skewness_S = mean(fluct_S.^3)/std_S^3;  % skewness
kurtosis_S = mean(fluct_S.^4)/std_S^4;  % kurtosis




