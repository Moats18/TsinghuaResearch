function [minVal, uopt] = simulatedAnnealing(varargin)
% Example command line usage: simulatedAnnealing(NN, x0)
% NN - number of time steps that make up the pulse sequence 
% x0 - initial point for the algorithm to start. Typically used in the case 
% where an optimization point from a different algorithm is validated 

%Default Inputs
if nargin == 0
    NN = 1000;
    x0 = randi([1 2], 1, NN);
% Assign Inputs
elseif nargin == 1
    NN = cell2mat(varargin(1));
    x0 = randi([1 2], 1, NN);
elseif nargin == 2
    NN = cell2mat(varargin(1));
    x0 = cell2mat(varargin(2));
end

%%%%%%%%%%%%%%% Calculation of U0 and U1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Default Variables 
w0 = 5;     %% qubit frequency (GHz);
chi = 0.2;  %% anharmonicity (GHz)
lb = 3;      % N qubit cut-off dimension (in consideration of leakage)
ub = 3;
Delta = 0.01;  %% clock period (ns) - each period includes one or no SFQ pulse  


for N = lb:ub
    U{N+1-lb} = calculateU(w0, chi, N, Delta);
end

%%%%%%%%%%%%%%%%%% Calculation of U0 and U1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%% Mixed integer programming using Simulated Annealing Algorithm %%%%%%%%%%%%%%%
options = optimoptions(@simulannealbnd, 'PlotFcn', @saplotbestf,'InitialTemperature',100);
G = [0 1;-1 0]; %% Target gate
myfit = @(ut) QfitAvg(ut,U,G); % sets the function myfit to be a single variable function with constants U and G 

[uopt, minVal] = simulannealbnd(myfit,x0, ones(1,NN),2*ones(1,NN), options);
uopt = round(uopt);

end



