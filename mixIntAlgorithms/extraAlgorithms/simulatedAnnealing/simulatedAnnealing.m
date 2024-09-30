function [minVal, uopt] = simulatedAnnealing(varargin)
% Example command line usage: simulatedAnnealing(L, NN, x0)
% L -  length of the segments that make up the pulse sequence
% NN - number of segments that make up the pulse sequence 
% x0 - initial point for the algorithm to start. Typically used in the case 
% where an optimization point from a different algorithm is validated 

%Default Inputs
if nargin == 0
    L = 100;
    NN = 20;
    x0 = randi([-L L], 1, NN);
% Assign Inputs
elseif nargin == 2
    L = cell2mat(varargin(1));
    NN = cell2mat(varargin(2));
    x0 = randi([-L L], 1, NN);
elseif nargin == 3
    L = cell2mat(varargin(1));
    NN = cell2mat(varargin(2));
    x0 = cell2mat(varargin(3));
end
%%%%%%%%%%%%%%% Calculation of U0 and U1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Default Variables 
w0 = 5;     %% qubit frequency (GHz);
chi = 0.2;  %% anharmonicity (GHz)
lb = 4;      % N qubit cut-off dimension (in consideration of leakage)
ub = 4;
Delta = 0.01;  %% clock period (ns) - each period includes one or no SFQ pulse  


for N = lb:ub
    U{N+1-lb} = calculateU(w0, chi, N, Delta);
end

%%%%%%%%%%%%%%%%%% Calculation of U0 and U1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%% Mixed integer programming using Simulated Annealing Algorithm %%%%%%%%%%%%%%%
options = optimoptions(@simulannealbnd, 'PlotFcn', @saplotbestf,'InitialTemperature', 1);
G = [0 1;-1 0]; %% Target gate
myfit = @(ut) QfitAvg(ut,U,G); % sets the function myfit to be a single variable function with constants U and G 

[uopt, minVal] = simulannealbnd(myfit,x0, -L*ones(1,NN),L*ones(1,NN), options);
uopt = round(uopt);

end



