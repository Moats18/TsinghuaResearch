function [minVal, uopt] = SAavgdimensions(varargin)
% Example command line usage: SAavgdimensions(NN, x0)
% NN - number of time steps used the pulse sequence 
% x0 - initial point for the algorithm to start. Used in the case where 
% an optimization point from a different algorithm is validated 

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
w0 = 5;     %% qubit frequency (GHz);
chi = 0.2;  %% anharmonicity (GHz)
lb = 2;      % N qubit cut-off dimension (in consideration of leakage)
ub = 4;
Delta = 0.01;  %% clock period (ns) - each period includes one or no SFQ pulse  


for N = lb:ub
    U{N+1-lb} = calculateU(w0, chi, N, Delta);
end
%%%%%%%%%%%%%%%%%% Calculation of U0 and U1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%% Mixed integer programming using Surrogate Solver %%%%%%%%%%%%%%%%
G = [0 1;-1 0]; %% Target gate
option =  optimoptions('surrogateopt', 'InitialPoints', x0, 'MaxFunctionEvaluations', 10000);

myfit = @(ut) QfitAvg(ut,U,G); % sets the function myfit to be a single variable function with constants U and G  
[uopt, minVal] = surrogateopt(myfit, 1*ones(1,NN), 2*ones(1,NN), 1:NN, [], [], [], [], option);


