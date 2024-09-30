function [paretoPoints, fVal] = GAbinary(varargin)
% Example command line usage: GeneticAlgorithm(L, NN, x0)
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
%Default Variables 
w0 = 5;     %% qubit frequencyu (GHz);
chi = 0.2;  %% anharmonicity (GHz)
lb = 4;      %% qubit cut-off dimension (in consideration of leakage)
ub = 4;
Delta = 0.01;  %% clock period (ns) - each period includes one or no SFQ pulse  


for N = lb:ub
    U{N+1-lb} = calculateU(w0, chi, N, Delta);
end
%%%%%%%%%%%%%%% Calculation of U0 and U1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%% Binary integer programming using Multi-Objective Genetic Algorithm %%%%%%%%%%%%%%%
option =  optimoptions('gamultiobj','InitialPopulationMatrix', x0,'PlotFcn', 'gaplotpareto', 'CrossoverFraction',0.9,...
    'PopulationSize', 100, 'ParetoFraction', 0.30);

G = [0 1;-1 0]; %% Target gate
myfit = @(ut) fitnessFunction(ut, U, G);

[paretoPoints, fVal] = gamultiobj(myfit,NN,[],[],[],[],ones(1,NN),2*ones(1,NN),[], 1:NN,option);

end

