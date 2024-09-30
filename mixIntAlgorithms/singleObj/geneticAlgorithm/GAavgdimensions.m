function [minVal, uopt] = GAavgdimensions(varargin)
% Example command line usage: GeneticAlgorithm(L, NN, x0)
% L -  length of the segments that make up the pulse sequence
% NN - number of segments that make up the pulse sequence 
% x0 - initial point for the algorithm to start. Used in the case where 
% an optimization point from a different algorithm is validated 

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

%%%%%%%%%%%%%% Mixed integer programming using Genetic Algorithm %%%%%%%%%%%%%%%
crossoverFunction = @crossoverintermediate;
option =  optimoptions('ga','InitialPopulationMatrix', x0,'CrossoverFcn' , crossoverFunction, ...
    'PlotFcn', 'gaplotbestf', 'crossoverFraction', 0.9, 'EliteCount', 25, 'FitnessScalingFcn', ...
    'fitscalingshiftlinear', 'PopulationSize', 1000, "MaxGenerations", 1000);

G = [0 1;-1 0]; %% Target gate
myfit = @(ut) QfitAvg(ut,U,G); % sets the function myfit to be a single variable function with constants U and G 

[uopt, minVal] = ga(myfit,NN,[],[],[],[],-L*ones(1,NN),L*ones(1,NN),[],1:NN, option);
%disp(['Average Minimized Gate Error = ',num2str(minVal)]);

end

