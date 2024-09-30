function [minVal, uopt] = GAcrossoverFcn(crossoverFunction, varargin)
% Used to determine which crossover function produced the lowest gate errors
%
% Example command line usage: GAcrossoverFcn(crossoverFunction, L, NN, x0)
% L -  length of the segments that make up the pulse sequence
% NN - number of segments that make up the pulse sequence 
% x0 - initial point for the algorithm to start. Used in the case where 
% an optimization point from a different algorithm is validated

%Default Inputs
if nargin == 1
    L = 100;
    NN = 20;
    x0 = randi([-L L], 1, NN);
% Assign Inputs
elseif nargin == 3
    L = cell2mat(varargin(1));
    NN = cell2mat(varargin(2));
    x0 = randi([-L L], 1, NN);
elseif nargin == 4
    L = cell2mat(varargin(1));
    NN = cell2mat(varargin(2));
    x0 = cell2mat(varargin(3));
end

%%%%%%%%%%%%%%% Calculation of U0 and U1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Default Variables 
w0 = 5;     %% qubit frequencyu (GHz);
chi = 0.2;  %% anharmonicity (GHz)
N = 7;      %% qubit cut-off dimension (in consideration of leakage)
Delta = 0.01;  %% clock period (ns) - each period includes one or no SFQ pulse  

U = calculateU(w0, chi, N, Delta);
%%%%%%%%%%%%%%% Calculation of U0 and U1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%% Mixed integer programming using Genetic Algorithm %%%%%%%%%%%%%%%
option =  optimoptions('ga','InitialPopulationMatrix', x0,'CrossoverFcn' , crossoverFunction, ...
    'PlotFcn', 'gaplotbestf', 'crossoverFraction', 0.9, 'EliteCount', 25, 'FitnessScalingFcn', ...
    'fitscalingshiftlinear', 'PopulationSize', 1000);

maxPulses = 2000; % maximum total number of clock periods
constraint = @(x) constraintFunction(maxPulses, x); % total number is less than or equal to maxPulses
G = [0 1;-1 0]; %% Target gate
myfit = @(ut) Qfit(ut,U,G); % sets the function myfit to be a single variable function with constants U and G 

[uopt, minVal] = ga(myfit,NN,[],[],[],[],-L*ones(1,NN),L*ones(1,NN),constraint,1:NN, option);
%disp(['Minimized Gate Error = ',num2str(J0)]);

