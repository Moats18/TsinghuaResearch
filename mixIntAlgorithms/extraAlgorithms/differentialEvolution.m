function differentialEvolution(varargin)

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

% Set title
optimInfo.title = 'Using Differenial EvolutionAlgorithm';

% Specify objective function
objFctHandle = @(ut) QfitAvg(ut,U,G); 

% Define parameter names, ranges and quantization:

% 1. column: parameter names
% 2. column: parameter ranges
% 3. column: parameter quantizations
% 4. column: initial values (optional)

paramName = 'Optimized Pulse Sequence';
range = [-L L];

%sets the quantization to 1 indicating that the values are integers
quant = 1;

paramDefCell = {paramName, repmat(range, NN), ones(quant, NN)' , x0'};

% Set single additional function parameter
objFctSettings = 100;

% Get default DE parameters
DEParams = getdefaultparams;

% Set number of population members (often 10*D is suggested) 
DEParams.NP = 20;

% Do not use slave processes here. If you want to, set feedSlaveProc to 1 and
% run startmulticoreslave.m in at least one additional Matlab session.
DEParams.feedSlaveProc = 0;

% Set times
DEParams.maxiter  = 20;
DEParams.maxtime  = 30; % in seconds
DEParams.maxclock = [];

% Set display options
DEParams.infoIterations = 1;
DEParams.infoPeriod     = 10; % in seconds

% Do not send E-mails
emailParams = [];

% Start differential evolution
[bestmem, bestval, bestFctParams, nrOfIterations, resultFileName] = differentialevolution(...
	DEParams, paramDefCell, objFctHandle, objFctSettings, objFctParams, emailParams, optimInfo); %#ok

disp(' ');
disp('Best parameter set returned by function differentialevolution:');
disp(bestFctParams);

end

