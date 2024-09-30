
% Default Values
w0 = 5;     %% qubit frequencyu (GHz);
chi = 0.2;  %% anharmonicity (GHz)
N = 3; % qubit cut-off dimension (in consideration of leakage)
Delta = 0.01;  %% clock period (ns) - each period includes one or no SFQ pulse 
G = [0 1;-1 0]; %% Target gate


U = calculateU(w0, chi, N, Delta);

% Specify these values
sigma = 0.001; % standard deviation of timing error 
type = "Internal"; % type of clock 
error = Qfit_timingError(uopt, U, G, sigma, type);
