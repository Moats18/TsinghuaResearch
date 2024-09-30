function [Ud1, Ud2] = FreeEvolutionTimeError(TE)
% Calculates the forward and backward time evolution of the quantum state
% used to simulate a scenario where a "jitter"/time error - variable TE -
% occurs in the application of a SFQ pulse 

%Default Variables 
w0 = 5;     %% qubit frequencyu (GHz);
chi = 0.2;  %% anharmonicity (GHz)
N = 3;      % N qubit cut-off dimension (in consideration of leakage)

A = zeros(N,N);
for k=1:N-1
    A(k,k+1) = sqrt(k);  %%% anhilation operator
end

H0 = w0*A'*A - chi/2*A'*A'*A*A;  %%% nonlinear oscillator 

Ud1 = expm(-1i*H0*TE); 
Ud2 = expm(-1i*H0*(-TE));

end

