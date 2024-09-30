function U = squarePulse(w0, chi, N, Delta)

%SQUAREPULSE 

% w0 - qubit frequencyu (GHz);
% chi - anharmonicity (GHz)
% N - qubit cut-off dimension (in consideration of leakage)
% Delta - clock period (ns) - each period includes one or no SFQ pulse  

A = zeros(N,N);

for k=1:N-1
    A(k,k+1) = sqrt(k);  %%% anhilation operator
end

H0 = w0*A'*A - chi/2*A'*A'*A*A;  %%% nonlinear oscillator
H1 = (A'+A)/2;                   %%% control Hamiltonian  

%%% Pulse-off  propagator U{1}
U{1} = expm(-1i*H0*(19*Delta));

area = pi/300;
height = area/Delta;

U{2} = expm(-1i*(H0+(height*H1))*(Delta));

end

