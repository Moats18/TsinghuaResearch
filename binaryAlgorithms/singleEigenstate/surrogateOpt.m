NN = 1000; %%% number of clock periods
x0 = randi([1 2], 1, NN);

%%%%%%%%%%%%%%% Calculation of U0 and U1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
w0 = 5;     %% qubit frequencyu (GHz);
chi = 0.2;  %% anharmonicity (GHz)
N = 2;      %% qubit cut-off dimension (in consideration of leakage)
Delta = 0.01;  %% clock period (ns) - each period includes one or no SFQ pulse  

A = zeros(N,N);
for k=1:N-1
    A(k,k+1) = sqrt(k);  %%% anhilation operator
end
H0 = w0*A'*A - chi/2*A'*A'*A*A;  %%% nonlinear oscillator
H1 = (A'+A)/2;                   %%% control Hamiltonian  

%%% Pulse-off propagator U{1}
U{1} = expm(-1i*H0*Delta);
%%% Pulse-on  propagator U{2}
mu = 0.4*Delta;   % SFQ pulse width
sigma= 0.1*Delta; % SFQ pulse center
dtheta = pi/20;  % SFQ pulse area
M = 100;
dt = Delta/M;
t = (1:M)*dt;
sfq = dtheta*exp(-(t-mu).^2/2/sigma^2)/sqrt(2*pi)/sigma;  %% Gaussian shape
U{2} = eye(N);
for k = 1:length(t)
    U{2} = expm(-1i*(H0+sfq(k)*H1)*dt)*U{2};
end
%%%%%%%%%%%%%%%%%% Calculation of U0 and U1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%% Binary integer programming using Genetic Algorithm %%%%%%%%%%%%%%%
G = [0 1;-1 0]; %% Target gate
option =  optimoptions('surrogateopt', 'InitialPoints', x0, 'MaxFunctionEvaluations', 1000);

myfit = @(ut) Qfit_Original(ut,U,G);
[uopt, minVal] = surrogateopt(myfit, 1*ones(1,NN), 2*ones(1,NN), 1:NN, [], [], [], [], option);