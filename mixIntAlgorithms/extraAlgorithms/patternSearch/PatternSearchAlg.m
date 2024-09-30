function [minVal, uopt] = PatternSearchAlg(l, NN, x0)
%%%%%%%%%%%%%%% Calculation of U0 and U1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
w0 = 5;     %% qubit frequencyu (GHz);
chi = 0.2;  %% anharmonicity (GHz)
N = 4;      %% qubit cut-off dimension (in consideration of leakage)
Delta = 0.01;  %% clock period (ns) - each period includes one or no SFQ pulse  

A = zeros(N,N);

for k=1:N-1
    A(k,k+1) = sqrt(k);  %%% anhilation operator
end

H0 = w0*A'*A - chi/2*A'*A'*A*A;  %%% nonlinear oscillator
H1 = (A'+A)/2;                   %%% control Hamiltonian  

%%% Pulse-off  propagator U{1}
U{1} = expm(-1i*H0*Delta); % matrix exponential (if the matrix is diagonal then replace the diagonal elements with the exp(element) ) 
                           % matrix H0 is diagonal 

%%% Pulse-on  propagator U{2}
mu = 0.4*Delta;   % SFQ pulse center 
sigma= 0.1*Delta; % SFQ pulse width/standard deviation
dtheta = pi/20;  % SFQ pulse area

M = 100;
dt = Delta/M;
t = (1:M)*dt; %creates an array of times ranging from dt to delta with time step dt 

% Gaussian shape with a curve area of dtheta instead of 1
% The pulse is split into smaller chunks by the interval dt

sfq = dtheta*exp(-(t-mu).^2/2/sigma^2)/sqrt(2*pi)/sigma ;  

U{2} = eye(N); %matrix with N rows and cols with 1 along the main diagonal

% each loop mulitplies the current time evolution operator by the previous one 
% each time evolution operator is a chunk of the total time evolution
% caused by the discrete nature of the pulse 

for k = 1:length(t)
    U{2} = expm(-1i*(H0+sfq(k)*H1)*dt)*U{2}; 
end
%%%%%%%%%%%%%%%%%% Calculation of U0 and U1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%% Mixed integer programming using Pattern Search Algorithm %%%%%%%%%%%%%%%
%maxPulses = 1000; % maximum total number of clock periods
%constraint = @(x) constraintFunction(maxPulses, x); 
G = [0 1;-1 0]; %% Target gate
option =  optimoptions('patternsearch','PlotFcn','psplotbestf', "ScaleMesh", "Off", "TolMesh", 0.9);

myfit = @(ut) Qfit(ut,U,G); % sets the function myfit to be a single variable function with constants U and G 
[uopt, minVal] = patternsearch(myfit,x0, [],[],[],[],-l*ones(1,NN),l*ones(1,NN),[], option);

end

