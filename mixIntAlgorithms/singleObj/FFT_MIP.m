function Z = FFT_MIP(uopt)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

NN = length(uopt);
Delta = 0.01;  %% clock period (ns) - each period includes one or no SFQ pulse  

mu = 0.4*Delta;   % SFQ pulse width
sigma= 0.1*Delta; % SFQ pulse center
dtheta = pi/300;  % SFQ pulse area
M = 100;
dt = Delta/M;
t = (1:M)*dt;       % Time vector
sfq = dtheta/(sigma*sqrt(2*pi))*(exp(-(t-mu).^2/(2*sigma^2)));

%%%%%%%%%%%%%%%%%%
figure()
ut = [];
for i = 1:NN
    if uopt(i) <= 0 
        ut = [ut,zeros(1,M*abs(uopt(i)))];
    else
        ut = [ut,repmat(sfq, 1, uopt(i))];
    end
end

plot((1:length(ut))*dt, ut);

%L = length(ut);      % Signal length
L = 4000;

%%%%%%%%%%%%%%%
title("Gaussian Pulses in Time Domain")
xlabel("Time (t)")
ylabel("X(t)")

figure()
n = 2^nextpow2(L);
Z = fft(ut,n);
f = M*(0:(n/2))/n;
P = abs(Z/n).^2;

plot(f,P(1:n/2+1)) 
title("Gaussian Pulse in Frequency Domain")
xlabel("f (Hz)")
ylabel("|P(f)|^2")


end


