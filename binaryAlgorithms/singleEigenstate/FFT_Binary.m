function FFT_Binary(uopt)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

figure()
NN = length(uopt);
Delta = 0.01;  %% clock period (ns) - each period includes one or no SFQ pulse  

mu = 0.4*Delta;   % SFQ pulse width
sigma= 0.1*Delta; % SFQ pulse center
dtheta = pi/20;  % SFQ pulse area
M = 100;
dt = Delta/M;
t = (1:M)*dt;       % Time vector
sfq = dtheta/(sigma*sqrt(2*pi))*(exp(-(t-mu).^2/(2*sigma^2)));

%%%%%%%%%%%%%%%%%%
ut = [];
M = length(sfq);
dt = Delta/M;
for k = 1:NN
    switch uopt(k)
        case 1
            ut = [ut,zeros(size(sfq))];
        case 2
            ut = [ut,sfq];
    end
end
plot((1:M*NN)*dt,ut);
L = M*NN;      % Signal length
%%%%%%%%%%%%%%%

title("Gaussian Pulses area " + num2str(dtheta) + " in Time Domain")
xlabel("Time (t)")
ylabel("X(t)")

figure()
n = 2^nextpow2(L);
Z = fft(ut,n);
f = M*(0:(n/2))/n;
P = abs(Z/n).^2;

plot(f,P(1:n/2+1)) 
title(" Gaussian Pulses of area " + num2str(dtheta) + " in Frequency Domain")
xlabel("f (Hz)")
ylabel("|P(f)|^2")


end

