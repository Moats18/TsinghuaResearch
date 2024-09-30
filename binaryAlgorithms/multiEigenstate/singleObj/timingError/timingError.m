function timeError = timingError(sigma, n)
% generates a random timing error for an applied pulse
% from a normally distributed random time 
% described by standard dev sigma and expected value mu
% External clock n = 1;
% Internal clock n = applied pulse number

mu = 0; % expected timing error 
timeError = normrnd(mu,sqrt(n)*sigma);

end

