function totalLeakage = greatestLeakProb(uopt, U)
% determines whether the leakage probabilities given an optimized pulse sequence
% are bounded by the probabilities representing the definite ground state
% or the definite first eigenstate. And returns the average total leakage -
% calculated by a summation of the probabilites at each time step divided by 
% number of time steps. 

sum = zeros(2,1);
Dim = length(U{1});

for n = 1:2
    if n == 1
    probAmp = [1, zeros(1, Dim-1)]';
    else
    probAmp = [0, 1, zeros(1, Dim-2)]';
    end

TimeUf{1} = eye(Dim);
NN = length(uopt);
i = 1;

for s = 1:NN
    if uopt(s) <= 0
        for k = 1:abs(uopt(s))
        TimeUf{i+1} = U{1}*TimeUf{i};
        i = i + 1;
        end

    else
        for k = 1:uopt(s)
        TimeUf{i+1} = U{2}*TimeUf{i};
        i = i + 1;
        end
    end
end

len = length(TimeUf);
propability = zeros(Dim, 1, len-1);

for t = 1:len
propability(:, 1, t) = TimeUf{t}*probAmp;
end

propability = abs(propability);

for z = 1:len
sum(n) = sum(n) + propability(3, :, z);
end

sum(n) = sum(n)/len;
end 

% determining which probability produces the highest total leakage
if sum(1) > sum(2)
    totalLeakage = sum(1);
else 
    totalLeakage = sum(2);
end

end

