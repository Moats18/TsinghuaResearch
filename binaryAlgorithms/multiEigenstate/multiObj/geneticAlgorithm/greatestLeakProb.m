function totalLeakage = greatestLeakProb(uopt, U)
% determines whether the leakage probabilities given an optimized pulse sequence
% are bounded by the probabilities representing the definite ground state
% or the definite first eigenstate. And returns the average total leakage -
% calculated by a summation of the probabilites at each time step divided by 
% number of time steps. 

sum = zeros(2,1);
Dim = length(U{1});
len = length(uopt);

for n = 1:2
    if n == 1
    probAmp = [1, zeros(1, Dim-1)]';
    else
    probAmp = [0, 1, zeros(1, Dim-2)]';
    end

TimeUf{1} = eye(Dim);
propability = zeros(Dim, 1, len+1);

for k = 1:len
    switch uopt(k)
        case 1
        TimeUf{k+1} = U{1}*TimeUf{k};
        case 2
        TimeUf{k+1} = U{2}*TimeUf{k}; 
    end
end

for l = 1:len+1
propability(:, 1, l) = TimeUf{l}*probAmp;
end

propability = abs(propability);

% summation of the leakage probablitity for higher order eigenstates 
% (n >= 3)

for i = Dim-2:Dim
 for z = 1:len
 sum(n) = sum(n) + propability(i, :, z);
 end
end 

sum(n) = sum(n)/((Dim-2)*len);
end 

% determining which probability produces the highest total leakage
if sum(1) > sum(2)
    totalLeakage = sum(1);
else 
    totalLeakage = sum(2);
end

end

