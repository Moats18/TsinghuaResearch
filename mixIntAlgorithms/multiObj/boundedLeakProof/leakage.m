function [best, avg, state0, state1] = leakage(uopt)
% from a given optimized pulse sequence (uopt) the average and highest total leakage from a 
% random starting probability of the first two eigenstates - determined through 
% a summation of the leakage probability at each time step - is calculated along 
% with the leakage probability of the probability where the state is definitely the ground
% state or the first excited state. 

% illustrates that for any uopt the highest total leakage occurs either when the
% state is definitely the ground state or the first excited state. 
 
N = 1000; % number of total random probabilities that are calcualted

sum = zeros(N, 1);

for n = 1:N
    if n <=N-2
    p = rand(1,1);
    probAmp = [p, 1-p, 0]';
    elseif n == N-1
    probAmp = [1, 0, 0]';
    else 
    probAmp = [0, 1, 0]';
    end

Dim = length(probAmp);
TimeUf{1} = eye(Dim);
NN = length(uopt);
i = 1;

U = calculateU(5, 0.2, Dim, 0.01);

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
end 

state0 = sum(N-1);
state1 = sum(N);

best = sum(1);
s = 0;
 
for j = 1:N
        if sum(j) > best
            best = sum(j);
        end
        s = s + sum(j);
 end 
 avg = s/N;

end

