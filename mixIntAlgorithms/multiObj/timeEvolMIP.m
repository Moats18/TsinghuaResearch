function timeEvolMIP(uopt, probAmp)
% plot the time evolution of the varying eigenstates based on an optimized 
% pulse sequence (uopt) and a starting probability amplitude (probAmp)
% inputted as a row vector

Dim = length(probAmp);
probAmp = probAmp';
TimeUf{1} = eye(Dim);
NN = length(uopt);
Delta = 0.01;
i = 1;

U = calculateU(5, 0.2, Dim, 0.01);

figure()
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
t = 0:Delta:((len-1)*Delta);

for i = 1:Dim
plot(t,propability(i,:));
hold on
end

xlabel('Time in nanoseconds');
ylabel('Probability Amplitude');
legend({'Ground State','Excited 1', 'Excited 2'}, 'Location', 'northwest');

end

