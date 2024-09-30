function timeEvolBin(uopt, probAmp)
% plot the time evolution of the varying eigenstates based on an optimized 
% pulse sequence (uopt) and a starting probability amplitude (probAmp)
% inputted as a row vector

probAmp = probAmp';
NN = length(probAmp); 
TimeUf{1} = eye(NN);
len = length(uopt);
propability = zeros(NN, 1, len+1);
U = calculateU(5, 0.2, NN, 0.01);

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
t = 0:0.01:(len*0.01);
figure()
for i = 1:NN
plot(t,propability(i,:));
hold on
end
xlabel('Time in nanoseconds');
ylabel('Probability Amplitude');
legend({'Ground State','Excited 1', 'Excited 2', 'Excited 3'}, 'Location', 'northwest');

end

