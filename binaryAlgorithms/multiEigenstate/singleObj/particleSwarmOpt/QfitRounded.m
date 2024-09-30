function IF = QfitRounded(ut,U,G)
% Function that evaluates the fitness of a certain sequence i.e. the error
% that the sequence has with the gate transformation G
% Splits the pulse sequence up into the multiplication of
% indivdual time evloution operateros which correspond to a pulse being
% applied U{2} or not applied U{1} 

% rounds values because Particle Swarm and Simulated Annealing sweep
% continuous values

NN = length(ut);
N = size(U{1},1);
N0 = size(G,1);
Uf = eye(N);
ut = round(ut);

for k = 1:NN
    switch ut(k)
        case 1
        Uf = U{1}*Uf;
        case 2
        Uf = U{2}*Uf; 
    end
end

IF = 1-abs(trace(G'*Uf(1:N0,1:N0)))^2/4;

end