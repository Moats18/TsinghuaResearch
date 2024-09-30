% Function that evaluates the fitness of a certain sequence i.e. the error
% that the sequence has with not gate transformation
% Splits the pulse sequence up into the multiplication of
% indivdual time evloution operateros which correspond to a pulse being
% applied U{2} or not applied U{1} 
% rounds values because the Particle Swarm sweeps continuous values

function [IF, Uf] = QfitPartSwarm(ut,U,G)

NN = length(ut);
N = size(U{1},1);
N0 = size(G,1);
Uf = eye(N);
ut = round(ut);

for k = 1:NN
    if ut(k) <= 0
        Uf = (U{1}^abs(ut(k)))*Uf;
    else 
        Uf = (U{2}^ut(k))*Uf; 
    end
end
%IF = 1-abs(trace(G'*Uf(1:N0,1:N0)))^2/4;
IF = immse(G, Uf(1:N0,1:N0));

end