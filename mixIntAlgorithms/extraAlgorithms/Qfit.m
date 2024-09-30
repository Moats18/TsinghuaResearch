function [IF, Uf] = Qfit(ut,U,G)
% Used for Mixed Integer Problem
% Function that evaluates the fitness of a certain sequence i.e. the error
% that the sequence has with the gate transformation
% Splits the pulse sequence up into the multiplication of
% indivdual time evolution operaters which correspond to a pulse being
% applied U{2} or not applied U{1} 


NN = length(ut);
N = size(U{1},1);
N0 = size(G,1);
Uf = eye(N);

for k = 1:NN
    if ut(k) <= 0
        Uf = (U{1}^abs(ut(k)))*Uf;
    else 
        Uf = (U{2}^ut(k))*Uf;
    end
end

IF = 1-abs(trace(G'*Uf(1:N0,1:N0)))^2/4;

end




