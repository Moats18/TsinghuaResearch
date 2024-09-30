function [IF, Uf] = Qfit_Original(ut,U,G)

NN = length(ut);
N = size(U{1},1);
N0 = size(G,1);
Uf = eye(N);

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






