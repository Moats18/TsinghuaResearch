function IF = Qfit_timingError(ut,U,G, sigma, type)

% calculates the error from an optimized pulse sequence ut when a random 
% timing error is applied to each pulse. The timing error is randomly
% chosen from a normal probability distribution determined by sigmal. 
% A type of clock need to be specified: external clock where the 
% error distribution is constant or an internal clock where the error 
% distributed is varied over time. 


NN = length(ut);
N = size(U{1},1);
N0 = size(G,1);
Uf = eye(N);

if type == "External"
    n = 1;
for k = 1:NN
    TE = timingError(sigma, n);  
    [Udt1, Udt2] = FreeEvolutionTimeError(TE);
    switch ut(k)
        case 1
        Uf = U{1}*Uf;
        case 2
        Uf = Udt2*U{2}*Udt1*Uf; 
    end
end
end

if type == "Internal"
    n = 1;
    for k = 1:NN
        TE = timingError(sigma, n);  
        [Udt1, Udt2] = FreeEvolutionTimeError(TE);
        switch ut(k)
            case 1
            Uf = U{1}*Uf;
            case 2
            Uf = Udt2*U{2}*Udt1*Uf; 
            n = n + 1;
        end
    end
end

IF = 1-abs(trace(G'*Uf(1:N0,1:N0)))^2/4;

end






