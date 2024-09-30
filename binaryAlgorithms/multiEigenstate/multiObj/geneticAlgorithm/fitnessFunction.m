function y = fitnessFunction(ut, U, G)
% Multiobjective FitnessFunction that minimizes the average error 
% between the gate transforation G across multiple eigenstates and 
% the total leakage probability over time 

% Initialize for two objectives 
y = zeros(1,2);

% Compute first objective
y(1) = QfitAvg(ut,U,G);

%Compute second objective
y(2) = leakProbAvg(ut, U); 

end

