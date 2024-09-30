function IFAvg = QfitAvg(ut,U,G)
% Used for Mixed Integer Problem 
% Evaluates the average error that the sequence has with the gate transformation
% across a defined number of eigenstates within the unitary time operator 
% database U

sum = 0;
size = length(U);

for i = 1:size
sum = sum + QfitRounded(ut, U{i}, G);
end

IFAvg = sum/size;

end

