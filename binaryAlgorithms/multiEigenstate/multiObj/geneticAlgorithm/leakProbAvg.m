function leakProbAvg  = leakProbAvg(ut, U)
sum = 0;
size = length(U);

for i = 1:size
sum = sum + greatestLeakProb(ut, U{i});
end

leakProbAvg = sum/size;

end

