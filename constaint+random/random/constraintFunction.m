function [c,ceq] = constraintFunction(maxPulses,x)
counter = 0;

for i = 1:length(x)
    counter = counter + abs(x(i)); 
end
c = counter-maxPulses;
ceq = [];

end

