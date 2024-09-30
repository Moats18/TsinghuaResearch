function c = maxPulseConstr(x)
maxPulses = 1000;
counter = 0;

for i = 1:length(x)
    counter = counter + abs(x(i)); 
end
c = counter-maxPulses;

end

