sum = 0;
cycles = 10;
k = 100;
NN = 30;
x0 = zeros(1, NN);


GA = zeros(1, cycles);
for z = 1:cycles
   GA(z) = GAcrossoverFcn(k,NN, );
 end

 best = GA(1);
   for z = 1:cycles
        if GA(z) < best
            best = GA(z);
        end
     sum = sum + GA(z); 
   end

