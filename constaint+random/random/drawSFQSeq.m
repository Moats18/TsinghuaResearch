function drawSFQSeq(M, dt, sfq, uopt)
figure();
NN = length(uopt);
ut = [];

for i = 1:NN
    if uopt(i) <= 0 
        ut = [ut,zeros(1,M*abs(uopt(i)))];
    else
        ut = [ut,repmat(sfq, 1, uopt(i))];
    end
end

%{
xline(0.01,'Color','black','LineWidth',2);
hold on
xline(0.02,'Color','black','LineWidth',2);
hold on
xline(0,'Color','black','LineWidth',2);
hold on
xline(0.03,'Color','black','LineWidth',2);
hold on
%}
plot((1:length(ut))*dt,ut, 'Color','blue','LineWidth',1);

title("Pulse Sequence", 'FontSize', 14);
xlabel('Time in nanoseconds');
ylabel('Amplitude u rad/ns');

end

