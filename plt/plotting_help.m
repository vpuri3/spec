% Plotting
close all % Close all currently open figures
hold on % Plot on the same figure
% Colors: b = blue, k = black, r = red, g = green, y = yellow
% Styling: o = circles, - = solid line, -- = dashed line
plot(x,y,'bo',var,P,'r--')  % Plot...add as many x and y combinations as needed 
legend('Data Points','Langrange Approximation') % names in order of plotting above


loglog , semilogx , semilogy % logarithmic, semi-logarithmic plots
hold on % between plot calls would create different figures.


%% Plotting with loops.
hold on
for c = 1:d
    len = 5/h(c)+1;
    plot(T(1:len,c),er_F(1:len,c),'-o','Linewidth',1.75,'DisplayName',['h=' num2str(h(c))])
end
title('Absolute Error in Forward Euler with varied grid spacing (h)','FontSize',14)
xlabel('Time','FontSize',14)
ylabel('Absolute Error','FontSize',14)
axis sqaure
legend('show')