b0 = beta*(0 :99);
%for i = 1:54
%    k0(i) = 2*(pi/2-0.01*i);
    
%end
%surf(b0,k0,abs(freq(1:54,1:53)))
figure
surf(b0,kinput,abs(freq(:,1:100)))
xlabel("kbeta")
ylabel("kd")
title("100 scatterers refractive index 3")


%spercific frequency plots
figure
plot(b0,abs(freq(40,1:100)))
xlabel("kbeta")
ylabel("response")
title("y slice at 32")
