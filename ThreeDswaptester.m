
%addpath '/home/rhys/Documents/Tmatirix/TMATROM_OBJECT_ORIENTED_CORE'
%kvalues I see bloch waves for

%definetly need to write its as to only calculate the Tmatrix once.

%maybe should shorten the sampling area to where the bloch wave is more
%significant.
%what is the sampling frequency should be linxpace delta 
% the right peak for high resolution(one past pi) is larger, it might be
% because of the signs

for i = 3:33
tic
    [KBetaarray(i), KspArray(i)] = twodmultiscatterinputable3Dsawp(35,0.5,4,0.6952*0.05*i);%+1000*(i-1));%+0.79);%0.6952);
disp("This is the " + i + " iteration")
toc
end
%dont know dow to deal with the dropping out

KspArray(any(KBetaarray > 3.3, 1)) = [];
KBetaarray(any(KBetaarray > 3.3, 1)) = [];
figure
plot(2*KBetaarray,KspArray)
title("That one")
%I thing I have to multiply y by 2 to normalise to the scatter wideth and
%distance