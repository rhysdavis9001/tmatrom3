
%addpath '/home/rhys/Documents/Tmatirix/TMATROM_OBJECT_ORIENTED_CORE'
%kvalues I see bloch waves for

%definetly need to write its as to only calculate the Tmatrix once.

%maybe should shorten the sampling area to where the bloch wave is more
%significant.
%what is the sampling frequency should be linxpace delta 
% the right peak for high resolution(one past pi) is larger, it might be
% because of the signs betascale = 0.0923
%take alot of slice with alot of scatterers near kd = pi (the cuttoff)
%so start at k = pi/2 and slowly go down
numscat = 100;
scatterwidth = 0.7;
centertocenter = 2;

for i = 1:101
tic
    %[KBetaarray(i), KspArray(i) ,freq(i,:)] = twodmultiscatterinputable3Dsawp(35,0.5,4,0.6952*0.05*i);%+1000*(i-1));%+0.79);%0.6952);
    kinput(i) = (pi/centertocenter*(0.01*i));
    [KBetaarray(i), KspArray(i) ,freq(i,:),beta,totalfield(i)] = twodmultiscatterinputable3Dsawp(numscat,scatterwidth,centertocenter,kinput(i));
    %kvalue = 1000*pi/2-0.01*i;
    save(numscat + "scatterers"+ scatterwidth+"length"+ centertocenter+ "distancebetweenscatters" +i +"iteration")
    
   
    disp("This is the " + i + " iteration")
toc
end
%dont know dow to deal with the dropping out

KspArray(any(KBetaarray > 3.3, 1)) = [];
KBetaarray(any(KBetaarray > 3.3, 1)) = [];
figure
plot(KBetaarray,KspArray)
%title("penatrable spheres density = 10 refractive index = 10")
xlabel("kbeta" ...
    )
ylabel("kinc")
for i= 1:25

    freq(:,i) = freq(:,i)*beta;
end
b0 = 0.0923*(0 :39);
k0 = 0.6952*0.05*(1:54);
surf(b0,k0,abs(freq(1:54,1:40)))
%I thing I have to multiply y by 2 to normalise to the scatter wideth and
%distance