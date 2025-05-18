% Copyright 2014, 2015, 2016, 2017, 2018, 2022, 2023, 2024 Stuart C. Hawkins and M. Ganesh.

function [KBeta ,Ksp ,Freq,Betascale] = twodmultiscatterinputable3Dsawp(AmmountOfScatters,ScatterWidth,DistRadiusRatio,KInput)
%-----------------------------------------
% set main parameters
%-----------------------------------------
%ScatterWidth = 0.5;now an input
%kwave = 0.6952/ScatterWidth;%=1.3904
kwave = KInput/ScatterWidth;
% number of scatterers
numscat = AmmountOfScatters;
%kwave = 2.5888/2;
Ksp = KInput/ScatterWidth;
% setup incident plane wave
%p = plane_wave(pi/2,kwave);
% have seperate point source
p = point_source([-2 0 0],kwave);


%-----------------------------------------
% set up scatterers
%-----------------------------------------
FourierMeasureHeight = ScatterWidth*1.2;

Radius = ScatterWidth;
Density = 1.1;
RefractiveIndex = 3;
%-----------------------------------------
% setup the solvers
%-----------------------------------------

% setup solvers with empty incident field for now
% for each type of shape
SolvePeram = 30;
Location = [0,0,0];
AmmountUnique = 1;

%disp("Solver calulation started")
for k=1:AmmountUnique
    % create solver object
    %solver{k} = solverEllipsoidPenetrable(kwave,[],Location+[k,0,0],[Radius Radius Radius],RefractiveIndex,[1 Density],SolvePeram,2*SolvePeram);%solverNystromRobin(kwave,[],scatterer{k},0);
    solver{k} = solver_mie_penetrable2(kwave,[],RefractiveIndex,[1 Density],Radius);%solverNystromRobin(kwave,[],scatterer{k},0);
    % set Nystrom parameter
    solver{k}.setup
    
end

%disp("solver calulation ended")

%-----------------------------------------
% derived parameters
%-----------------------------------------
nmax = suggestedorder(kwave,Radius);
%assignin("base","nmax",nmax)
%-----------------------------------------
% setup the T-matrices
%-----------------------------------------
for k=1:AmmountUnique
    % setup T-matrix
    tmat{k} = ghtmatrix(nmax,kwave,solver{1},[0,0,0]);
    %assignin("base","Tobject",tmat)
    % print the T-matrix error estimate... based on Symmetry condition
    %fprintf('T-matrix %d (%s) error estimate %0.2e\n',k,...
    %    class(solver{k}.scatterer),tmat{k}.error());
end

%disp('Tmatrix calculated')

% set up arrays for scatterer types and positions... type(j) is a pointer
% into the cell array of possible shapes and sets the shape of the jth
% scatterer.
% pos(j) is the location of the jth scatterer (real part gives the x
% coordinate and imaginary part gives the y coordinate)
%type = zeros(numscat,1);
pos = zeros(numscat,3);%[];%= [zeros(numscat,1) zeros(numscat,1) zeros(numscat,1)];
xPosStart = 2;
 for q = 1:numscat
     pos(q,:) = [((ScatterWidth*DistRadiusRatio)*(q-1)+xPosStart) 0 0];
 end
 
 %assignin("base","Posistions",pos)
for j=1:numscat
    a{j} = regularwavefunctionexpansion(nmax,pos(j,:),p);
    %assignin("base","a",a)
end

% apply the T-matrices to the incident field
for j=1:numscat

    tmat{1}.origin = pos(j,:);%tmat{1}.setOrigin(pos(j,:));

    b{j} = tmat{1} * a{j};    
end

% Note: GMRES works with vectors so we need to convert our cell array of
% wavefunction expansions into a vector...
% setup an array right hand side coefficients (to pass to GMRES)
rhs = pack(b);
%assignin("base","rhs",rhs)
%assignin("base","b",b)

% set number of GMRES iterations
nitns = min(100,floor(numscat*((nmax+1)^2)/2));
%disp(nitns + "this is nitis")
% Solve the linear system using GMRES (the reshape functions are reshaping
% to 2D not 3D at the moment.
[x,flag,relres,iter,reshist] = gmres(@matrix_product,rhs,nitns,1e-8,1);

% convert coefficients into wavefunction expansions
c = unpack(x);


disp('Computed: Iteratively solved the multiple scattering model')

%-----------------------------------------
% visualize the scattered field
%-----------------------------------------
%figure(1)

% setup a grid
xstart = pos(1,1);

xstop = pos(numscat,1);
xsteps = 1000;
yheight = 2;
ysteps = 40;
zheight = 2; 
zsteps = 40;

tx=linspace(xstart,xstop,xsteps);
%assignin("base","Xspace",tx)
ty = linspace(-yheight,yheight,ysteps);
tz = linspace(-zheight,zheight,zsteps);
[x,y,z1] = meshgrid(tx,ty,tz);
Domain(1,:,:,:) = x;
Domain(2,:,:,:) = y;
Domain(3,:,:,:) = z1;
%maxrad = Radius;

function MaskedDomain = ThreeDmasking(x,y,z,SphereRadius,SphereCentreX,SphereCentreY,SphereCentreZ,Domain)
    [a1, b1, e] = meshgrid(x, y, z);
    Distance = sqrt((a1 - SphereCentreX).^2 + (b1 - SphereCentreY).^2 + (e - SphereCentreZ).^2);
    assignin("base","Distance",Distance);
    Mask = Distance >= SphereRadius;
    MaskedDomain = Mask;
end
%need to apply the mask correctly over every scatterer

PosStepsStart = 1;
PosStepsEnd = xsteps;

SampleHeight = round(FourierMeasureHeight/yheight*ysteps);
%assignin("base","linspace",tx)
%assignin("base","End",PosStepsEnd)
%assignin("base","Height",SampleHeight)
% get the largest radius of the scatterers

%for j=2:numscat
%    maxrad = max(maxrad,solver{type(j)}.getRadius());
%end

% get a mask for the scatterers

% get the scattered field... this is just the sum of the radiating fields
% from each scatterer
%disp(z)%
mask = ThreeDmasking(tx,ty,tz,Radius*0.9,pos(1,1),pos(1,2),pos(1,3),Domain);
scatfield = c{1}.evaluate(Domain,mask);

for j=2:numscat
      scatfield = scatfield + c{j}.evaluate(Domain,ThreeDmasking(tx,ty,tz,Radius*0.9,pos(j,1),pos(j,2),pos(j,3),Domain));
end

%-----------------------------------------
% visualize the total field
%-----------------------------------------
%disp('Evaluating and Visualizing (Fig. 2):')
%disp ('Output total field at 25,000 grid points (Fig. 2)')



% get the total field
totalfield = scatfield + p.evaluate(Domain,mask);
%assignin("base","Totalfield",totalfield)
%totalfield(round(SampleHeight/2)+round(ysteps/2),PosStepsStart:PosStepsEnd,round(zsteps/2)) = 100;

% plot the total field
%assignin("base","x",x)
%assignin("base","y",y)
%assignin("base","z1",z1)

%assignin("base","totalfield",totalfield)
surf(x(:,:,round(ysteps/2)),y(:,:,round(zsteps/2)),real(totalfield(:,:,round(zsteps/2))))
view([0 90]);
shading interp;

colorbar
title('Total field scatterers')
figure(2)
%-----------------------------------------
% indented function to implement the matrix
% product in GMRES
%-----------------------------------------
    function val = regularzero(n,x0,k)
% its actually (n+1)^2 throughout
        val = regularwavefunctionexpansion(n,x0,k,zeros((n+1)^2,1));
        %assignin("base","vec",vec)
        %assignin("base","n",n1)
        %val = regularwavefunctionexpansion(n,x0,k,zeros(2*n+1,1));
    end

    function y = matrix_product(x)
        
        % convert vector of coefficients into wavefunction expansions
        c = unpack(x);
        assignin("base","here",c)
        % apply matrix
        for j=1:numscat
            
            % we temporarily set the origin for the T-matrix object to the position
            % of the jth scatterer... this allows the T-matrix to interact with
            % wavefunction expansions with the same origin
            tmat{1}.origin = pos(j,:);
            
            
            % initialize sum
            csum = regularzero(nmax,pos(j,:),kwave);
            
            
            
            % sum contributions from the other scatterers
            for i=1:numscat
                
                if i~=j
                    
                    % get the expansion of c{i} at pos{j}
                    %disp(pos(i,:))
                    %disp(i)
                    %disp(pos(j,:))
                    %disp(j)
                    csum = csum + regularwavefunctionexpansion(c{i},pos(j,:));
                    
                end
                
            end
            
            % apply the T-matrix to the sum
            d{j} = c{j} - tmat{1} * csum;
            
        end
        
        % convert coefficients into a vector
        y = pack(d);
        
    end
%-----------------------------------------
% indented function to pack wavefunction
% expansion coefficients into a vector
%-----------------------------------------

    function vec = pack(a)
        
        % create an array to hold the coefficients
        vec = zeros((nmax+1)^2,numscat);

        % copy the coefficients into the array
        for j=1:numscat
            vec(:,j) = a{j}.getCoefficients();
        end
        
        % reshape the array into a vector
        vec = reshape(vec,[],1);
        
    end

%-----------------------------------------
% indented function to extract wavefunction
% expansion coefficients from a vector
%-----------------------------------------

    function a = unpack(vec)
        % reshape the vector into an array
        vec = reshape(vec,(nmax+1)^2,numscat);
        % create radiating wavefunction expansions from the columns of the
        % array
        for j=1:numscat
            a{j} = radiatingwavefunctionexpansion(nmax,pos(j,:),kwave,vec(:,j));
            %a{j} = wavefunctionexpansion(nmax,pos(j,:),kwave,vec(:,j));
        end
        
    end
%Need to make this 3D (might not being in the right place yet
    function Slice = FourierSliceTaker(Field,FirstPoint,LastPoint,HeightAbove,xsteps,ysteps,zsteps)
        Slice = Field(round(HeightAbove/2)+round(ysteps/2),FirstPoint:LastPoint,round(HeightAbove/2)+round(ysteps/2));
    end
%fourier(real(totalfield),10,450,1,0)
    function FreqeuencySpectrum = fourier(Sampleline)
        FreqeuencySpectrum = fft(Sampleline);
        Betascale = (-pos(1)+pos(2))*(2*pi*((xsteps)-1)/(xsteps*(pos(numscat)-pos(1))));
        n = length(FreqeuencySpectrum);
        %disp((1:n)*Betascale)
        %fshift = (-n/2:n/2-1)*((xsteps/2)/n);
        %yshift = fftshift(FreqeuencySpectrum);
        plot((0:30-1)*Betascale,abs(FreqeuencySpectrum(1:30)))
        %start from zero since fft starts from zero frequency
        figure(3)
        %plot((1:length(Sampleline))*((xsteps)/length(Samp(line)),abs(FreqeuencySpectrum))
        %xsteps/2 works and maybe since the transform is mirrored.
        %plot((0:length(FreqeuencySpectrum)-1)*(xsteps/2)/length(FreqeuencySpectrum),abs(FreqeuencySpectrum))
        %plot(fshift,abs(yshift)
        %xlabel("Frequency (scaled)")
        %ylabel("Frequency promanence this one")
        %ax = axis;
        %ax(1:2) = [0 2*pi];
        %axis(ax)
        %figure(3)

    end

FLine = FourierSliceTaker(real(totalfield),PosStepsStart,PosStepsEnd,SampleHeight,xsteps,ysteps,zsteps);
%assignin("base","FLine",FLine)
%figure(3)
%plot((1:length(FLine))*((pos(numscat)-pos(1)))/length(FLine),FLine)
%title("Offset sample")
%xlabel("X Location")
%ylabel("Wave amplitde")

FrequencySpec = fourier(FLine);

%test = BetaCal(PosStepsStart,PosStepsEnd,tx);
BetaScale = (2*pi*((xsteps)-1)/(xsteps*(pos(numscat,1)-pos(1,1))));
    function KValue = KCal(Frequency)
        BetaScale = (2*pi*((xsteps)-1)/(xsteps*(pos(numscat,1)-pos(1,1))));
        Positive = Frequency;%(150:end);
        [mag, index] = max(Positive);
        %disp(index)
        
        KValue = (index-1)*BetaScale;%1/krepi;
    end
KBeta = KCal(FrequencySpec);
Freq = FrequencySpec;
%figure(5)
%plot(test,abs(FrequencySpec))
end

