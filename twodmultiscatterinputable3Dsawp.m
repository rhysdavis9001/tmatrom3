% Copyright 2014, 2015, 2016, 2017, 2018, 2022, 2023, 2024 Stuart C. Hawkins and M. Ganesh.
% 	
% This file is part of TMATROM.
% 
% TMATROM is free software: you can redistribute it and/or modify	
% it under the terms of the GNU General Public License as published by	
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% TMATROM is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% a change
% You should have received a copy of the GNU General Public License
% along with TMATROM.  If not, see <http://www.gnu.org/licenses/>.

function [KBeta ,Ksp] = twodmultiscatterinputable(AmmountOfScatters,ScatterWidth,DistRadiusRatio,KInput)

%-----------------------------------------
% set main parameters
%-----------------------------------------
%ScatterWidth = 0.5;now an input
% wavenumber
%kwave = 0.6952/ScatterWidth;%=1.3904
kwave = KInput/ScatterWidth;
%kwave = 2.5888/2;
Ksp = KInput/ScatterWidth;
% setup incident plane wave
%p = plane_wave(pi/2,kwave);
% have seperate point source
p = point_source([-10 0 0],kwave);

%-----------------------------------------
% set up scatterers
%-----------------------------------------
FourierMeasureHeight = ScatterWidth*1.2;
%scatterer{1} = obstacleCircle(ScatterWidth); It appears thattheopject is
%defined in the solver

Radius = 1;
Density = 1;
RefractiveIndex = 1;
%-----------------------------------------
% setup the solvers
%-----------------------------------------

% setup solvers with empty incident field for now
% for each type of shape
SolvePeram = 30;
Location = [0,0,0];
AmmountUnique = 1;
for k=1:AmmountUnique
    % create solver object
    %solver{k} = solverEllipsoidPenetrable(kwave,[],Location+[k,0,0],[Radius Radius Radius],RefractiveIndex,[1 Density],SolvePeram,2*SolvePeram);%solverNystromRobin(kwave,[],scatterer{k},0);
    solver{k} = solver_mie_penetrable2(kwave,[],RefractiveIndex,[1 Density],Radius);%solverNystromRobin(kwave,[],scatterer{k},0);
    % set Nystrom parameter
    solver{k}.setup
    
end
%-----------------------------------------
% derived parameters
%-----------------------------------------
nmax = suggestedorder(kwave,Radius);
assignin("base","nmax",nmax)
%-----------------------------------------
% setup the T-matrices
%-----------------------------------------
for k=1:AmmountOfScatters
    % setup T-matrix
    tmat{k} = ghtmatrix(nmax,kwave,solver{1},[0,0,0]);
    assignin("base","Tobject",tmat)
    % print the T-matrix error estimate... based on Symmetry condition
    %fprintf('T-matrix %d (%s) error estimate %0.2e\n',k,...
    %    class(solver{k}.scatterer),tmat{k}.error());
end

disp('Computed: Input/output/multiple-configuration  independent')
disp('ROM characterization of several scatterers')

% A three scatterer multiple particle configuration 
% with centers (-2, -1), (2, -1), and (0,2)
%disp('Multiple config.: Choosing three scatterers (3,4,5) and their locations')


% number of scatterers
numscat = AmmountOfScatters;

% set up arrays for scatterer types and positions... type(j) is a pointer
% into the cell array of possible shapes and sets the shape of the jth
% scatterer.
% pos(j) is the location of the jth scatterer (real part gives the x
% coordinate and imaginary part gives the y coordinate)
%type = zeros(numscat,1);
pos = zeros(numscat,3);%[];%= [zeros(numscat,1) zeros(numscat,1) zeros(numscat,1)];
xPosStart = 2;
 for q = 1:numscat
     %type(q) = 1;
     pos(q,:) = [((ScatterWidth*DistRadiusRatio)*(q)+xPosStart) 0 0];
 end
 assignin("base","Posistions",pos)
for j=1:numscat
    a{j} = regularwavefunctionexpansion(nmax,pos(j,:),p);
    assignin("base","a",a)
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
assignin("base","rhs",rhs)
assignin("base","b",b)

% set number of GMRES iterations
nitns = min(100,floor(numscat*(2*nmax+1)/2));

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
xstart = pos(numscat);

xstop = pos(1);
xsteps = 30;
yheight = 5;
ysteps = 75;
zheight = 5; 
zsteps = 75;

tx=linspace(xstart,xstop,xsteps);
%assignin("base","Xspace",tx)
ty = linspace(-yheight,yheight,ysteps);
tz = linspace(-zheight,zheight,zsteps);
[x,y,z1] = meshgrid(tx,ty,tz);
Domain(1,:,:,:) = x;
Domain(2,:,:,:) = y;
Domain(3,:,:,:) = z1;%= meshgrid(tx,ty,tz);
%z = x+y*1i;
maxrad = Radius;

function MaskedDomain = ThreeDmasking(x,y,z,SphereRadius,SphereCentreX,SphereCentreY,SphereCentreZ,Domain)
    [a, b, e] = meshgrid(x, y, z);
    Distance = sqrt((a - SphereCentreX).^2 + (b - SphereCentreY).^2 + (e - SphereCentreZ).^2);
    assignin("base","Distance",Distance);
    Mask = Distance >= SphereRadius;

    %MaskedDomain = Domain.*Mask;
    MaskedDomain = Mask;
end

%mask = abs(z-pos(1)) < 1.1 * maxrad;

mask = ThreeDmasking(tx,ty,tz,Radius*1.1,pos(1,1),pos(1,2),pos(1,3),Domain);

%for j=1:numscat
%    mask = mask | abs(z-pos(j)) < 1.1 * maxrad;
%end

PosStepsStart = 1;%round(pos(1)/xstart*xsteps);
PosStepsEnd = xsteps;%xsteps;%round(pos(numscat)/xstart*xsteps);

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
z = class(c(1));
disp("here")
assignin("base", "C",c)
disp(z)%scatfield = c{1}.evaluate(z,~mask);
scatfield = c{1}.evaluate(Domain,mask);
disp("This")
%BE CAREFUL SINCE ONLY ONE OBJECT curly braces do not wor yet.

%for j=2:numscat
%    scatfield = scatfield + c.evaluate(Domain,mask);
%end

%-----------------------------------------
% visualize the total field
%-----------------------------------------
disp('Evaluating and Visualizing (Fig. 2):')
disp ('Output total field at 25,000 grid points (Fig. 2)')



% get the total field
totalfield = scatfield + p.evaluate(Domain,mask);
%assignin("base","Totalfield",totalfield)
%totalfield(round(SampleHeight/2)+round(ysteps/2),PosStepsEnd:PosStepsStart) = 100;

% plot the total field
surf(x(:,:,34),y(:,:,34),real(totalfield(:,:,34)))
view([0 90]);
shading interp;
figure(2)
%colorbar
%title('Total field ')
% add visualisation of scatterers
%hold on
%pp = 2*pi*(0:99)/100;
%for j=1:numscat
%    [sx,sy,qx,qy]=scatterer{type(j)}.geom(pp);
%    obj = plot(qx+real(pos(j)),qy+imag(pos(j)),'k-');
%    set(obj,'linewidth',2);
%end
%SliceCheck = line([pos(numscat),pos(1)],[FourierMeasureHeight,FourierMeasureHeight])
%hold off
%SliceMatrixStart = pos(1)
%-----------------------------------------
% indented function to implement the matrix
% product in GMRES
%-----------------------------------------
    function val = regularzero(n,x0,k)
% its actually (n+1)^2 throughout
        val = regularwavefunctionexpansion(n,x0,k,zeros((n+1)^2,1));
    end
    function y = matrix_product(x)
        
        % convert vector of coefficients into wavefunction expansions
        c = unpack(x);
        
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
        end
        
    end

    function Slice = FourierSliceTaker(Field,FirstPoint,LastPoint,HeightAbove,xsteps,ysteps)
        %disp(FirstPoint)
        %disp(LastPoint)
        Slice = Field(round(HeightAbove/2)+round(ysteps/2),FirstPoint:LastPoint);
    end
%fourier(real(totalfield),10,450,1,0)
    function FreqeuencySpectrum = fourier(Sampleline)
        FreqeuencySpectrum = fft(Sampleline);
        %disp("here")
        Betascale = (pos(1)-pos(2))*(2*pi*((xsteps)-1)/(xsteps*(-pos(numscat)+pos(1))));
        %disp(Betascale)
        %disp(FreqeuencySpectrum)
        n = length(FreqeuencySpectrum);
        %disp((1:n)*Betascale)
        
        %fshift = (-n/2:n/2-1)*((xsteps/2)/n);
        %yshift = fftshift(FreqeuencySpectrum);
        plot((0:n-1)*Betascale,abs(FreqeuencySpectrum))
        %start from zero since fft starts from zero frequency
        %figure 
        %plot((1:length(Sampleline))*((xsteps)/length(Sampleline)),abs(FreqeuencySpectrum))
        %xsteps/2 works and maybe since the transform is mirrored.
        %plot((0:length(FreqeuencySpectrum)-1)*(xsteps/2)/length(FreqeuencySpectrum),abs(FreqeuencySpectrum))
        %plot(fshift,abs(yshift))
        %title('Fourier Transform')
        xlabel("Frequency (scaled)")
        ylabel("Frequency promanence")
        ax = axis;
        ax(1:2) = [0 2*pi];
        axis(ax)
        figure(3)

    end

FLine = FourierSliceTaker(real(totalfield),PosStepsStart,PosStepsEnd,SampleHeight,xsteps,ysteps);
%assignin("base","FLine",FLine)
%figure(3)
%plot((1:length(FLine))*((pos(numscat)-pos(1)))/length(FLine),FLine)
%title("Offset sample")
%xlabel("X Location")
%ylabel("Wave amplitde")

FrequencySpec = fourier(FLine);
%assignin("base","FrequencySpec",FrequencySpec)

%    function Beta = BetaCal(spacestart,spaceend,xlinspace)
        %so beta is something that actually plots x axis against the
        %fft.need slice
        %numsample = spacestart-spaceend;
%        i = 1;
%        for n = 1:length(xsteps)
%            
%            Beta(i) = 2*pi*((xlinspace(n)*(xsteps-1))/((pos(numscat)-pos(1))*xsteps));
%            i = i +1;
%        end
%    end 
%test = BetaCal(PosStepsStart,PosStepsEnd,tx);
Betascale = (2*pi*((xsteps)-1)/(xsteps*(-pos(numscat)+pos(1))));
    function KValue = KCal(Beta,Frequency)
        Positive = Frequency;%(150:end);
        [mag, index] = max(Positive);
        disp(index)
        krepi = Beta(index);
        KValue = krepi;%1/krepi;
    end
KBeta = 1;%KCal(Betascale ,FrequencySpec);
%figure(5)
%plot(test,abs(FrequencySpec))
end

