% Point source incident field.
%
%  u = point_source(x0,k) returns a point source object u with 
%  wavenumber k and source location x0 (a vector with length 3).
%
% Also:
%
%   f = u.evaluate(p) returns the values f of the point source at points p.
%   Here p must be a 3 x n matrix.
%
%   f = u.evaluate(z,mask) returns the values f of the point source at
%   points z for which mask==1 and NaN elsewhere.
%
%   [dx,dy,dz] = u.evaluateGradient(p) returns dx, dy and dz the partial 
%   derivatives of the point source in the x, y and z directions respectively
%   at the points p. Here p must be a 3 x n matrix.
%
%   [dx,dy,dz] = u.evaluateGradient(z,mask) returns dx, dy and dz the partial 
%   derivatives of the point source in the x, y and z directions respectively
%   at the points p for which mask==1 and NaN elsewhere.
%
%   cof = u.get_coefficients(x0,n) returns the vector cof of regular
%   wavefunction expansion coefficients of the point source field with 
%   wavefunction origin x0 and order n.
%
% See also: plane_wave, incident.
%
% Stuart C. Hawkins - 20 April 2021

classdef point_source < incident
    
    properties
        location
    end
    
    methods
       
        %-----------------------------------------
        % constructor
        %-----------------------------------------

        function self = point_source(location,kwave)           
            
            % call parent constructor
            self = self@incident(kwave);
            
            % set location (as column vector)
            self.location = location(:);
            
        end
        
        %-----------------------------------------
        % evaluate
        %-----------------------------------------

        function val = evaluate(self,points,mask)
            
            % get size of points
            n = size(points);
            
            % reshape points
            points = reshape(points,3,[]);
            
            % reshape mask if given
            if nargin>2
                mask = reshape(mask,1,[]);
            end
            
            % initialise return array
            val = zeros(1,size(points,2));
            
            % apply mask if necessary
            if nargin > 2
                points = points(:,mask);
            end
            
            % get distance of points from source
            r = sqrt(sum((points - repmat(self.location,1,size(points,2))).^2,1));
            
            % evaluate incident field
            v = (0.25/pi)*exp(1i*self.kwave*r)./r;
            
            % insert values into the return array
            if nargin>2
                val(:,mask) = v;
            else
                val = v;
            end
            
            % reshape the return array to match points
            if length(n)==2
                val = reshape(val,n(end),1);
            else
                val = reshape(val,n(2:end));
            end
            
        end
        
        %-----------------------------------------
        % evaluate gradient
        %-----------------------------------------

        function [dx,dy,dz] = evaluateGradient(self,points,mask)
            
            % get size of points
            n = size(points);
            
            % reshape points
            points = reshape(points,3,[]);
            
            % reshape mask if given
            if nargin>2
                mask = reshape(mask,1,[]);
            end
            
            % initialise return array
            dx = zeros(1,size(points,2));
            dy = zeros(1,size(points,2));
            dz = zeros(1,size(points,2));
            
            % apply mask if necessary
            if nargin > 2
                points = points(:,mask);
            end
            
            % get direction of points from source
            d = points - repmat(self.location,1,size(points,2));
            
            % get distance of points from source
            r = sqrt(sum(d.^2,1));

            % evaluate incident field
            v = (0.25/pi)*exp(1i*self.kwave*r)...
                .*((1i*self.kwave)*r-1)./r.^3;
            
            % insert values into the return array
            if nargin>2
                dx(:,mask) = v.*d(1,:);
                dy(:,mask) = v.*d(2,:);
                dz(:,mask) = v.*d(3,:);
            else
                dx = v.*d(1,:);
                dy = v.*d(2,:);
                dz = v.*d(3,:);
            end
            
            % reshape the return array to match points
            if length(n)==2
                dx = reshape(dx,1,n(end));
                dy = reshape(dy,1,n(end));
                dz = reshape(dz,1,n(end));
            else
                dx = reshape(dx,n(2:end));
                dy = reshape(dy,n(2:end));
                dz = reshape(dz,n(2:end));
            end
            
        end
        
        %-----------------------------------------
        % evaluate far field
        %-----------------------------------------

        function val = evaluateFarField(self,points)
        
            % compute the far field using (2.15) in Colton and Kress,
            % Inverse Acoustic and Electromagnetic Scattering Theory, 4th
            % edition.
            dp = sum(repmat(self.location,1,size(points,2)) .* points,1);            
            val = (0.25/pi)*exp(-1i*self.kwave*dp);
            
        end
        
        
        %-----------------------------------------
        % get coefficients
        %-----------------------------------------

        function cof = get_coefficients(self,centre,nmax)

            % compute translation direction
            x = self.location(:)-centre(:);

            % convert to spherical polar coordinates
            r = sqrt(sum(x.^2));
            theta = acos(x(3)/r);
            phi = atan2(x(2),x(1));

            % get vector of orders
            n = 0:nmax;

            % get Bessel function
            J = sphbesselh(n,self.kwave*r);

            % compute the spherical harmonic part
            Y = associatedLegendre(nmax,cos(theta));

            % loop through the orders
            for n=0:nmax

                % compute coefficients using Colton and Kress, Inverse Acoustic
                % and Electromagnetic Scattering Theory, Equation (2.43)
                % Note: the conjugate of the spherical harmonic is included
                % here
                tmp{n+1} = 1i * self.kwave * J(n+1) * Y.get(n) .* exp(-1i*(-n:n)*phi);

            end

            % convert cell to vector
            cof = cell2vec(tmp);
            
        end
        
    end
    
end
       