classdef solver_mie_penetrable2 < solver

    % ** Note: the class name above should match the filename and the name
    % of the constructor below **
    
    properties
        radius
        cof
        nmax
        refin
        rho
    end
    
    methods
        

        %-----------------------------------------
        % constructor
        %-----------------------------------------
        
        function self = solver_mie_penetrable2(kwave,incidentField,refin,rho,radius)

            % ** Note: the solver_name above should be replaced with the
            % name of your class. You may add other parameters as required.
            % **
            
            % call parent constructor
            self = self@solver(kwave,incidentField);
            
            self.refin = refin;
            
            self.radius = radius;
            self.rho = rho;
            
            self.nmax = suggestedorder(self.kwave,1);
            
        end
                 
        %===============================================================
        % these methods must be provided 
        %===============================================================

        %-----------------------------------------
        % setup
        %-----------------------------------------
        
        % This methods sets up your solver, eg assembles discretisation
        % matrices etc
        
        function setup(self)                      
            
            % ** your code goes here **
            
        end
        
        %-----------------------------------------
        % solve
        %-----------------------------------------
        
        % This method solves the scattering problem for every right hand
        % side specified in the self.incidentField cell array.
        
        function solve(self)

            
            self.refin
            
            self.cof = zeros((self.nmax+1)^2,self.numIncidentField);

            for j=1:self.numIncidentField

                v = zeros((self.nmax+1)^2,1);
                
                c = vec2cell(v);

                for n = 0:self.nmax

                    A = [ 
                        sphbesselj(n,self.refin*self.kwave*self.radius) -sphbesselh(n,self.kwave*self.radius)
                        (self.rho(1)/self.rho(2))*self.refin*self.kwave*sphbesseljd(n,self.refin*self.kwave*self.radius) -self.kwave*sphbesselhd(n,self.kwave*self.radius)
                        ];
                    
                    B = [
                        sphbesselj(n,self.kwave*self.radius)
                        self.kwave*sphbesseljd(n,self.kwave*self.radius)
                        ];
                    
                    X = A \ B;
                    
                    c{n+1} = X(2) * ones(2*n+1,1);

                end

                v = cell2vec(c);

                tmp = self.incidentField{j}.get_coefficients([0;0;0],self.nmax);
                
                self.cof(:,j) = tmp(:) .* v;
                
            end

        end
        
        %-----------------------------------------
        % get far field
        %-----------------------------------------
        
        % This method should compute the far field for the incident fields
        % self.incidentField{k} for each k in the array index. The return
        % value val should be an array with the column val(:,j) containing
        % the far field for self.incidentField{index(j)}.
        
        function val = getFarField(self,points,index)
            
            if nargin<3 || isempty(index)
                index = 1:self.numIncidentField;
            end
            
            r = sqrt(sum(points.^2,1));
            theta = acos(points(3,:)./r);
            phi = atan2(points(2,:),points(1,:));

            % evaluate spherical harmonic
            Y = associatedLegendre(self.nmax,cos(theta));            

            val = zeros(size(points,2),1);

            c = vec2cell(self.cof);

            for n=0:self.nmax

                j = -n:n;

                val = val + ((-1i).^(n+1)/self.kwave) * (Y.get(n) ...
                    .* exp(1i*phi(:)*j)) * c{n+1}(:,index);

            end
                    
            
        end
        
        %-----------------------------------------
        % description
        %-----------------------------------------
        
        function val = description(self)
            
            val = 'Mie series (penetrable).';
            
        end
        
        %===============================================================
        % you may provide other methods required to implement your solver
        % or help the user
        %===============================================================

    end % end methods
    
end