% Incident field class
%
%  This abstract class provides the key interfaces for objects of class 
%  incident. These are used by objects of solver class for creating T-matrices.
%
%  Note: this is an abstract class, it cannot be instantiated.
%
% See also: plane_wave, point_source.
%
% Stuart C. Hawkins - 20 April 2021


classdef (Abstract) incident < handle
    
    properties
        kwave
    end
    
    methods
        
        %-------------------------------------------------
        % constructor
        %-------------------------------------------------

        function self = incident(kwave)

            % store wavenumber
            self.kwave = kwave;

        end

        %-------------------------------------------------
        % overloaded +
        %-------------------------------------------------
        
        function obj = plus(self,other)
            
            % check that input object is of class incident
            if ~isa(self,'incident')
                error('self must be of class incident')
            end
        
            % check that input object is of class incident
            if ~isa(other,'incident')
                error('other must be of class incident')
            end

            % create instance of incidentplus class
            obj = incidentplus(self,other);
            
        end
    
        %-------------------------------------------------
        % overloaded -
        %-------------------------------------------------
        
        function obj = minus(self,other)
            
            % check that input object is of class incident
            if ~isa(self,'incident')
                error('self must be of class incident')
            end
        
            % check that input object is of class incident
            if ~isa(other,'incident')
                error('other must be of class incident')
            end

            % create instance of incidentminus class
            obj = incidentminus(self,other);
            
        end
        
        function obj = uminus(self)
            
            % check that input object is of class incident
            if ~isa(self,'incident')
                error('self must be of class incident')
            end
        
            % create instance of incidentminus class
            obj = incidentuminus(self);
            
        end
        
        %-------------------------------------------------
        % overloaded *
        %-------------------------------------------------

        % This provides facility to multiply by a scalar
        
        function obj = mtimes(obj1,obj2)
            
            % check that one of the input objects is of class incident
            if ~isa(obj1,'incident') && ~isa(obj2,'incident')
                % This should never happen... this method will only be
                % called if one of the arguments is of class incident
                error('one of arguments must be of class incident')
            end
        
            % check that one of the input object is a double
            if ~isa(obj1,'double') && ~isa(obj2,'double')
                error('one of the arguments must be of class double')
            end

            % create instance of incidenttimes class
            if isa(obj1,'incident')
                obj = incidenttimes(obj1,obj2);
            else
                obj = incidenttimes(obj2,obj1);
            end
            
        end

    end
    
    methods(Abstract=true)

        % these must be overloaded in the child class
        cof = get_coefficients(self,centre,nmax);
        val = evaluate(self,points,mask);
                
    end
    
end