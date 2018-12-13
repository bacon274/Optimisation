

classdef window
    properties
        transparency
        area
        g 
        u 
        power
        cost
    end
    methods
        function obj = window(A,u)
            if nargin == 0
                obj.area = 1;
                obj.u = 1;

            else
                obj.area = A;
                obj.u = u;
            end
        
        end
    end
    
end