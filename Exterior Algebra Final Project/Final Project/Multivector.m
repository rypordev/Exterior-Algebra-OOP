classdef Multivector < handle
    properties (SetAccess = private, GetAccess = public)
        components
        degree
        dimension
        type
    end

    methods
        function self = Multivector(ExAlgebra, Degree, Components, Type)
            assert(isa(ExAlgebra,"ExteriorAlgebra"), ...
                "Ensure first argument is of type ExteriorAlgebra")
            
            s = size(ExAlgebra.permuteBasis{Degree},1);
            c = size(Components);
            if((c(1) ~= s|c(2)~=1) & ((c(1) ~= 1|c(2)~=s)))
                error("ERROR:\nComponents size for Degree %d in Dimension %d" + ...
                    " should be %d",Degree, ExAlgebra.dimension, s)
            end

            self.components = Components;
            if(nargin < 4)
                Type = MultivectorType.Contravector;
            end

            if(Type == MultivectorType.Covector)
                if(~isrow(self.components))
                    self.components = self.components';
                end
            elseif(Type == MultivectorType.Contravector)
                if(isrow(self.components))
                    self.components = self.components';
                end
            else
                error("ERROR:\nLabel Multivector as " + ...
                    "MultivectorType.Covector or MultivectorType.Contravector")
            end
            
            addlistener(ExAlgebra, 'ChangeBasisEvent', @self.basisChangeCallback);
            self.degree = Degree;
            self.dimension = ExAlgebra.dimension;
            self.type = Type;

        end

        function basisChangeCallback(self, ~, event)
            if(self.type == MultivectorType.Covector)
                self.components = (self.components/event.oldBasis{self.degree})*event.newBasis{self.degree};
            elseif(self.type == MultivectorType.Contravector)
                self.components = event.newBasis{self.degree}\(event.oldBasis{self.degree}*self.components);
            end
        end
    end
end