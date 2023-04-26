classdef ExteriorAlgebra < handle
    properties (SetAccess = private, GetAccess = public)
        gradeBasis
        gradeMetric
        permuteBasis
        dimension
        hodgeArr
    end
    properties (SetAccess = public, GetAccess = public)
    end

    methods (Access = public)
        %Input Basis and Metric are for grade 1 vector space
        function self = ExteriorAlgebra(Dim, Basis, Metric)
            self.dimension = Dim;
            if(self.dimension > 32)
                error("ERROR:\nDimension %s is too large. Currently " + ...
                    "only supports up to Dim 32", self.dimension);
            end

            if(nargin < 2)
                Basis = eye(Dim);
            else
                Basis = ExteriorAlgebra.CheckBasis(Basis, Dim);
            end

            if(nargin < 3)
                Metric = eye(Dim);
            else
                Metric = ExteriorAlgebra.CheckMetric(Metric, Dim);
            end
            
            self.permuteBasis = ExteriorAlgebra.CreateGradeBasisPermute(self.dimension);
            self.gradeBasis = self.CreateGradeBasis(Basis);
            self.gradeMetric = self.CreateGradeMetrics(Metric);
            self.hodgeArr = self.CreateHodgeArrays();
        end

        function obj = CreateMultivectorCurrentBasis(self, Degree, Components, Type)
            assert(isa(Type,"MultivectorType"), ...
                "Ensure third argument is of type MultivectorType")

            %Check Components Length matches permuteBasis{Degree} length
            s = size(self.permuteBasis{Degree},1);
            c = size(Components);
            if((c(1) ~= s|c(2)~=1) & ((c(1) ~= 1|c(2)~=s)))
                error("ERROR:\nComponents size for Degree %d in Dimension %d" + ...
                    " should be %d",Degree,self.dimension,s)
            end

            obj = Multivector(self, Degree, Components, Type);
        end
        function obj = CreateMultivectorStandardBasis(self, Degree, Components, Type)
            assert(isa(Type,"MultivectorType"), ...
                "Ensure third argument is of type MultivectorType")

            %Check Components Length matches permuteBasis{Degree} length
            s = size(self.permuteBasis{Degree},1);
            c = size(Components);
            if((c(1) ~= s|c(2)~=1) & ((c(1) ~= 1|c(2)~=s)))
                error("ERROR:\nComponents size for Degree %d in Dimension %d" + ...
                    " should be %d",Degree,self.dimension,s)
            end

            %Convert Components from Standard Basis
            cellWrap = cell(Degree);
            cellWrap{Degree} = eye(s);
            obj = Multivector(self, Degree, Components, Type);
            obj.basisChangeCallback('', ChangeBasisData(cellWrap, self.gradeBasis));
        end
        function val = EvalMultivector(self, mv)
            if(mv.type == MultivectorType.Covector)
                val = mv.components/self.gradeBasis{mv.degree};
            elseif(mv.type == MultivectorType.Contravector)
                val = self.gradeBasis{mv.degree}*mv.components;
            end
        end

        function ChangeBasis(self, newBasis)
            newBasis = ExteriorAlgebra.CheckBasis(newBasis, self.dimension);
            
            newGradeBasis = self.CreateGradeBasis(newBasis);
            
            notify(self, 'ChangeBasisEvent', ChangeBasisData(self.gradeBasis, newGradeBasis));
            self.gradeBasis = newGradeBasis;
            self.CreateHodgeArrays();
        end
        function ChangeMetric(self, newMetric)
            newMetric = ExteriorAlgebra.CheckMetric(newMetric, self.dimension)

            newGradeMetric = self.CreateGradeMetrics(newMetric);
            self.hodgeArr = self.CreateHodgeArrays();
            self.gradeMetric = newGradeMetric;
        end

        function DispBasis(self)
            disp("Grade Basis:")
            g = self.gradeBasis;
            for i=1:size(g,1)
                fprintf("%d Basis\n", i)
                disp(g{i})
            end
        end
        % 'std', 'co', 'con'
        function DispStdMetric(self)
            m = self.gradeMetric;
            disp("Standard Basis Grade Metrics:")
            for i=1:size(m,1)
                fprintf("%d Metric\n", i)
                disp(m{i})
            end
        end
        function m = GetMetric(self, grade, mv1Type, mv2Type)
            m = self.gradeMetric{grade};
            Q = self.gradeBasis{grade};
            if mv1Type == MultivectorType.Covector
                m = Q\m;
            else
                m = Q'*m;
            end

            if mv2Type == MultivectorType.Covector
                m = m/(Q');
            else
                m = m*Q;
            end
        end
    
        function mv3 = Wedge(self, mv1, mv2)
            assert(isa(mv1,"Multivector"), ...
                "Wedge takes Multivectors as inputs")
            assert(isa(mv2,"Multivector"), ...
                "Wedge takes Multivectors as inputs")
            assert(mv1.dimension == self.dimension, ...
                "Multivector A dimension does not match %d",self.dimension)
            assert(mv2.dimension == self.dimension, ...
                "Multivector B dimension does not match %d",self.dimension)
            assert(mv1.type == mv2.type, ...
                "Ensure Multivector types match")

            newDegree = mv1.degree + mv2.degree;
            if(newDegree > self.dimension)
                mv3 = 0;
                warning("Wedge Product Resulted in Multivector Degree %d, returning zero.", newDegree)
                return
            end

            xb = self.permuteBasis{mv1.degree};
            yb = self.permuteBasis{mv2.degree};
            zb = self.permuteBasis{newDegree};

            x = mv1.components;
            y = mv2.components;
            z = zeros(size(zb,1),1);

            I = eye(newDegree);

            for n = 1:size(xb,1)
                for m=1:size(yb,1)
                    [bnm, perm] = sort([xb(n,:) yb(m,:)]);
                    if prod(diff(bnm)) ~= 0
                        [member,k] = ismember(bnm,zb,'rows');
                        if member
                            z(k) = z(k)+x(n)*y(m)*det(I(:,perm));
                        end
                    end
                end
            end
            mv3 = Multivector(self, newDegree, z, mv1.type);
        end
        function result = InnerProduct(self, mv1, mv2)
            assert(isa(mv1,"Multivector"), ...
                "InnerProduct takes Multivectors as inputs")
            assert(isa(mv2,"Multivector"), ...
                "InnerProduct takes Multivectors as inputs")
            assert(mv1.dimension == self.dimension, ...
                "Multivector A dimension does not match %d",self.dimension)
            assert(mv2.dimension == self.dimension, ...
                "Multivector B dimension does not match %d",self.dimension)
            assert(mv2.degree == mv1.degree, ...
                "Multivector degrees do not match")
    
            x = mv1.components;
            y = mv2.components;
            m = self.GetMetric(mv1.degree, mv1.type, mv2.type);

            n = size(m,1);
            result = 0;

            if ~isrow(x)
                x = x';
            end
            if isrow(y)
                y = y';
            end

            result = x*m*y;
        end
        
        function result = HodgeStar(self, mv)
            newDegree = self.dimension - mv.degree;
            %Go Component-wise, multiply 
            if(mv.type == MultivectorType.Covector)
                h = self.hodgeArr{1, mv.degree};
                cStar = mv.components*h;
                result = self.CreateMultivectorCurrentBasis(newDegree, cStar, mv.type);
            else
                h = self.hodgeArr{2, mv.degree};
                cStar = h*mv.components;
                result = self.CreateMultivectorCurrentBasis(newDegree, cStar, mv.type);
            end
        end
    end

    %Create Grade Basis & Metric
    methods (Access = protected)
        function metric = CreateGradeMetrics(self, m)
            metric = cell([self.dimension 1]);
            metric{1} = m;

            for k=2:self.dimension
                kbasis = self.permuteBasis{k};
                s2 = size(kbasis,1);
                kmetric = zeros(s2);
                for x=1:s2
                    for y=x:s2
                        a = kbasis(x,:);
                        b = kbasis(y,:);
                        submatrix = zeros(k);
                        for i=1:k
                            for j=1:k
                                submatrix(i,j) = m(a(i),b(j));
                            end
                        end
                        d = det(submatrix);
                        kmetric(x,y) = d;
                        kmetric(y,x) = d;
                    end
                end
                metric{k} = kmetric;
            end
        end
        function basis = CreateGradeBasis(self, b)
            basis = cell([self.dimension 1]);
            basis{1} = b;

            %Future: Vectorize CoB matrix Computation
            for i=2:self.dimension
                iBasis = self.permuteBasis{i};
                n = size(iBasis,1);
                Q = zeros(n);
                for x=1:n
                    for y = 1:n
                        %Find P_xy Here
                        subMatrix = zeros(i);
                        pRow = iBasis(x,:);
                        pCol = iBasis(y,:);
                        for x1 = 1:i
                            for y1 = 1:i
                                subMatrix(x1,y1) = b(pRow(x1),pCol(y1));
                            end
                        end
                        Q(x,y) = det(subMatrix);
                    end
                end
                basis{i} = Q;
            end
        end
        function hodge = CreateHodgeArrays(self)
            %Hodge star is invariant of basis. These hodge arrays are
            %calculated from the metric in standard basis

            %Hodge cells is 2xn
            % first row is covariant hodge values
            % second row is contravariant hodge values
            % new index after hodge is n-i+1
            dim = self.dimension;
            hodge = cell(2,dim-1);
            options = 1:dim;

            for k = 1:dim-1
                basis = self.permuteBasis{k};
                n = size(basis,1);

                hCo = zeros(n);
                hContra = zeros(n);

                mCo = self.GetMetric(dim-k, MultivectorType.Covector, MultivectorType.Covector);
                mContra = self.GetMetric(dim-k, MultivectorType.Contravector, MultivectorType.Contravector);

                I = eye(dim);
                for i=1:n
                    starIndex = n-i+1;
                    sb = setdiff(options, basis(i,:));
                    [~, perm] = sort([basis(i,:) sb]);
                    %Hodge star indexes are reversed of normal index
                    sign = det(I(:,perm));

                    YVector = zeros(1,n);
                    YVector(starIndex) = 1;
                    hCo(i,:) = sign.*(YVector/mCo);
                    hContra(:,i) = sign.*(YVector/mContra);
                end
                hodge{1, k} = hCo;
                hodge{2, k} = hContra;
            end
        end
    end

    methods (Static, Access=protected)
        function kbasis = CreateGradeBasisBinary(dim)
            kbasis = cell([dim 1]);

            % Pre-allocate cell array lengths
            %for i=1:dim
            %    nchoosek(dim,i)
            %    kbasis{i} = NaN(1,nchoosek(dim,i))
            %end
            

            helper = 1;
            indices = 1;
            % Generate indexing
            for i=1:dim-1
                helper = [helper, helper+ones([1,2^(i-1)])];
                indices = [indices, helper];
            end

            % Populate k-th basis using above indexing
            for i=1:2^(dim)-1
                x = indices(i);
                kbasis{x} = [kbasis{x} uint32(i)];
            end
        end
        function kbasis = CreateGradeBasisPermute(dim)
            kbasis = cell([dim 1]);
            for i=1:dim
                kbasis{i} = nchoosek(1:dim,i);
            end
        end

        function basis = CheckBasis(nbasis, dim)
            if(issparse(nbasis))
                nbasis = full(nbasis);
            end
            basis = nbasis;

            s1 = size(basis,1);
            s2 = size(basis,2);
            
            if(s1 ~= dim)
                error("ERROR:\nBasis size should match dimension\n" + ...
                    "Current matrix is [%d by %d]",s1, s2)
            end
            
            if(s1 ~= s2)
                error("ERROR:\nPlease use a square matrix as your basis\n" + ...
                    "Current matrix is [%d by %d]",s1, s2)
            end

            if(rank(basis) ~= dim)
                error("ERROR:\nBasis must have rank = dimension\n" + ...
                    "Current matrix has rank %d and dim %d",rank(basis), s1)
            end
        end
        function metric = CheckMetric(m, dim)
            if(issparse(m))
                m = full(m);
            end
            metric = m;
            
            assert(size(metric,2) == size(metric,1), "Please use a square matrix as " + ...
                "your metric\nCurrent matrix is [%d by %d]",size(metric,1), size(metric,2))
            assert(size(metric,1) == dim, "Metric must have size = dimension\n" + ...
                    "Current matrix has size %d, Algebra dim is %d",size(metric,1),dim)
            assert(issymmetric(metric), "Metric must be symmetric.")
        end

        %Unused for fear of reprecussions
        function sign = GetWedgeSign(b1, b2)
            count = 0;
            flips = 0;
            k = 1;
            s1 = size(b1,2);
            s2 = size(b2,2);
            for i=1:s1
                while b1(i) > b2(k)
                    count = count+1;

                    % If we reached the end of second vector's basis
                    if k == s2

                        % Check rest of first vector's basis for matching
                        % basis
                        for z=i:s1
                            if b1(z) == b2(k)
                                sign = 0;
                                return
                            end
                        end

                        %Calculate remaining flips and sign
                        flips = flips + (count * (s1-i+1));
                        sign = (-1)^flips;
                        return
                    end
                    k = k+1;
                end

                %If basis matches, component is zeroed out
                if b1(i) == b2(k)
                    sign = 0;
                    return
                end
                flips = flips + count;
            end
            flips = flips;
            sign = (-1)^flips;
        end
    end

    events
        ChangeBasisEvent
    end
end