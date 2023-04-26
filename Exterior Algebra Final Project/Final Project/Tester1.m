clear
clc

Metric4 = [
1 0 0 0;
0 1 0 0;
0 0 1 0;
0 0 0 -1;
];

%Constructing an Exterior Algebra
E = ExteriorAlgebra(4, eye(4), Metric4);
E.DispBasis();
E.DispStdMetric();

%Metric adjusted to basis depends on multivector types. First arg is grade
CoToCo_Metric = E.GetMetric(1, MultivectorType.Covector, MultivectorType.Contravector)

%Construction Multivectors
m1 = E.CreateMultivectorCurrentBasis(1,[1 1 1 1], MultivectorType.Covector);
m2 = E.CreateMultivectorCurrentBasis(2,randn(1,6), MultivectorType.Contravector);


%Display Multivector components, then their true value
m1Data = [m1.components;E.EvalMultivector(m1)]
m2Data = [m2.components E.EvalMultivector(m2)]


%Change Basis
E.ChangeBasis(randn(4))
E.DispBasis();

%Metric is also adjusted
CoToCo_Metric = E.GetMetric(1, MultivectorType.Covector, MultivectorType.Contravector)

%Display Multivector components, then their true value
%Components change after 
m1Data2 = [m1.components;E.EvalMultivector(m1)]
m2Data2 = [m2.components E.EvalMultivector(m2)]


%Show evaluated components do not change when basis changes
EvalDiff = [sum(m1Data(2,:) - m1Data2(2,:)) sum(m2Data(:,2) - m2Data2(:,2))]

