Metric4 = [
1 0 0 0;
0 1 0 0;
0 0 1 0;
0 0 0 -1;
];

%Constructing an Exterior Algebra
E = ExteriorAlgebra(4, eye(4), Metric4);
E.DispBasis();


%Simple Wedge Product Test
m11 = E.CreateMultivectorCurrentBasis(1,[1 0 0 0], MultivectorType.Contravector);
m12 = E.CreateMultivectorCurrentBasis(1,[0 0 0 1], MultivectorType.Contravector);
m13 = E.Wedge(m11,m12)
[m13.components E.EvalMultivector(m13)]

m21 = E.CreateMultivectorCurrentBasis(1,randn(1,4), MultivectorType.Covector);
m22 = E.CreateMultivectorCurrentBasis(1,randn(1,4), MultivectorType.Covector);
m23 = E.Wedge(m21,m22)
[m23.components;E.EvalMultivector(m23)]


E.ChangeBasis(randn(4));
E.DispBasis();
[m13.components E.EvalMultivector(m13)]
[m23.components;E.EvalMultivector(m23)]