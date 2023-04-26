clear
clc
Metric4 = [
1 0 0 0;
0 1 0 0;
0 0 1 0;
0 0 0 -1;
];

E = ExteriorAlgebra(4)

m1 = E.CreateMultivectorCurrentBasis(2,randn(1,6), MultivectorType.Contravector);
m2 = E.CreateMultivectorCurrentBasis(2,randn(1,6), MultivectorType.Covector);


MVEval11 = E.EvalMultivector(m1);
MVEval21 = E.EvalMultivector(m2);

%Show inner products for contra and co multivectors
a = [E.InnerProduct(m1,m1) E.InnerProduct(m1,m2) E.InnerProduct(m2,m1) E.InnerProduct(m2,m2)]

NewBasis = randn(4);
E.ChangeBasis(NewBasis);
E.DispBasis();

MVEval12 = E.EvalMultivector(m1);
MVEval22 = E.EvalMultivector(m2);

%Show inner products do not change with respect to basis
b = [E.InnerProduct(m1,m1) E.InnerProduct(m1,m2) E.InnerProduct(m2,m1) E.InnerProduct(m2,m2)]

Results = [sum(MVEval11-MVEval12) sum(MVEval21-MVEval22) a-b]

