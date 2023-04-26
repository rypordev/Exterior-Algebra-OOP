clear
clc
dim = 3;
m = randn(dim);
m = m*m';
k = 1;


E = ExteriorAlgebra(dim, eye(dim), m);
a = E.CreateMultivectorStandardBasis(1,randn([1,3]), MultivectorType.Covector);
b = E.CreateMultivectorStandardBasis(1,randn([1,3]), MultivectorType.Covector);
ab = E.Wedge(a,b);
E.ChangeBasis([3,1,2;1,2,4;3,1,1]);
abd = E.Wedge(a,b);
ab.components
abd.components

E.GetMetric(1,MultivectorType.Contravector,MultivectorType.Contravector);
u = E.CreateMultivectorCurrentBasis(1,randn([1 3]), MultivectorType.Covector);
v = E.CreateMultivectorCurrentBasis(1,randn([1 3]), MultivectorType.Covector);

bw = E.Wedge(u,v);
w = E.HodgeStar(bw);

ww = E.InnerProduct(w,w) * det(m)
s = E.InnerProduct(u,u)*E.InnerProduct(v,v)-E.InnerProduct(u,v)^2