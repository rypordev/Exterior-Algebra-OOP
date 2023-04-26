clear
clc
Metric4 = [
1 0 0 0;
0 1 0 0;
0 0 1 0;
0 0 0 -1;
];

%E = ExteriorAlgebra(4, eye(4), Metric4);
dim = 4;
m = randn(dim);
m = m*m';
E = ExteriorAlgebra(dim, eye(dim), m);

%Use Covector for nicer looking output
mv1 = E.CreateMultivectorStandardBasis(1,ones(1,4), MultivectorType.Covector);
mv2 = E.CreateMultivectorStandardBasis(2,ones(1,6), MultivectorType.Covector);
mv3 = E.CreateMultivectorStandardBasis(3,ones(1,4), MultivectorType.Covector);

star1 = E.HodgeStar(mv1);
star2 = E.HodgeStar(mv2);
star3 = E.HodgeStar(mv3);

disp("Hodge Star in Lorentz, grade 1")
disp(star1.components)

disp("Hodge Star in Lorentz, grade 2")
disp(star2.components)

disp("Hodge Star in Lorentz, grade 3")
disp(star3.components)

%Change basis changes components after hodge star
E.ChangeBasis(eye(4)*2)
E.DispBasis();

%Display new Covector metrics for clarity
E.GetMetric(1,MultivectorType.Covector,MultivectorType.Covector)
E.GetMetric(2,MultivectorType.Covector,MultivectorType.Covector)
E.GetMetric(3,MultivectorType.Covector,MultivectorType.Covector)

star1 = E.HodgeStar(mv1);
star2 = E.HodgeStar(mv2);
star3 = E.HodgeStar(mv3);

disp("Hodge Star in CoB Lorentz, grade 1")
disp(star1.components)

disp("Hodge Star in CoB Lorentz, grade 2")
disp(star2.components)

disp("Hodge Star in CoB Lorentz, grade 3")
disp(star3.components)
