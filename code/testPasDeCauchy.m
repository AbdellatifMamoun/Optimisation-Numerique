g1=[0;0];
H1=[7 0;0 2];
g2=[6;2];
H2=[7 0;0 2];

g3=[-2;1];
H3=[-2 0;0 10];


disp('===================test du pas de Cauchy===========================')
disp('============= avec g1 et H1 ====================')
x1 = PasDeCauchy(g1,H1,6)


disp('============= avec g2 et H2 ====================')
x2 = PasDeCauchy(g2,H2,6)


disp('============= avec g3 et H3 ====================')
x3 = PasDeCauchy(g3,H3,6)