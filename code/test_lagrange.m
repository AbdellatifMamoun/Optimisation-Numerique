
%les points de tests
xc11 = [0;1;1];
xc12 = [0.5;1.25;1];
xc22 = [sqrt(3)/2;sqrt(3)/2];
xc21=[1;0];

%les paramètres d'entree
phi = 0;
lambda0 = 0.75;
mu0 = 0.75;






disp('------------------------------------------test du Lagrangien pour les problèmes d poptimisations avec contraintes  -------------------------------------------------------------')
disp('------------------------la fonction f1  avec x0 = Xc11 -------------------------------------');xc11
[xsol_f11,lambda1,u1] = Lagrangien(@f,@gradientf,@hesiennef,@c,@gradC,@hesiC,xc11,lambda0,@normc,@jacobienneC,phi)
f_xsol = f(xsol_f11)
gradient_xsol = gradientf(xsol_f11)

%disp('------------------------la fonction f1  avec x0 = Xc12  -------------------------------------');xc12
%[xsol_f12,lambda2,u2] =  Lagrangien(@f,@gradientf,@hesiennef,@c,@gradC,@hesiC,xc12,lambda0,@normc,@jacobienneC,phi);
%f_xsol = f(xsol_f12)
%gradient_xsol = gradientf(xsol_f12)

disp('------------------------la fonction f2  avec x0 = Xc21  -------------------------------------');xc21
[xsol_f21,lambda3,u3] = Lagrangien(@f2,@gradient2,@hesienne2,@C2,@gradientC2,@hesienneC2,xc21,lambda0,@normeC2,@jacobienC2,phi)
f_xsol = f2(xsol_f21)
gradient_xsol = gradient2(xsol_f21)


disp('------------------------la fonction f2  avec x0 = Xc22-------------------------------------');xc22
[xsol_f22,lambda4,u4] = Lagrangien(@f2,@gradient2,@hesienne2,@C2,@gradientC2,@hesienneC2,xc22,lambda0,@normeC2,@jacobienC2,phi)
f_xsol = f2(xsol_f22)
gradient_xsol = gradient2(xsol_f22)


