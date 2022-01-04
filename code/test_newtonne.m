

x011=[1 ;0 ;0];
x012=[10;3;-2.2];
x021 = [-1.2;1];
x022 = [10;0];
x023 = [0;1/200+1/(10^12)];

tol = 10^(-8);
iterMax = 8;

%[xsol_f11,iter1,flag11] = Newtonne_locale(@f,x011,@gradientf,@hesiennef,tol,iterMax);
%[xsol_f12,iter2,flag12] = Newtonne_locale(@f,x012,@gradientf,@hesiennef,tol,iterMax);
%[xsol_f21,iter3,flag21] = Newtonne_locale(@f2,x021,@gradient2,@hesienne2,tol,iterMax);
%[xsol_f22,iter4,flag22] = Newtonne_locale(@f2,x022,@gradient2,@hesienne2,tol,iterMax);
%[xsol_f23,iter5,flag23] = Newtonne_locale(@f2,x023,@gradient2,@hesienne2,tol,iterMax);



disp('------------------------------------------test de newtonne local -------------------------------------------------------------')
disp('------------------------la fonction f1  avec x0 = X011 -------------------------------------');x011
[xsol_f11,iter1,flag11] = Newtonne_locale(@f,x011,@gradientf,@hesiennef,tol,iterMax)
f_xsol = f(xsol_f11)
gradient_xsol = gradientf(xsol_f11)

disp('------------------------la fonction f1  avec x0 = X012 -------------------------------------');x012
[xsol_f12,iter2,flag12] = Newtonne_locale(@f,x012,@gradientf,@hesiennef,tol,iterMax)
f_xsol = f(xsol_f12)
gradient_xsol = gradientf(xsol_f12)

disp('------------------------la fonction f2  avec x0 = X021 -------------------------------------');x021
[xsol_f21,iter3,flag21] = Newtonne_locale(@f2,x021,@gradient2,@hesienne2,tol,iterMax)
f_xsol = f2(xsol_f21)
gradient_xsol = gradient2(xsol_f21)


disp('------------------------la fonction f2  avec x0 = X022 -------------------------------------');x022
[xsol_f22,iter4,flag22] = Newtonne_locale(@f2,x022,@gradient2,@hesienne2,tol,iterMax)
f_xsol = f2(xsol_f22)
gradient_xsol = gradient2(xsol_f22)

disp('------------------------la fonction f2  avec x0 = X023 -------------------------------------');x023
[xsol_f23,iter5,flag23] = Newtonne_locale(@f2,x023,@gradient2,@hesienne2,tol,iterMax)
f_xsol = f2(xsol_f23)
gradient_xsol = gradient2(xsol_f23)












