
x011=[1 ;0 ;0];
x012=[10;3;-2.2];
x021 = [-1.2;1];
x022 = [10;0];
x023 = [0;1/200+1/(10^12)];


delta_max = 10000;
delta0 = 200;
y1 = 0.5;
y2=2;
n1 = 0.25;
n2=0.75;
tol = 10^(-9);
iterMax = 70;

%[s1,iter1]= Region_confiance(@f,x011,@gradientf,@hesiennef,tol,y1,y2,n1,n2,delta_max,delta0,iterMax,"Cauchy");
%[s2,iter2] = Region_confiance(@f,x012,@gradientf,@hesiennef,tol,y1,y2,n1,n2,delta_max,delta0,iterMax,"Cauchy");
%[s3,iter3] = Region_confiance(@f2,x021,@gradient2,@hesienne2,tol,y1,y2,n1,n2,delta_max,delta0,iterMax,"Cauchy");
%[s4,iter4] = Region_confiance(@f2,x022,@gradient2,@hesienne2,tol,y1,y2,n1,n2,delta_max,delta0,iterMax,"Cauchy");
%[s5,iter5] = Region_confiance(@f2,x023,@gradient2,@hesienne2,tol,y1,y2,n1,n2,delta_max,delta0,iterMax,"Cauchy");


%[x1,k1]= Region_confiance(@f,x011,@gradientf,@hesiennef,tol,y1,y2,n1,n2,delta_max,delta0,iterMax,"Gradient");
%[x2,k2]= Region_confiance(@f,x012,@gradientf,@hesiennef,tol,y1,y2,n1,n2,delta_max,delta0,iterMax,"Gradient");
%[x3,k3]= Region_confiance(@f2,x021,@gradient2,@hesienne2,tol,y1,y2,n1,n2,delta_max,delta0,iterMax,"Gradient");
%[x4,k4]= Region_confiance(@f2,x022,@gradient2,@hesienne2,tol,y1,y2,n1,n2,delta_max,delta0,iterMax,"Gradient");
%[x5,k5]= Region_confiance(@f2,x023,@gradient2,@hesienne2,tol,y1,y2,n1,n2,delta_max,delta0,iterMax,"Gradient");


disp('------------------------------------------test de la région de confinace   -------------------------------------------------------------')
disp('région de confiance avec pas de cauchy ')
disp('------------------------la fonction f1  avec x0 = X011  et pas de cauchy-------------------------------------');x011
[xsol_f11,iter1,flag11] = Region_confiance(@f,x011,@gradientf,@hesiennef,tol,y1,y2,n1,n2,delta_max,delta0,iterMax,"Cauchy")
f_xsol = f(xsol_f11)
gradient_xsol = gradientf(xsol_f11)

disp('------------------------la fonction f1  avec x0 = X012 et pas de cauchy -------------------------------------');x012
[xsol_f12,iter2,flag12] =  Region_confiance(@f,x012,@gradientf,@hesiennef,tol,y1,y2,n1,n2,delta_max,delta0,iterMax,"Cauchy")
f_xsol = f(xsol_f12)
gradient_xsol = gradientf(xsol_f12)

disp('------------------------la fonction f2  avec x0 = X021 et pas de cauchy -------------------------------------');x021
[xsol_f21,iter3,flag21] = Region_confiance(@f2,x022,@gradient2,@hesienne2,tol,y1,y2,n1,n2,delta_max,delta0,iterMax,"Cauchy")
f_xsol = f2(xsol_f21)
gradient_xsol = gradient2(xsol_f21)


disp('------------------------la fonction f2  avec x0 = X022 et pas de cauchy-------------------------------------');x022
[xsol_f22,iter4,flag22] = Region_confiance(@f2,x022,@gradient2,@hesienne2,tol,y1,y2,n1,n2,delta_max,delta0,iterMax,"Cauchy")
f_xsol = f2(xsol_f22)
gradient_xsol = gradient2(xsol_f22)

disp('------------------------la fonction f2  avec x0 = X023 et pas de cauchy -------------------------------------');x023
[xsol_f23,iter5,flag23] =Region_confiance(@f2,x023,@gradient2,@hesienne2,tol,y1,y2,n1,n2,delta_max,delta0,iterMax,"Cauchy")
f_xsol = f2(xsol_f23)
gradient_xsol = gradient2(xsol_f23)

disp('==============================================================================================================================')

disp('région de confiance avec gradient conjugué ')
disp('------------------------la fonction f1  avec x0 = X011  et gradient conjugué-------------------------------------');x011
[xsol_f11_conj,iter1_conj,flag11_conj] = Region_confiance(@f,x011,@gradientf,@hesiennef,tol,y1,y2,n1,n2,delta_max,delta0,iterMax,"Gradient")
f_xsol_conj = f(xsol_f11_conj)
gradient_xsol_conj = gradientf(xsol_f11_conj)

disp('------------------------la fonction f1  avec x0 = X012 et gradient conjugué -------------------------------------');x012
[xsol_f12_conj,iter2_conj,flag12_conj] =  Region_confiance(@f,x012,@gradientf,@hesiennef,tol,y1,y2,n1,n2,delta_max,delta0,iterMax,"Gradient")
f_xsol_conj = f(xsol_f12_conj)
gradient_xsol_conj = gradientf(xsol_f12_conj)

disp('------------------------la fonction f2  avec x0 = X021 et gradient conjugué -------------------------------------');x021
[xsol_f21_conj,iter3_conj,flag21_conj] = Region_confiance(@f2,x022,@gradient2,@hesienne2,tol,y1,y2,n1,n2,delta_max,delta0,iterMax,"Gradient")
f_xsol_conj = f2(xsol_f21_conj)
gradient_xsol_conj = gradient2(xsol_f21_conj)


disp('------------------------la fonction f2  avec x0 = X022 et gradient conjugué-------------------------------------');x022
[xsol_f22_conj,iter4_conj,flag22_conj] = Region_confiance(@f2,x022,@gradient2,@hesienne2,tol,y1,y2,n1,n2,delta_max,delta0,iterMax,"Gradient")
f_xsol_conj = f2(xsol_f22_conj)
gradient_xsol_conj = gradient2(xsol_f22_conj)

disp('------------------------la fonction f2  avec x0 = X023 et gradient conjugué -------------------------------------');x023
[xsol_f23_conj,iter5_conj,flag23_conj] =Region_confiance(@f2,x023,@gradient2,@hesienne2,tol,y1,y2,n1,n2,delta_max,delta0,iterMax,"Gradient")
f_xsol_conj = f2(xsol_f23_conj)
gradient_xsol_conj = gradient2(xsol_f23_conj)




