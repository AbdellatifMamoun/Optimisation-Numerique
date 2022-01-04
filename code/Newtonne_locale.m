%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Newtonne_locale fonction  qui calcule le minimum local d'une fonction 
% par l'algorithme de Newtonne local
%[xsol,iter,flag] = Newtonne_locale(f,x0,gradientf,hesiennef,tolerance,itermax)
%Sorties
%iter : le nombre d'iteration atteint
%flag : indicateur de déroulement de l'algorithme :
        %-0 : CN de convergence de f
        %-1 : stagnation des itérés , faible variation de xk
        %-2: stagnation , faible variation de f(xk)
        %-3: nombre d'itération max atteint 
% xsol : le minimum local , iteration :  le nombre d'itérations effectuées
% par l'algorithme.
%Entrees
% f : la fonction à optimiser                          
%x0 : le point de départ de l'algorithme 
%gradientf : le gradient de f 
%hesiennef : la hesienne de f 
%tolerance : la tolerance de convergence du gradient
%itermax : le nombre d'itération maximale
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xsol,iter,flag] = Newtonne_locale(f,x0,gradientf,hesiennef,tolerance,itermax)

tol= tolerance;
eps = 0.000001;
xk=x0;
k=0;
iterMax=itermax;
nostop = true;
while (nostop) % && norm(xk1-xk) > tol * (norm(x0)+eps)
     k=k+1;
    dk1 = hesiennef(xk)\(-gradientf(xk));
    xk1=xk+dk1;
   
    flag1 = (norm(xk1-xk))< tol*norm(xk+sqrt(eps));
    flag2 = (norm(f(xk1)-f(xk))< tol*norm(f(xk)+eps));
    flag0 =  (norm(gradientf(xk)) < tol* norm(gradientf(x0)+eps));
    flag3 = k > iterMax;
    nostop = ((~flag0) && (~flag1) && (~flag2) && (~flag3));
    xk=xk1;
    
   
end
xsol=xk;
iter=k-1;
if (flag0) 
    flag =0;
elseif(flag1)
    flag =1;
elseif(flag2)
    flag=2;
else
    flag=3;
end
    


