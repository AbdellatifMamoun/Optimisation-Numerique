
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Lagrangien fonction  qui calcule le minimum local pour un probleme
%d'optimisation avec contrainte d'égalité
%[x,lambda,u] = Lagrangien(f,gradif,hesif,c,gradC,hesiC,x0,lambda0,normc,jacobienC,phi)
%Sorties :
% x : le minimum local , iteration :  le nombre d'itérations effectuées
% par l'algorithme.
% lambda : multiplicateur de lagrange  
% u : multiplicateur de lagrange
%Entrees
%f:la fonction à optimiser
%x0 : le point de départ de l'algorithme 
%gradif : le gradient de f 
%hesif : la hesienne de f
%c : la fonction de contrainte
%gradC : gradient de c
%hesiC : gradient de c
%lambda0 : multiplicateur de lagrange initial
%normc : norme de c
%jacobienC : jacobienC
%phi : constante pour approximer la hesienne du Lagrangien
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [x,lambda,u] = Lagrangien(f,gradif,hesif,c,gradC,hesiC,x0,lambda0,normc,jacobienC,phi)

u0 = 25;
tol = 10^(-12);
n_0 = 0.1258925;
alpha = 0.1;
beta = 0.9;
betak=beta;
eps0 = 1/u0; 
n0=n_0/(u0^(alpha));
xk = x0;
lambdak = lambda0;
uk = u0;
epsk = eps0;
nk=n0;
k=0;
t = 2;
eps = 10^(-12);
cond = true;

while ( cond || normc(xk) > eps ) && k < 200
    
    %calcul du lagrangien augumenté du problème 
    L = @(x) (f(x)+lambdak'*c(x)+(uk/2)*(normc(x)^2));
    gradiL = @(x) (gradif(x) + gradC(x)*lambdak + uk*(jacobienC(x)')*c(x));
    hesiL =@(x) ( hesif(x) + lambdak'*hesiC(x) + uk*(jacobienC(x)')*c(x)+phi);
    
    %calcul du solution du problème sans cantrainte
    [xk1,~,~] = Region_confiance(L,xk,gradiL,hesiL,epsk,0.5,2,0.25,0.75,10000,200,200,"Cauchy");
    cond = (norm(gradif(xk1)) > tol* norm(gradif(x0)+eps));

    %si on a convergence de l'algorithme global on s'arrete 
    if (norm(gradiL(xk1)) <= tol* norm(gradiL(x0)+eps) && normc(xk1) <= eps) 
        x = xk;
        return
    end
    xk =xk1;
    %mise à jour des multiplicateurs
    if norm(c(xk))<= nk
        lambdak=lambdak + uk*c(xk);
        epsk = epsk / uk;
        nk = nk/(uk^betak);
        k=k+1;

    %mise à jour des pénalités.
    else
        %lambdak=lambdak;
        uk=t*uk;
        epsk = eps0 / uk;
        nk = n_0/(uk^alpha);
        k=k+1;
    end
    
end

x=xk;
lambda=lambdak;
u=uk;

