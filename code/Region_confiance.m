%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Region_confiance fonction  qui calcule le minimum local d'une fonction 
% par l'algorithme de region de confiance	
%[xsol,iteration] = Region_confiance(f,x0,gradi,hesi,epsk,gamma1,gamma2,neta1,neta2,deltaMax,Delta0,iterM,option)
%Sorties :
% xsol : le minimum local , iteration :  le nombre d'itérations effectuées
% par l'algorithme.
% iteration : le nombre d'itération 
%flag : indicateur de déroulement de l'algorithme :
        %-0 : CN de convergence de f
        %-1 : stagnation des itérés , faible variation de xk
        %-2: stagnation , faible variation de f(xk)
        %-3: nombre d'itération max atteint 
%Entrees
%x0 : le point de départ de l'algorithme 
%gradi : le gradient de f 
%hesi : la hesienne de f
%epsk : la tolérance pour la convergence du gradient 
%gamma1 et gamma 2 : deux constantes tel que 0 < gamma1 < gamma2 < 1
%neta1 et neta2 : deux constantes tels que  0 < neta1 < neta2 < 1
%deltaMax : région de confiance maximale
%Delta0 : région de confiance initiale
%iterM : nombre d'itération maximale 
%option : option pour la region de confiance , "Cauchy" pour le pas de
%cauchy ou "Gradient pour le gradient conjugué      

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xsol,iteration,flag] = Region_confiance(f,x0,gradi,hesi,epsk,gamma1,gamma2,neta1,neta2,deltaMax,Delta0,iterM,option)

%les parametres de l'algortihme
delta_max = deltaMax;
delta0 = Delta0;
y1 = gamma1;
y2=gamma2;
n1 = neta1;
n2=neta2;
tol= epsk;
eps = epsk;
xk=x0;
deltak=delta0;

iterMax = iterM;

xk1=xk;%Xk+1

nostop = true; % condition d'arret de l'algorithme on s'arrete si nostop = false 

iter=0;

%on détermine si l'option choisie est celle du pas de cauchy ou gradient conjugé
 if (option == "Cauchy")
       op="Cauchy";
 else
        if(option == "Gradient")
         op = "Gradient";
        else
            disp("erreur de choix d'option")
            return
        end
 end
    
while (nostop)
    xk =xk1;
    if (op == "Cauchy")
        
        %calcul de la direction de descente par le pas de cauchy
        sk=PasDeCauchy(gradi(xk),hesi(xk),deltak);
        
    else
     
         %calcul de la direction de descente par gradient conjugué
         sk=Grad_Conj(f,xk,deltak,gradi(xk),hesi(xk));
        
    end
    
    
    phiK = (f(xk)-f(xk+sk))/(-gradi(xk).'*sk+0.5*sk.'*hesi(xk)*sk);
    
    %mise à jour de l'itéré courant 
    if phiK >= n1
        xk1=xk+sk;
        
        flag1 = (norm(xk1-xk))< tol*norm(xk+sqrt(eps));
        flag2 = (norm(f(xk1)-f(xk))< tol*norm(f(xk)+eps));
            
     
    else
        
        xk1=xk;
        flag1=false;
        flag2=false;
    end
    
    %mise à jour de la région de confiance
    if phiK >= n2
        deltak=min(y2*deltak,delta_max);
    elseif phiK <= n2 && phiK >= n1 
        deltak=deltak;
    else
        deltak = y1*deltak;
    end
    iter = iter+1;
    flag0 =  (norm(gradi(xk)) < tol* norm(gradi(x0)+eps));
    flag3 = iter >= iterMax;
    
    %vérifier la condition d'arret
    nostop = ((~flag0) && (~flag1) && (~flag2) && (~flag3));
end
xsol=xk1;
iteration=iter;

%affectation du flag
if (flag0) 
    flag =0;
elseif(flag1)
    flag =1;
elseif(flag2)
    flag=2;
else
    flag=3;
end

