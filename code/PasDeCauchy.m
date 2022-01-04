%Calcul du pas de Cauchy d'une fonction pour un point
%[x] = PasDeCauchy(gradi,hesi,D)
%sortie : x le pas de Cauchy
%Entrees : 
%gradi : le gradient de f en un point 
%hesi : la hesienne de f en un point
%D : la région de confiance

function [x] = PasDeCauchy(gradi,hesi,D)
delta = D/norm(gradi);

m=norm(gradi)^2/((gradi.')*hesi*gradi);

if gradi.'*hesi*gradi >0 
    if  m > 0 && m <= delta 
        x = -m*gradi; 
    elseif m > delta 
        x=-delta*gradi;
    elseif m < 0
        x= -0.0000000000001*gradi;
    end
else
    x = -delta*gradi ;
end
  
end

