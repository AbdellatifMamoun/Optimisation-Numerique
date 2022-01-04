function [x] = Grad_Conj(f,xk,deltak,gradf,hesf)
%RC Summary of this function goes here
%   Detailed explanation goes here

 sj=0;
 s=0;
 sigmaj = 0;
 epsilon = 0.001;
 
 gj=gradf;
 
 pj=-gradf;
 
 sortir = true ;
 
 kj = (pj.').*hesf.*pj;

 betaj=0;
 iter =0;
 while (sortir) 
      iter=iter+1;
      kj = (pj.')*hesf*pj;
      
      if kj <= 0 
          
          x1=-(norm(sj)+deltak)/norm(pj);

          x2=(-norm(sj)+deltak)/norm(pj);
          
          q1 = f(xk)+(gj.')*(sj+x1*pj)+0.5*((sj+x1*pj).')*hesf*(sj+x1*pj) ;
          q2 = f(xk)+(gj.')*(sj+x2*pj)+0.5*((sj+x2*pj).')*hesf*(sj+x2*pj) ;
          
          if q1 <= q2
              sigmaj=x1;
          else
              sigmaj=x2;
          end
          
         s=sj+sigmaj*pj;
         x=s;
         return 
         

          
      end
      
      alphaj = (gj.')*gj/kj;
      if norm(sj+alphaj*pj) >= deltak
          
          sigmaj = (-norm(sj)+deltak)/norm(pj);
          
          s=sj+sigmaj*pj;
          x=s;
          return
         
      end
      
       sj = sj+alphaj*pj;
       gj1 = gj + alphaj*hesf*pj;
       betaj = (norm(gj1) / norm(gj))^2;
       gj=gj1;
       pj = -gj1+betaj*pj;
      
       if norm(gj) < norm(gradf) * epsilon
           s=sj;
           sortir = false;
       end
       
       s=sj;
     
 end

 x=s;


end

