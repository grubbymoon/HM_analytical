#include<stdio.h>
#include<string.h>
#include<math.h>
int main ( int argc, char *argv[] )
{         
     char buff[120];
     char buff1[120];
     double tol=1.0e-10;
     double E,nu, k, lambda, G; 
     double H, pinit, t;
     double M, x; 
     double u, s, p;
     double u0, s0, p0;
     const int nn=20;

     const double pi = 4.0*atan(1.0); 

     printf("File name \n" ); 
     scanf(" %s%*[^\n]%*c",buff1);
      if(argc>1)
     {  
       strcpy(buff,argv[1]);
     }
     else
     {
       strcpy(buff,buff1);
     }

     FILE *in=fopen(buff, "r");
     strcat(buff, ".out");
     FILE *out=fopen(buff, "wt");
    
     fscanf(in, " %lf %lf  %lf",&E, &nu, &k);
     fscanf(in, " %lf %lf  %lf",&H, &pinit, &t);
     
     lambda = E * nu / ((1. + nu) * (1. - 2. * nu));
     G = 0.5 * E / (1. + nu);

 
     t *= (lambda+2.0*G)*k/(H*H);

      double dH=H/nn;
      for(int i=0; i<nn+2; i++)
      {

         x=dH*i/H;
         
        
         if(i==nn+1) x=0.080814600000000;
        
         int n=0;
         double err=1.0e+10; 
         
         u0=1.0e10;
         p0=1.0e10;
         s0=1.0e10;  
         u=1.0-x;
         p=0.0;
         s=-1.0;  
            
         while(err>tol)
         {
            M=0.5*pi*(2.0*n+1); 
            s +=  2.0*sin(M*x)*exp(-M*M*t)/M;
            err = fabs(s0-s);
            s0=s;
            n++;               
         } 
         
         err=1.0e+10; 
         n=0;
         while(err>tol)
         {
            M=0.5*pi*(2.0*n+1); 
            u -=  2.0*cos(M*x)*exp(-M*M*t)/(M*M);
            err = fabs(u0-u);
            u0=u;
            n++;                        
         }
     
         err=1.0e+10;
         n=0;  
         while(err>tol)
         {
            M=0.5*pi*(2.0*n+1); 
            p += 2.0*sin(M*x)*exp(-M*M*t)/M;
            err = fabs(p0-p);
            p0=p;
            n++;
                        
         } 
         
         s *= pinit;
         u *= (pinit*H)/(lambda+2.0*G);
         p *= pinit;
         
         fprintf(out, "%12.5g  %12.5g  %12.5g  %12.5g \n", x*H, u, s, p);
//         fprintf(out, "%12.5g  %12.5g  \n", x*H, u);
 
      }            



     fclose(out);
     return 1;
}
