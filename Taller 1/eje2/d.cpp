#include <iostream>
#include <fstream>
#include <cmath>



double fun(double t, double alpha, double x);
double fun1(double x1, double x2, double t, double L);
double fun2(double x1, double x2, double t, double L);
void pasoRungeKutta4(double & x1, double & x2, double & t, double dt, double L);
double simpson (double a, double b, int n, double alpha, double x);
double bessel (double a, double b, double n, double alpha, double x);
template <typename F>
double bisection(double xl0, double xu0, double eps, F f);





//---------------------------- functores ----------------------------------


//----------------------------Functor de Runge-Kutta----------------------
struct Functor{
  double x0=0;
  double xf=0;
  double t=0;
  double dt=0;
 
  double operator() (double z){
   
    for(t=0.01, x0=0, xf=1; t<=1;){ //colocar las condiciones iniciales
      pasoRungeKutta4(x0, xf, t, dt, z);    
    }
    // std::cout<< xf << std::endl;
    return xf;
  };
};
//-----------------------------Functor de Bessel-------------------------
struct Functor1{

  double a=0, b=M_PI, nmedios=2000;   //parametros de la función de Bessel 
  double alpha=0, x=0;
 
  double operator() (double z){
   
    return bessel(a,b,nmedios, alpha, z);
  };
};


//-----------------------------funcion principal -----------------------------------

int main(void)
{
  std::cout.precision(8);
  std::cout.setf(std::ios::scientific);

//x2 es R y x1 es la derivada de R  
  double t=0;
  double x1=0, x2=1; //condiciones iniciales
  double dt=0.001;
  double L=1;
  double eps=1e-10;
  
  double L1=0, L2=0, L3=0, L4=0, L5=0;


  Functor fun;
  fun.dt=dt;
  
  L1= bisection(2,4,eps,fun);
  L2= bisection(6,8,eps,fun);   //Lambda encontrado con Runge-Kutta
  L3= bisection(9,11,eps,fun);
  L4= bisection(14,16,eps,fun);
  
  Functor1 Bessel;
  
  L5= bisection(4,8,eps,Bessel);
    

  std::cout << L5 << std::endl;

  
  return 0;
}

//----------------------- ecuacion diferencial ----------------------------

double fun1(double x1, double x2, double t, double L) 
{
  return (-x1/t-L*L*x2);
}
double fun2(double x1, double x2, double t, double L)
{
  return x1;
}
void pasoRungeKutta4(double & x1, double & x2, double & t, double dt, double L)
{
  double dx11, dx21, dx31, dx41;                      double dx12, dx22, dx32, dx42;
  dx11=dt*fun1(x1,x2,t,L);                             dx12=dt*fun2(x1,x2,t,L);
  dx21=dt*fun1(x1+0.5*dx11,x2+0.5*dx12,t+0.5*dt,L);    dx22=dt*fun2(x1+0.5*dx11,x2+0.5*dx12,t+0.5*dt,L);
  dx31=dt*fun1(x1+0.5*dx21,x2+0.5*dx22,t+0.5*dt,L);    dx32=dt*fun2(x1+0.5*dx21,x2+0.5*dx22,t+0.5*dt,L);
  dx41=dt*fun1(x1+dx31,x2+dx32,t+dt,L);                dx42=dt*fun2(x1+dx31,x2+dx32,t+dt,L);

  x1+=(dx11+2*dx21+dx31+dx41)/6;                       x2+=(dx12+2*dx22+dx32+dx42)/6;
  
  t+=dt;
}

//------------------------encontrar las raices-------------------------------------
template <typename F>
double bisection(double xl0, double xu0, double eps, F f)
{
  int iter = 1;
  double mult=f(xl0)*f(xu0);
  double xl=xl0;
  double xu=xu0;
  
  while (mult < 0){
    double xr = (xl + xu)/2;
    //std::cout << iter << "\t"<<xl<<"\t" <<xu<<"\t" <<xr<<"\t"<< f(xr)<<"\n";
    if(std::fabs(f(xr))<=eps){
      return xr;
    }
    if(f(xr)*f(xu)<0){
      xl=xr;
    } else {
      xu=xr;
    }
    mult=f(xu)*f(xl);
    iter++;
  }
  return 0.5*(xl+xu);
}


//----------------------------- funcion de Bessel------------------------------ 


//funcion para integrar
double fun(double t, double alpha, double x)
{
  return (std::cos(alpha*t-x*(std::sin(t))));
}
//integracion
double simpson (double a, double b, int n, double alpha, double x)
{
  n*=2;
  double h=(b-a)/n;
  double suma=0, t=0;
  
  for(int ii=0; ii<=n; ii++){
    t=a+h*ii;
    
    if(ii==0 || ii==n){
      suma+=fun(t,alpha,x);
    }
    else{ if(ii%2==0){suma+=2*fun(t,alpha,x);}
      else{suma+=4*fun(t,alpha,x);}
    }
  }
    return suma*h/3;
  }

double bessel (double a, double b, double n, double alpha, double x)
{
  return simpson(a, b, n, alpha, x)*1.0/M_PI;
}
