#include <iostream>
#include <fstream>
#include <cmath>

double fun1(double x1, double x2, double t, double gama, double beta);
double fun2(double x1, double x2, double t, double gama, double beta);
double fun3(double x1, double x2, double t, double gama, double beta);
void pasoRungeKutta4(double & x1, double & x2,double & x3, double & t, double dt, double gama, double beta);

template <typename F>
double bisection(double xl0, double xu0, double eps, F f);

const double alpha=1.038;
const double N=8e6;
const double pico=159; //dia del pico



struct Functor{
  double s, i, r, t, told, tnew, iold, inew, beta,  dt=0.01;

  double operator() (double z){

    beta=std::log(alpha)+z;   
        
    for(i=1143/N, r=45/N, s=1-i-r, t=0; t<=500;){ //colocar las condiciones iniciales

      told=t-dt;
      pasoRungeKutta4( s, i, r, told, dt, z, beta);
      iold=i; 
      
      pasoRungeKutta4( s, i, r, t, dt, z, beta);
      inew=i; tnew=t;
      
      if(iold>inew)
	{
	  return told-pico;}
    }
    // std::cout<< xf << std::endl;
    return 0;
  };
};

int main(void)
{
  std::cout.precision(8);
  std::cout.setf(std::ios::scientific);
 
  double t=0;
  double i=1143/N, r=45/N;   //condiciones iniciales
  double s=1-i-r;
  double dt=0.01;
  double gama;
  double eps=1e-08;
  double beta;
  
Functor fun;
gama= bisection(0,1,eps,fun);
beta=std::log(alpha)+gama;

std::cout << " # Gama y beta son" <<"\t" << gama <<"\t" << beta << std::endl; 
    
  //std::ofstream myfile("datos.dat");
  
    for(i=1143/N, r=45/N, s=1-i-r, t=0; t<=500; ){

      // myfile << t <<"\t"<< x2 <<std::endl;
       std::cout << t <<"\t" << s <<"\t" << i <<"\t" << r << std::endl;
       pasoRungeKutta4( s, i, r, t, dt, gama, beta);    
  }

    // myfile.close();
  
  return 0;
}

double fun1(double x1, double x2, double t, double gama, double beta) 
{
  return -beta*x2*x1;
}
double fun2(double x1, double x2, double t, double gama, double beta)
{
  return beta*x1*x2-gama*x2;
}
double fun3(double x1, double x2, double t, double gama, double beta)
{
  return gama*x2;
}
void pasoRungeKutta4(double & x1, double & x2,double & x3, double & t, double dt, double gama, double beta)
{
  double dx11, dx21, dx31, dx41;                      double dx12, dx22, dx32, dx42;                            double dx13, dx23, dx33, dx43;
  dx11=dt*fun1(x1,x2,t,gama,beta);                             dx12=dt*fun2(x1,x2,t,gama,beta);                                     dx13=dt*fun3(x1,x2,t,gama,beta);
  dx21=dt*fun1(x1+0.5*dx11,x2+0.5*dx12,t+0.5*dt,gama,beta);    dx22=dt*fun2(x1+0.5*dx11,x2+0.5*dx12,t+0.5*dt,gama,beta);            dx23=dt*fun3(x1+0.5*dx11,x2+0.5*dx12,t+0.5*dt,gama,beta);
  dx31=dt*fun1(x1+0.5*dx21,x2+0.5*dx22,t+0.5*dt,gama,beta);    dx32=dt*fun2(x1+0.5*dx21,x2+0.5*dx22,t+0.5*dt,gama,beta);            dx33=dt*fun3(x1+0.5*dx21,x2+0.5*dx22,t+0.5*dt,gama,beta);
  dx41=dt*fun1(x1+dx31,x2+dx32,t+dt,gama,beta);                dx42=dt*fun2(x1+dx31,x2+dx32,t+dt,gama,beta);                   dx43=dt*fun3(x1+dx31,x2+dx32,t+dt,gama,beta);
  
  x1+=(dx11+2*dx21+dx31+dx41)/6;                     x2+=(dx12+2*dx22+dx32+dx42)/6;                        x3+=(dx13+2*dx23+dx33+dx43)/6;
  t+=dt;
}
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
