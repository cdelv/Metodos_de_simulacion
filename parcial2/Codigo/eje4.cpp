
#include <iostream>
#include <fstream>
#include <cmath>
#include "Random64.h"
using namespace std;

const int Lx=256;
const int Ly=64;

const int Q=9;

const double tau=0.55;
const double Utau=1.0/tau;
const double UmUtau=1-Utau;

//--------------------- Clase LatticeBoltzmann ------------
class LatticeBoltzmann{
private:
  double w[Q];
  int V[2][Q]; //V[0][i]=Vi_x ,  V[1][i]=Vi_y
  double f[Lx][Ly][Q], fnew[Lx][Ly][Q]; // f[ix][iy][i]
public:
  LatticeBoltzmann(void);
  double rho(int ix,int iy,bool UseNew);
  double Jx(int ix,int iy,bool UseNew);
  double Jy(int ix,int iy,bool UseNew);
  double feq(double rho0,double Ux0,double Uy0,int i);
  void Colisione(void);
  void Adveccione(void);
  void Inicie(void);
  void Imprimase(const char * NombreArchivo);
  double Varianza(void);
  double Detector(void);
};  
LatticeBoltzmann::LatticeBoltzmann(void){
  //cargar los pesos
  w[0]=4.0/9; w[1]=w[2]=w[3]=w[4]=1.0/9; w[5]=w[6]=w[7]=w[8]=1.0/36;
  //cargar los vectores velocidad
  V[0][0]=0;  V[0][1]=1;  V[0][2]=0;  V[0][3]=-1; V[0][4]=0;
  V[1][0]=0;  V[1][1]=0;  V[1][2]=1;  V[1][3]=0;  V[1][4]=-1;

              V[0][5]=1;  V[0][6]=-1; V[0][7]=-1; V[0][8]=1;
              V[1][5]=1;  V[1][6]=1;  V[1][7]=-1; V[1][8]=-1;
}
double LatticeBoltzmann::rho(int ix,int iy,bool UseNew){
  double suma; int i;
  for(suma=0,i=0;i<Q;i++)
    if(UseNew) suma+=fnew[ix][iy][i]; else suma+=f[ix][iy][i];
  return suma;
}  
double LatticeBoltzmann::Jx(int ix,int iy,bool UseNew){
  
  double deltaP=0.001, eta=10, Jx;

  Jx=0.5*rho(ix,iy,false)*(deltaP/eta)*iy*(Ly-1-iy);

  return Jx;
}  
double LatticeBoltzmann::Jy(int ix,int iy,bool UseNew){
  
  return 0;
}  
double  LatticeBoltzmann::feq(double rho0,double Ux0,double Uy0,int i){
  double UdotVi=Ux0*V[0][i]+Uy0*V[1][i]; double U2=Ux0*Ux0+Uy0*Uy0;
  return w[i]*rho0*(1+3*UdotVi+4.5*UdotVi*UdotVi-1.5*U2);
}  
void LatticeBoltzmann::Colisione(void){
  int ix,iy,i; double rho0,Ux0,Uy0;
  for(ix=0;ix<Lx;ix++) //para cada celda
    for(iy=0;iy<Ly;iy++){
      //calcular los campos macroscópicos en la celda
      rho0=rho(ix,iy,false); Ux0=Jx(ix,iy,false)/rho0; Uy0=Jy(ix,iy,false)/rho0;
      for(i=0;i<Q;i++) //en cada dirección
	fnew[ix][iy][i]=UmUtau*f[ix][iy][i]+Utau*feq(rho0,Ux0,Uy0,i);
    }  
}

void LatticeBoltzmann::Adveccione(void){
  for(int ix=0;ix<Lx;ix++) //para cada celda
    for(int iy=0;iy<Ly;iy++)
      for(int i=0;i<Q;i++) //en cada dirección
	f[(ix+V[0][i]+Lx)%Lx][(iy+V[1][i]+Ly)%Ly][i]=fnew[ix][iy][i]; //fronteras periódicas
}
void LatticeBoltzmann::Inicie(void){
  double sigma=Ly/9; double mu=Lx/8; double escala=1/(sigma*sqrt(2*M_PI));
  double rho=0;
  for(int ix=0;ix<Lx;ix++) //para cada celda
    for(int iy=0;iy<Ly;iy++)
      for(int i=0;i<Q;i++){ //en cada dirección
	rho=escala*exp(-0.5*pow((ix-mu)/sigma,2));
	f[ix][iy][i]=feq(rho,0,0,i); //Ux=Uy=0
      }
}  
void LatticeBoltzmann::Imprimase(const char * NombreArchivo){
  ofstream MiArchivo(NombreArchivo); double rho0;
  for(int ix=0;ix<Lx;ix+=4){
    for(int iy=0;iy<Ly;iy+=4){
      rho0=rho(ix,iy,true); 
      MiArchivo<<ix<<"\t"<<iy<<"\t"<<rho0<<endl;
    }
    MiArchivo<<endl;
  }
  MiArchivo.close();
}
double LatticeBoltzmann::Varianza(void)
{
  int ix, iy;
  double N, ixprom, var, aux;

  //calculo de N
  for(ix=0; ix<Lx; ix++)
    for(iy=0; iy<Ly; iy++)
      aux+=rho(ix,iy,true);
  N=aux;

  //calculo de ixprom
  for(ix=0; ix<Lx; ix++)
    for(iy=0; iy<Ly; iy++)
      aux+=rho(ix,iy,true)*ix;
  ixprom=aux/N;

  //calculo de la varianza
  for(ix=0; ix<Lx; ix++)
    for(iy=0; iy<Ly; iy++)
      aux+=rho(ix,iy,true)*pow(ix-ixprom,2);
  var=aux/N;

  return var;
}
double LatticeBoltzmann::Detector(void)
{
  double Detector, aux;

  for(int iy=0; iy<Ly-1; iy++)
    aux+=rho(150,iy,true);
  Detector=aux;
  return Detector;
}
//------------------- Funciones principal ------------

int main(void){
  LatticeBoltzmann Aire;
  int t,tmax=5000;
  
  //Inicie
  Aire.Inicie();
  //Corra
  ofstream MiArchivo1("eje4_detector.dat");
  
  for(t=0;t<tmax;t++){
    Aire.Colisione();
    Aire.Adveccione();
    MiArchivo1 << t<<"\t"<<Aire.Detector()<<std::endl;
  }
  
  MiArchivo1.close();
  Aire.Imprimase("eje4_5000.dat");
   
  return 0;
}  
