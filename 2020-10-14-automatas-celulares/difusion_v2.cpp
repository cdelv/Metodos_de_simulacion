

#include <cmath>
#include "Random64.h"
using namespace std;

const int Lx=256;
const double p=0.5;

//--------------------- Clase LatticeGas ------------
class LatticeGas{
private:
  int V[2]; //V[i], i=0 (derecha), i=1 (izquierda)
  int n[Lx][2];  int nnew[Lx][2];
public:
  LatticeGas(void);
  void Inicie(int N,double mu,double sigma,Crandom & ran64);
  void Show(bool ShowNew);
  void Colisione(Crandom & ran64);
  void Adveccione(void);
  double GetSigma2(void);
};
LatticeGas::LatticeGas(void){
  V[0]=1;  V[1]=-1;
  for(int ix=0;ix<Lx;ix++)
    n[ix][0]=n[ix][1]=nnew[ix][0]=nnew[ix][1]=0;
}
void LatticeGas::Inicie(int N,double mu,double sigma,Crandom & ran64){
  int ix,i;
  do{
    ix=(int) ran64.gauss(mu,sigma); if(ix<0) ix=0; if(ix>(Lx-1)) ix=Lx-1;
    i=(int) 2*ran64.r();
    if(n[ix][i]==0)
      {n[ix][i]=1; N--;}
    }while(N>0);
}  
void LatticeGas::Show(bool ShowNew){
  for(int i=0;i<2;i++){
    for(int ix=0;ix<Lx;ix++)
      if(ShowNew) cout<<nnew[ix][i]; else cout<<n[ix][i];
    cout<<endl;
  }
  cout<<endl;
}  
void LatticeGas::Colisione(Crandom & ran64){
  for(int ix=0;ix<Lx;ix++) //para cada celda
    if (ran64.r()<p) //con probabilidad p
      {nnew[ix][0]=n[ix][0]; nnew[ix][1]=n[ix][1];}//dejelo igual;
    else  //con probabilidad 1-p
      {nnew[ix][0]=n[ix][1]; nnew[ix][1]=n[ix][0];}//intercámbielos;      
}
void LatticeGas::Adveccione(void){
  for(int ix=0;ix<Lx;ix++) //para cada celda
    for(int i=0;i<2;i++)
      n[(ix+V[i]+Lx)%Lx][i]=nnew[ix][i]; //fronteras periódicas
}
double LatticeGas::GetSigma2(void){
  double ixprom,suma,N; int ix;
  //Calcular N
  for(N=0,ix=0;ix<Lx;ix++)
    N+=(n[ix][0]+n[ix][1]);
  //Calcular ixprom
  for(suma=0,ix=0;ix<Lx;ix++)
    suma+=(n[ix][0]+n[ix][1])*ix;
  ixprom=suma/N;
  //Calcular sigma2
  for(suma=0,ix=0;ix<Lx;ix++)
    suma+=(n[ix][0]+n[ix][1])*pow(ix-ixprom,2);
  return suma/(N-1);
}  
//------------------- Funciones Globales ------------

int main(void){
  LatticeGas Difusion;
  Crandom ran64(2);
  int t,tmax=10000;
  int N=400; double mu=Lx/2,sigma=Lx/8;

  //Inicie
  Difusion.Inicie(N,mu,sigma,ran64);
  //Corra
  for(t=0;t<tmax;t++){
    Difusion.Colisione(ran64);
    Difusion.Adveccione();
    cout<<t<<" "<<Difusion.GetSigma2()<<endl;
  }
  
  return 0;
}  
