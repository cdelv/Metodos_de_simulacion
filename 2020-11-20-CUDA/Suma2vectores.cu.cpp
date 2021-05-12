//Mi Primer Programa en CUDA
//scp MiPrimerProgramaEnCuda.cu cdelv@168.176.8.34:CUDA2020
//nvcc -arch=sm_30  para darle un nombre -o ArchivoCUDA
//time ./ExplosiveSOC_LIST_SSF
#include<iostream>
#include <fstream>
#include <cmath>
#include <GL/glew.h>
#include <GL/glut.h>
#include <cuda_runtime_api.h>
#include <cuda_gl_interop.h>
using namespace std;


#define Lx=16
#define Nx=8 //hilos por bloque
const int Mx=(Lx+Nx-1)/Nx; //numero de bloques

//--------------- KERNELS ----------------
__global__ void Suma2Vectores(float * d_a, float * d_b, float * d_c){
  
  int ix=blockIdx.x*blockDim.x+threadIdx.x;
  d_c[ix]=d_a[ix]+d_b[ix];
}

int main(){
  //DECLARAR LAS MATRICES
  float h_a[Lx], h_b[Lx], h_c[Lx];
  float *d_a;  cudaMalloc((void**) &d_a,Lx*sizeof(float));
  float *d_b;  cudaMalloc((void**) &d_b,Lx*sizeof(float));
  float *d_c;  cudaMalloc((void**) &d_c,Lx*sizeof(float));

  //INICIALIZAR LOS DATOS
  //Cargarlos en el Host
  for(ix=0;ix<Lx;ix++){
    h_a[ix]=ix; h_b[ix]=2*ix;}
  
  //Enviarlos al Device
  cudaMemcpy(d_a,h_a,Lx*sizeof(float),cudaMemcpyHostToDevice);  //enviar a
  cudaMemcpy(d_b,h_b,Lx*sizeof(float),cudaMemcpyHostToDevice); // enviar b
  
  //PROCESAR EN LA TARJETA GRAFICA
  dim3 ThreadsPerBlock(Nx,1,1);
  dim3 BlocksPerGrid(Mx,1,1);
  Suma2Vectores<<<BlocksPerGrid,ThreadsPerBlock>>>(d_a,d_b,d_c);

  //IMPRIMIR LOS DATOS
  //Devolverlos al Host
  cudaMemcpy(h_c,d_c,Lx*sizeof(float),cudaMemcpyDeviceToHost);

  //Imprimirlos
  for(ix=0; ix<Lx; ix++)
    cout<<h_c[ix]<<"\t";
  cout<<endl;

  //LIBERAR MEMORIA
  cudaFree(d_a);
  cudaFree(d_b);
  cudaFree(d_c);

  return 0;
}
