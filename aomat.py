import math
import numpy   as np

# --- my module
import int_gto as gto

# === Parameter
NX=3
NA=2

# === Function
def calcS(d): 
   return (1+d+d*d/3.)*np.exp(-d)

def calcT(d): 
   return 0.5*np.exp(-d)*(1.+d-d*d/3.)

def calcV(d,r):
   v=np.zeros((NA,2,2),dtype='float64')
   for ia in range(NA): 
      for iao in range(2): 
         for jao in range(2): 
            v[ia][iao][jao]=gto.calcVG(r,ia,iao,jao) 

   return v

def calcP(d,r):
   p=np.zeros((NX,2,2),dtype='complex128')
   for ix in range(3): 
      for iao in range(2): 
         for jao in range(2): 
            p[ix][iao][jao]=gto.calcPG(r,iao,jao,ix)

   return p

def calcX(d,r): 
   x=np.zeros((2,NX,2,2),dtype='float64')
   for ia in range(2): 
      for ix in range(NX): 
         for iao in range(2): 
            for jao in range(2): 
               x[ia][ix][iao][jao]=gto.calcXG(r,iao,jao,ia,ix)

   return x

def calcE(d,r): 
   e=np.zeros((2,NX,2,2),dtype='float64')
   for ia in range(2): 
      for ix in range(NX): 
         for iao in range(2): 
            for jao in range(2): 
               e[ia][ix][iao][jao]=gto.calcEG(r,iao,jao,ia,ix)

   return e

def calcA(d,rn,vn,v,x,e): 
   a=np.zeros((NX,2,2),dtype='float64')
   for  ix in range(NX): 
      for iao in range(2): 
         for jao in range(2): 
            for  ia in range(2): 
               a[ix][iao][jao]+=vn[ia][ix]*v[ia][iao][jao]

   return a*0.5
