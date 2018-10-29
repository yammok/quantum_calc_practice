import sys
import numpy as np
import math

from scipy.spatial import distance

# === Parameter
C_STO3G=(0.444635,0.535328,0.154329)
E_STO3G=(0.109818,0.405771,2.227660)
PI=math.pi

# === Function
def loadParameter(argv): 
   r=argv[1:4]+argv[7 :10]
   v=argv[4:7]+argv[10:13]

   r=np.array(r,dtype='float64').reshape(2,3)
   v=np.array(v,dtype='float64').reshape(2,3)

   return r,v

def calcS(d): 
   return (1+d+d*d/3.)*np.exp(-d)

def calcT(d): 
   return 0.5*np.exp(-d)*(1.+d-d*d/3.)

def F0(t): 
   if t<1.e-06: 
      return 1.-t/3.
   else:    
      return 0.5*math.sqrt(PI/t)*math.erf(math.sqrt(t))

def VInt(e_ig,e_jg,p,d_ab,d_pc): 
   return -2.*PI/p*math.exp(-e_ig*e_jg/p*d_ab*d_ab)*F0(p*d_pc*d_pc) \
          *(2.*e_ig/PI)**0.75*(2.*e_jg/PI)**0.75

def calcVG(r,i,j): 
   c=np.array(C_STO3G,dtype='float64')
   e=np.array(E_STO3G,dtype='float64')
   d_ab=np.linalg.norm(r[i]-r[j])

   v=0.
   for na in range(2): 
      for ig in range(3): 
         for jg in range(3): 
            p =e[ig]+e[jg]
            rp=(e[ig]*r[i]+e[jg]*r[j])/p
            d_pc=np.linalg.norm(r[na]-rp)
            v+=c[ig]*c[jg]*VInt(e[ig],e[jg],p,d_ab,d_pc)

   return v

def calcV(d,r):
   v=np.zeros((2,2),dtype='float64')
   for i in range(2): 
      for j in range(2): 
         v[i][j]=calcVG(r,i,j) 

   return v

def calcHamiltonian(r,v): 
   d_ij=distance.cdist(r,r,metric='euclidean')
   d   =np.linalg.norm(r[0]-r[1])

   s=calcS(d_ij)
   t=calcT(d_ij)
   v=calcV(d_ij,r)
   h_bare=t+v

   e=np.zeros((2),dtype='float64')
   e[0]=(h_bare[0][0]+h_bare[0][1])/(1.+s[0][1])+1./d
   e[1]=(h_bare[0][0]-h_bare[0][1])/(1.-s[0][1])+1./d

   h=e # Originally, the dimension of this matrix is 2*2

   return h,d

def printResults(d,e): 
   print(d,e[0],e[1])

def PESCalculator(argv):
   r_nu,v_nu=loadParameter(argv)
   h,d=calcHamiltonian(r_nu,v_nu)
   printResults(d,h)

# === Main
argv=sys.argv
PESCalculator(argv)
