import sys
import numpy   as np

from scipy.spatial import distance

# --- my module
import memory  as mem
import print   as prt
import momat   as mo
import aomat   as ao

# === Function
def calcHamiltonian(rn,vn): 
   d_ij=distance.cdist(rn,rn,metric='euclidean')
   d   =np.linalg.norm(rn[0]-rn[1])

   h=np.zeros((2,2),dtype='complex128')
   s=ao.calcS(d_ij)
   t=ao.calcT(d_ij)
   v=ao.calcV(d_ij,rn)

   h0,c=mo.diagH0(d,rn,s,t,v)
   p   =mo.calcP(d_ij,rn,c)
   x   =mo.calcX(d_ij,rn,c)
   e   =mo.calcE(d_ij,rn,c)
   # this term leads to non Hermite matrix, so we need to fix
   a   =mo.calcA(d_ij,rn,vn,v,x,e,c) 
#   print(a)

   # make a function to calculate Hamitonian matrix

   return h,d

def PESCalculator(argv):
   r_nu,v_nu=mem.loadParameter(argv)
   h,d=calcHamiltonian(r_nu,v_nu)
   prt.printResults(d,h)

# === Main
argv=sys.argv
PESCalculator(argv)
