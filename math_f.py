import math 
import numpy as np

from scipy import integrate

# === Parameter
PI=math.pi
X0=0
X1=1
NX=100

# === Function
def F0(t): 
   if t<1.e-06: 
      return 1.-t/3.
   else:    
      return 0.5*math.sqrt(PI/t)*math.erf(math.sqrt(t))

def F2(t):
   x=np.linspace(X0,X1,NX)
   y=x*x*np.exp(-t*x*x)
   return integrate.simps(y,x)
