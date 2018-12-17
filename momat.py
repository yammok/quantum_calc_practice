import math
import numpy as np

import aomat as ao

# === Parameter
NA=2
NX=3

# === Function
def calcC(s):
   # calculate AO coefficients
   s_ab=s[0][1]

   c=np.zeros((2,2),dtype='float64')
   c[0][0]= math.sqrt(2.*(1.+s_ab))
   c[0][1]= math.sqrt(2.*(1.+s_ab))
   c[1][0]= math.sqrt(2.*(1.-s_ab))
   c[1][1]=-math.sqrt(2.*(1.-s_ab))

   return c

def diagH0(d,rn,s,t,v): 
   h=np.zeros((2,2),dtype='float64')

   h0=t
   for ia in range(NA):
     h0+=v[ia]
   h[0][0]=(h0[0][0]+h0[0][1])/(1.+s[0][1])+1./d
   h[1][1]=(h0[0][0]-h0[0][1])/(1.-s[0][1])+1./d 

   return h,calcC(s)

def calcP(d_ij,rn,c): 
   p_ao=ao.calcP(d_ij,rn)

   p_mo=np.zeros((NX,2,2),dtype='complex128')
   for ix in range(NX): 
      for imo in range(2):
         for jmo in range(2):
            for iao in range(2):
               for jao in range(2):
                  p_mo[ix][imo][jmo]+=c[imo][iao]*p_ao[ix][iao][jao] \
                                     *c[jmo][jao]

   return p_mo

def calcX(d_ij,rn,c): 
   x_ao=ao.calcX(d_ij,rn)

   x_mo=np.zeros((2,NX,2,2),dtype='float64')
   for ia in range(NA): 
      for ix in range(NX): 
         for imo in range(2):
            for jmo in range(2):
               for iao in range(2):
                  for jao in range(2):
                     x_mo[ia][ix][imo][jmo]+=c[imo][iao]*c[jmo][jao] \
                                            *x_ao[ia][ix][iao][jao] 

   return x_mo

def calcE(d_ij,rn,c): 
   e_ao=ao.calcE(d_ij,rn)

   e_mo=np.zeros((2,NX,2,2),dtype='float64')
   for ia in range(NA): 
      for ix in range(NX): 
         for imo in range(2):
            for jmo in range(2):
               for iao in range(2):
                  for jao in range(2):
                     e_mo[ia][ix][imo][jmo]+=c[imo][iao]*c[jmo][jao] \
                                            *e_ao[ia][ix][iao][jao] 
   return e_mo

def calcV(v_ao,c): 
   v_mo=np.zeros((NA,2,2),dtype='float64')
   for ia in range(NA): 
      for imo in range(2): 
         for jmo in range(2): 
            for iao in range(2): 
               for jao in range(2):  
                  v_mo[ia][imo][jmo]+=c[imo][iao]*v_ao[ia][iao][jao] \
                                     *c[jmo][jao]

   return v_mo

def calcA(d_ij,rn,vn,v,x,e,c): 
   v_ao=v
   v_mo=calcV(v_ao,c)

   a_mo=np.zeros((NX,2,2),dtype='float64')
   a_mo+=firstTermA(vn,v_mo)
#   a_mo+=secndTermA(vn,v_mo)

   print(a_mo)

   return 0.5*a_mo

def firstTermA(vn,v_mo): 
   a_mo=np.zeros((NX,2,2),dtype='float64')
   for ix in range(NX): 
      for imo in range(2): 
         for jmo in range(2): 
            a_mo[ix][imo][jmo]+=dotRV(ix,imo,jmo,vn,v_mo)

   return a_mo

def dotRV(ix,imo,jmo,vn,v_mo): 
   term=0.
   for ia in range(NA):
      term+=vn[ia][ix]*v_mo[ia][imo][jmo]

   return term

def secndTermA(ix,imo,jmo,vn,x,e):
   tmp=np.zeros((NA,2,2),dtype='float64')
   term=0.

   return term
