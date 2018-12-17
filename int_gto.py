import math
import numpy   as np

# --- my module
import int_gto as gto
import math_f  as mf

# === Parameter
NG=3
PI=math.pi
C_STO3G=(0.444635,0.535328,0.154329)
E_STO3G=(0.109818,0.405771,2.227660)

# ==== Function
def calcVG(r,ia,i,j): 
   c=np.array(C_STO3G,dtype='float64')
   e=np.array(E_STO3G,dtype='float64')
   d_ab=np.linalg.norm(r[i]-r[j])

   v_g=0.
   for ig in range(NG): 
      for jg in range(NG): 
         p=e[ig]+e[jg]
         rp=(e[ig]*r[i]+e[jg]*r[j])/p
         d_pc=np.linalg.norm(r[ia]-rp)
         v_g+=c[ig]*c[jg]*gto.VInt(e[ig],e[jg],p,d_ab,d_pc)

   return v_g

#def calcVG(r,i,j): 
#   c=np.array(C_STO3G,dtype='float64')
#   e=np.array(E_STO3G,dtype='float64')
#   d_ab=np.linalg.norm(r[i]-r[j])
#
#   v_g=0.
#   for na in range(2): 
#      for ig in range(NG): 
#         for jg in range(NG): 
#            p=e[ig]+e[jg]
#            rp=(e[ig]*r[i]+e[jg]*r[j])/p
#            d_pc=np.linalg.norm(r[na]-rp)
#            v_g+=c[ig]*c[jg]*gto.VInt(e[ig],e[jg],p,d_ab,d_pc)
#
#   return v_g

def VInt(e_ig,e_jg,p,d_ab,d_pc): 
   return -2.*PI/p*math.exp(-e_ig*e_jg/p*d_ab*d_ab)*mf.F0(p*d_pc*d_pc) \
          *(2.*e_ig/PI)**0.75*(2.*e_jg/PI)**0.75

def calcPG(r,i,j,ix): 
   c=np.array(C_STO3G,dtype='float64')
   e=np.array(E_STO3G,dtype='float64')
   d_ab=np.linalg.norm(r[i]-r[j])

   p_g=0.
   for ig in range(NG): 
      for jg in range(NG): 
         p=e[ig]+e[jg]
         rp=(e[ig]*r[i]+e[jg]*r[j])/p
         d_pb_l=rp[ix]-r[j][ix]
         p_g+=c[ig]*c[jg]*gto.PInt(d_pb_l,e[ig],e[jg],p,d_ab)

   return -1j*p_g

def PInt(d_pb_l,e_ig,e_jg,p,d_ab): 
   return -2.*e_jg*math.exp(-e_ig*e_jg/p*d_ab*d_ab)*d_pb_l*(PI/p)**1.5 \
          *(2.*e_ig/PI)**0.75*(2.*e_jg/PI)**0.75

def calcXG(r,i,j,ia,ix): 
   c=np.array(C_STO3G,dtype='float64')
   e=np.array(E_STO3G,dtype='float64')
   d_ab=np.linalg.norm(r[i]-r[j])

   x_g=0.
   for ig in range(NG): 
      for jg in range(NG): 
         p=e[ig]+e[jg]
         rp=(e[ig]*r[i]+e[jg]*r[j])/p
         d_pc_l=rp[ix]-r[ia][ix]
         x_g+=c[ig]*c[jg]*gto.XInt(d_pc_l,e[ig],e[jg],p,d_ab)

   return  x_g

def XInt(d_pc_l,e_ig,e_jg,p,d_ab): 
   return math.exp(-e_ig*e_jg/p*d_ab*d_ab)*d_pc_l*(PI/p)**1.5 \
          *(2.*e_ig/PI)**0.75*(2.*e_jg/PI)**0.75

def calcEG(r,i,j,ia,ix): 
   c=np.array(C_STO3G,dtype='float64')
   e=np.array(E_STO3G,dtype='float64')
   d_ab=np.linalg.norm(r[i]-r[j])

   e_g=0.
   for ig in range(NG): 
      for jg in range(NG): 
         p=e[ig]+e[jg]
         rp=(e[ig]*r[i]+e[jg]*r[j])/p
         d_cp_l=r[ia][ix]-rp[ix]
         d_pc=np.linalg.norm(r[ia]-rp)
         e_g+=c[ig]*c[jg]*gto.EInt(d_cp_l,e[ig],e[jg],p,d_ab,d_pc)

   if i==j and i==ia: 
      e_g=e_g*2. 

   return e_g

def EInt(d_cp_l,e_ig,e_jg,p,d_ab,d_pc): 
   return -4.*PI*math.exp(-e_ig*e_jg/p*d_ab*d_ab)*d_cp_l*mf.F2(p*d_pc*d_pc) \
          *(2.*e_ig/PI)**0.75*(2.*e_jg/PI)**0.75
