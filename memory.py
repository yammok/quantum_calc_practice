import numpy as np

# === Function
def loadParameter(argv): 
   r=argv[1:4]+argv[7 :10]
   v=argv[4:7]+argv[10:13]

   r=np.array(r,dtype='float64').reshape(2,3)
   v=np.array(v,dtype='float64').reshape(2,3)

   return r,v
