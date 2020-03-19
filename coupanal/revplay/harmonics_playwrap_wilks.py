# To play with harmonics code

import numpy as np
import sys,os
import harmonics_code as hc
cwd=os.getcwd()
sys.path.append(cwd)


## Wilks data
y=[22.2,22.7,32.2,44.4,54.8,64.3,68.8,67.1,60.2,49.5,39.3,27.4]
y=np.asarray(y)
t=np.arange(1,13)


phi,ps_deg,C = hc.single_harmonic(y,t,h=1)
list_C, list_ps, ex_var_list, amp_var_list = hc.higher_harmonics(y,t)

