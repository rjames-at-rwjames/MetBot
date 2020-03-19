# To play with harmonics code

from netCDF4 import Dataset as NetCDFFile
import numpy as np
import sys,os
import harmonics_code as hc
cwd=os.getcwd()
sys.path.append(cwd)


### Get directories
bkdir=cwd+"/../../../../CTdata/"
botdir=bkdir+"metbot_multi_dset/"

infile=botdir+"/era/erai/samples/ncfiles/play/lev200.compmean.latsel.nc"

nc=NetCDFFile(infile)
lat = nc.variables['latitude'][:]
lon = nc.variables['longitude'][:]
nlon = len(lon)

gpth=nc.variables['z']
gpth=np.squeeze(gpth)

mn_gpth=np.mean(gpth)
anoms=gpth[:]-mn_gpth

lonnum=np.arange(1,nlon+1)

phi,ps_deg,C = hc.single_harmonic(gpth,lon,h=1)

list_C, list_ps, ex_var_list, amp_var_list = hc.higher_harmonics(gpth,lon)
list_C, list_ps, ex_var_list, amp_var_list = hc.higher_harmonics_fx(gpth,lon,nh=20)

list_C, list_ps, ex_var_list, amp_var_list = hc.higher_harmonics_fx(anoms,lonnum,nh=20)
list_C, list_ps, ex_var_list, amp_var_list = hc.higher_harmonics(anoms,lonnum)

