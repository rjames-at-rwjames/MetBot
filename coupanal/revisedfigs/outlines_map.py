# To plot map with domains used on top

import os
import sys
import math as mh

runoffline=True
if runoffline==True:
    import matplotlib
    matplotlib.use('Agg')

import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.basemap import cm
from matplotlib.patches import Polygon

cwd=os.getcwd()
sys.path.append(cwd)
sys.path.append(cwd+'/..')
sys.path.append(cwd+'/../..')
import MetBot.mytools as my
import MetBot.MetBlobs as blb

# Function to draw rectangle
def draw_screen_poly( lats, lons, m, edgecol='k', fill=False, facecol=None, ls='-', lw=1, transparency=None, hatch=False, zord=1):
    # facecol can be 'None'
    # transparency can be 'None'
    x, y = m( lons, lats )
    xy = zip(x,y)
    poly = Polygon( xy, edgecolor=edgecol, facecolor=facecol, fill=fill, ls=ls,\
                    lw=lw, alpha=transparency, hatch=hatch, zorder=zord)
    plt.gca().add_patch(poly)

# Function to draw arrow at angle
def draw_arrow(x,y,angle,length,cl='k',edgecol='k', zord=1):
    # Code adapted from: https://stackoverflow.com/questions/15971768/drawing-arrow-in-x-y-coordinate-in-python
    # x is the base of the arrow (lon)
    # y is the base of the arrow (lat)
    # angle is the angle in degrees
    # length is the desired length of the arrow

    # convert from our coordinates (+x,+y,-x,-y) as (-90,180,0,+90)
    # to Cartesian where (+x, +y, -x, -y) as (0, 90, 180, 270)
    cartesianAngleRadians = (180+angle)*mh.pi/180.0
    terminus_x = x + length * mh.cos(cartesianAngleRadians)
    terminus_y = y + length * mh.sin(cartesianAngleRadians)
    arlen_x=terminus_x-x
    arlen_y=terminus_y-y

    arlen_x2= -arlen_x
    arlen_y2= -arlen_y

    plt.arrow(x,y,arlen_x,arlen_y,width=.07, fc=cl, ec=edgecol, zorder=zord)
    plt.arrow(x,y,arlen_x2,arlen_y2,width=.07, fc=cl, ec=edgecol,zorder=zord)


figdim=[8,6]

lat=np.arange(-70.0,20.0,10.0)
lon=np.arange(-2.5,120.0,10.0)
latsp = 20.  # lat spacing
lonsp = 25.  # lon spacing

### Get directories
bkdir=cwd+"/../../../../CTdata/"
botdir=bkdir+"metbot_multi_dset/"
figdir=botdir+"/histpaper_v2_figs/outlines_map/"
my.mkdir_p(figdir)

# Set up plot
print "Setting up plot..."
g, ax = plt.subplots(figsize=figdim)

m = blb.AfrBasemap2(lat, lon, latsp, lonsp, drawstuff=True, prj='cyl', rsltn='l', \
                    fontdict={'fontsize': 12, 'fontweight': 'normal'})

# Redraw map
m.drawcountries()
m.drawcoastlines()

# Plot MetBot full domain
fulldom_lats = [0, 0, -60, -60]
fulldom_lons = [7.5, 100.0, 100.0, 7.5]
draw_screen_poly( fulldom_lats, fulldom_lons, m , ls='--', lw=3, edgecol='dodgerblue',zord=10)

# Plot tropical box
tropdom_lats = [0,0,-10,-10]
tropdom_lons = [0,100.0,100.0,0]
draw_screen_poly( tropdom_lats, tropdom_lons, m, edgecol='k', fill=False,\
                  facecol=None, ls='-', transparency=0.2,zord=2, hatch='/')

# Plot Congo region
congo_lats = [0,0,-10,-10]
congo_lons = [10,30,30,10]
draw_screen_poly( congo_lats, congo_lons, m, edgecol='greenyellow', fill=True, \
                  facecol='greenyellow', ls='-', lw=1, transparency=0.7, zord=9)

# Plot subtropical southern Africa
str_lats = [-20,-20,-35,-35]
str_lons = [15,40,40,15]
draw_screen_poly( str_lats, str_lons, m, edgecol='deepskyblue', fill=True,\
                  facecol='deepskyblue', ls='-', lw=1, transparency=0.3,zord=3)

# Plot subtropical full dom
strfl_lats = str_lats
strfl_lons = [7.5,100,100,7.5]
draw_screen_poly( strfl_lats, strfl_lons, m, edgecol='lightsteelblue', fill=True,\
                  facecol='lightsteelblue', ls='-', lw=1, transparency=0.2,zord=2)

# Plot Froude dom box
fr_lats = [-3,-3,-11,-11]
fr_lons = [36,42,42,36]
draw_screen_poly( fr_lats, fr_lons, m, ls='-', lw=2, edgecol='green',zord=10)


# Plot lines

# dividing line between continental and madagascan
latlims=np.asarray([fulldom_lats[0],fulldom_lats[2]])
lonmada=np.asarray([55,55])
m.plot(lonmada,latlims,ls='-',lw=2,color='orange',zorder=5)

# Plot line & dot for centroid
cen_lat_range=np.asarray([-22,-32])
cen_lon=33.0
cen_lon_line=np.asarray([cen_lon,cen_lon])
m.plot(cen_lon_line,cen_lat_range,ls='-', lw=2, color='violet',zorder=8)

# Plot star for centroid
best_lat=np.mean(cen_lat_range)
m.plot(cen_lon,best_lat,marker='*',markersize=15,markeredgecolor='k',markerfacecolor='fuchsia',markeredgewidth=1.5, zorder=10)

# Plot arrows (following NH method in MetBlobs)
cl='k'
angles=np.asarray([-60.0,-25.0])
draw_arrow(cen_lon,best_lat,angles[0],5,cl='fuchsia',edgecol='fuchsia',zord=9)
draw_arrow(cen_lon,best_lat,angles[1],5,cl='fuchsia',edgecol='fuchsia',zord=9)

# Now Madagascan
# Plot line & dot for centroid
cen_lat_range=np.asarray([-22,-32])
cen_lon=58.0
cen_lon_line=np.asarray([cen_lon,cen_lon])
m.plot(cen_lon_line,cen_lat_range,ls='-', lw=2, color='royalblue',zorder=8)

# Plot star for centroid
best_lat=np.mean(cen_lat_range)
m.plot(cen_lon,best_lat,marker='*',markersize=15,markeredgecolor='k',markerfacecolor='blue', markeredgewidth=1.5, zorder=10)

# Plot arrows (following NH method in MetBlobs)
cl='k'
angles=np.asarray([-50.0,-15.0])
draw_arrow(cen_lon,best_lat,angles[0],5,cl='blue',edgecol='blue',zord=9)
draw_arrow(cen_lon,best_lat,angles[1],5,cl='blue',edgecol='blue',zord=9)

print "Finalising plot..."
figname = figdir + 'outlines_map.png'
print 'saving figure as '+figname
plt.savefig(figname, dpi=600)
plt.close()