def makemap(stat_lats,stat_lons,stat_names,itq=False,iterate=0,i_lons=0,i_lats=0,time=0):
    import os
    os.environ["PROJ_LIB"] = "C:\\Users\\natha\\Anaconda3\\Library\\share"; #fixr
    from mpl_toolkits.basemap import Basemap#, shiftgrid, cm
#    import numpy as np
    import matplotlib.pyplot as plt
    import warnings
    import matplotlib.cbook
#    from ellipsefunction import ellipse
#    from matplotlib.patches import Ellipse
    warnings.filterwarnings("ignore",category=matplotlib.cbook.mplDeprecation)
#    Basemap.ellipse=Ellipse
    
#    stat_lats=np.array([40.20,43.56,38.35,43.92,36.95,42.77,37.00,44.40,
#    		37.63,37.05,36.47,40.83])
#    stat_lons=np.array([-112.81,-114.41,-115.59,-110.58,-116.25,-109.56,-117.37,-110.58,
#    		-118.84,-118.23,-116.86,-108.26])
#    stat_names=["DUG","HLID","R11A","I17A","TPNV","BW06","GRA","H17A",
#    		"MLAC","TIN","FUR","N20A"]
    m = Basemap(llcrnrlon=-122.,llcrnrlat=32.,urcrnrlon=-108.,urcrnrlat=45.,resolution='l')
    # add stations
    x, y = m(stat_lons, stat_lats)
    m.scatter(x, y, marker='D',color='m')
    # map cosmetics
    m.drawstates()
    m.drawcountries()
    #m.drawmapboundary(fill_color='aqua')
    m.drawcoastlines()
    #m.fillcontinents(color='coral',lake_color='aqua',zorder=0)
    m.drawlsmask(land_color='coral',ocean_color='aqua',grid=1.25)
    # add labels for stations
    x2, y2 = (-20,5)
    for i in range(0,len(stat_lats)):
    	plt.annotate(stat_names[i], xy=(x[i], y[i]),  xycoords='data',
                    xytext=(x2, y2), textcoords='offset points',
                    color='m')
    if itq:
        if iterate==1:
            x, y = m(i_lons[0],i_lats[0])
            m.scatter(x,y,marker='D',color='g')	
        elif iterate==2:
            x, y = m(i_lons[0],i_lats[0])
            m.scatter(x,y,marker='D',color='g')
            xx, yy = m(i_lons[1],i_lats[1])
            m.scatter(xx,yy,marker='D',color='r')
        elif iterate>2:
            x, y = m(i_lons[0],i_lats[0])
            m.scatter(x,y,marker='D',color='g')
            xx, yy = m(i_lons[1:iterate],i_lats[1:iterate])
            m.scatter(xx,yy,marker='D',color='r')
        xcurr, ycurr = m(i_lons[iterate],i_lats[iterate])
        m.scatter(xcurr,ycurr,marker='D',color='b')
        mytitle='After '+str(iterate)+' iterations'
        plt.title(mytitle)
        
#    x, y = m(-115.59,38.35)
#    x2,y2=m(-112.81,40.20)
##    m.scatter(x,y,marker='D',color='c')
##    m.ellipse(x,y,3,1.5,25, facecolor='green', zorder=10,alpha=0.5)  
#    el=Ellipse((x,y),2*np.abs(x2-x),2*np.abs(y2-y), facecolor='green', zorder=10,alpha=0.5) 
#    plt.gca().add_patch(el)   
#    m.scatter(x2,y,marker='D',color='g')
        
    plt.show()
    plt.pause(0.2)
    
    return m

def add_ellipse(lonc,latc,sem1,sem2,rot=0):
    from matplotlib.patches import Ellipse
    import matplotlib.pyplot as plt
    el=Ellipse((lonc,latc),sem1,sem2,rot,facecolor='green', zorder=10,alpha=0.5) 
    plt.gca().add_patch(el) 

import numpy as np
from math import pi
#from matplotlib.patches import Ellipse
#import matplotlib.pyplot as plt
#from ellipsefunction import Basemap

stat_lons=np.array([-112.81,-114.41,-115.59,-110.58,-116.25,-109.56,-117.37,-110.58,
        -118.84,-118.23,-116.86,-108.26])
stat_lats=np.array([40.20,43.56,38.35,43.92,36.95,42.77,37.00,44.40,
        37.63,37.05,36.47,40.83])
stat_els=np.array([1.477,1.772,1.756,2.134,1.600,2.224,0.689,2.400,
        2.162,1.197,-0.37,2.192])
stat_t=np.array([33.262,44.200,48.640,69.650,69.314,70.290,73.800,73.526,
        74.990,76.908,77.749,79.414])
stat_names=["DUG","HLID","R11A","I17A","TPNV","BW06","GRA","H17A",
        "MLAC","TIN","FUR","N20A"]

mymap = makemap(stat_lats,stat_lons,stat_names)
x, y = mymap(-115.59,38.35)
x2,y2=mymap(-112.81,40.20)
ew=2*np.abs(x2-x)
eh=2*np.abs(y2-y)
semm=np.sqrt(ew**2+eh**2)
my_ang=np.arctan(eh/ew)
my_deg=my_ang/(2*pi)*360
#    m.scatter(x,y,marker='D',color='c')
#    m.ellipse(x,y,3,1.5,25, facecolor='green', zorder=10,alpha=0.5)  
add_ellipse(x,y,semm,1,my_deg) 
mymap.scatter(x2,y,marker='D',color='g')
mymap.scatter(x,y2,marker='D',color='c')
#plt.gca().add_patch(el)  
#x, y = mymap(-115,42)
#mymap.ellipse(x,y,10,5)
