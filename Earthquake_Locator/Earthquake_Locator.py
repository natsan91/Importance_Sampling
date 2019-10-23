
# This short project finds the epicenter of an earthquake from observatory data.
# Originally created for Data Science 422 at Northwestern in Feb. 2017 by Nathan Sanford.
# Written for Python 2.7 and updated to Python 3. Requires basemap from mpl_toolkits 
# and pylab for computations.
# Uses least squares iterationto find the **location and origin time** of the earthquake based 
# Creates a map that shows the locations of the observatories as well as iterated guesses for 
# the earthquake.
# See ReadMe for more information.

############################# Import relevant Python utilities #############################
import numpy as np
import matplotlib.pyplot as plt
from math import pi
from map_utilities import makemap, final_map
from solver import ls_solver
# Uncomment the following block for plotting in separate window
#=============================================================================
#try:
#    import IPython
#    shell = IPython.get_ipython()
#    shell.enable_matplotlib(gui='qt') 
#except:
#    pass  
#plt.close() 
#=============================================================================


############################################################################################
# input data from class
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
NN=len(stat_lons)
#############################################################################################
# numerical parameters
numtries=10 # maximum number of iterations allowed
lons=np.zeros((numtries+1,1)) # array of longitude guesses
lats=np.zeros((numtries+1,1)) # array of latitude guesses
tol=2 # tolerance for exiting loop
# the iteration ends when the maximum of the adjustment vector is below the tolerance
# the tolerance can be thought of as seconds or hundreds of km 

# constants
alpha=5.7 # speed km/s of earthquake waves (assumed constant)
lonfac=85.18 # km to degree longitude in the region
latfac=111.19 # km to degree latitude in the region

# center grid on these coordinates and convert to km
lonc=-114.
latc=40.
# define bounds for region we're considering
lond = 8. # half-width
latd = 6. # half-height
stat_x=lonfac*(stat_lons-lonc)
stat_y=latfac*(stat_lats-latc)

# initial guess for epicenter, in coords
ep_lon=lonc + np.random.uniform(-lond,lond)
ep_lat=latc + np.random.uniform(-latd,latd)

# make map with initial guess
lons[0]=ep_lon
lats[0]=ep_lat
makemap(stat_lats,stat_lons,stat_names,lonc,latc,lond,latd,True,0,lons,lats)
plt.show()
plt.pause(0.5)
print("The random initial guess is {:.2f} degrees longitude and {:.2f} degrees latitude.".format(ep_lon,ep_lat))
# initial guess converted to km
xoe=lonfac*(ep_lon-lonc)
yoe=latfac*(ep_lat-latc)

# compute guess for initial time from distance to a research station
# compute distances from initial guess
D=np.matrix(np.ones((NN,1))) # distance
for i in range(0,NN):
    D[i]=np.sqrt((stat_x[i]-xoe)**2.+(stat_y[i]-yoe)**2.+(stat_els[i]+10.)**2.) # distance vector
# using research station nearest to initial guess
temp_stat=np.argmin(D)
# compute initial time guess 
toe=stat_t[temp_stat]-1./alpha*np.sqrt((stat_x[temp_stat]-xoe)**2.+
          (stat_y[temp_stat]-yoe)**2.+(stat_els[temp_stat]+10.)**2.)
answer=[toe,xoe,yoe]
print("The nearest station to these coordinates is ",stat_names[temp_stat],".")
print("The origin time guess based on time to this station is {:.2f}s relative to 14:16:00.".format(toe))

print("Beginning iteration")
# call the solver, can turn off plotting of intermediate solutions by 
# additionally supplying an argument of 0 at the end of the call
answer, GGinv, count = ls_solver(stat_lons,stat_lats,stat_els,stat_t,stat_names,
                                 lonc,latc,lond,latd,lons,lats,answer,tol,numtries)
# put epicenter back into coordinates
answer_lon=answer[1]/lonfac+lonc
answer_lat=answer[2]/latfac+latc


if (count>=numtries):
    # The solver is dependent on the initial guess
    # Most regions converge, but if the initial guess is outside of all of the
    # stations then the solver may not converge.
    print("***did not converge in requested number of iterates*** \n")
else:
    # compute covariance matrix for answer uncertainties
    covmat=0.8**2.*GGinv # assumes the standard deviation in the data is plus/minus 0.8s
    epmat=covmat[1:3,1:3]  # isolate lower right 2 by 2 sub-block for epicenter
    eval, evecs=np.linalg.eig(epmat) # compute for error ellipse
    
    # scale error eigenvectors with eigenvalues
    for i in range(0,len(eval)):
        evecs[:,i]=evecs[:,i]*eval[i]
    # scale to be in longitude and latitude    
    evecs[0,:] = evecs[0,:]/lonfac
    evecs[1,:] = evecs[1,:]/latfac
    
    stddevs = 3 # number of standard deviations for the error ellipse
    
    sem1 = stddevs*2*np.linalg.norm(evecs[:,0]) # semimajor axis diameter length 
    sem2 = stddevs*2*np.linalg.norm(evecs[:,1]) # semiminor axis diameter length 
    rot = np.arctan2(evecs[0,0],evecs[1,0])
    rot = rot/(2*pi)*360 # matplotlib wants rotation angle in degrees
    
    # plot the epicenter in a zoomed-in map with error ellipse
    final_map(answer_lon,answer_lat,sem1,sem2,rot,answer[0],np.sqrt(covmat[0,0])) 
    plt.show()
    
    print("This plot shows the coordinates of the epicenter as well as its origin time.")
    print("Uncertainty in the origin time as well as an error ellipse for the coordinates of the epicenter are also shown. ")
    print("The earthquake's epicenter is at {:.2f} degrees longitude and {:.2f} degrees latitude.".format(answer_lon,answer_lat))
    print("The earthquake occurred at {:.2f}s relative to 14:16:00.".format(answer[0]))