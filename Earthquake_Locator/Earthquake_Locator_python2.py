
# This short project finds the epicenter of an earthquake from observatory data.
# Originally created for Data Science 422 at Northwestern in Feb. 2017 by Nathan Sanford.
# Written for Python 2.7 clui. Requires basemap from mpl_toolkits and pylab for computations.
# Uses regularized least squares to find the **location and origin time** of the earthquake based 
# Creates a map that shows the locations of the observatories as well as iterated guesses for 
# the 

############################# Import relevant Python utilities #############################
import pylab
from mpl_toolkits.basemap import Basemap
import numpy as np
import matplotlib.pyplot as plt
import warnings
import matplotlib.cbook
warnings.filterwarnings("ignore",category=matplotlib.cbook.mplDeprecation)

# input data
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
#############################################################################################
# initial guess for epicenter
ep_lon=-120.
ep_lat=43.
# constants
alpha=5.7
lonfac=85.18
latfac=111.19
NN=len(stat_lats)
# center grid on these coordinates and convert to km
lonc=-114.
latc=40.
stat_x=lonfac*(stat_lons-lonc)
stat_y=latfac*(stat_lats-latc)
# initial guess converted to km
xoe=lonfac*(ep_lon-lonc)
yoe=latfac*(ep_lat-latc)
# compute guess for initial time from distance to a research station
temp_stat=0 # using first station
# compute initial time guess 
toe=stat_t[temp_stat]-1./alpha*np.sqrt((stat_x[temp_stat]-xoe)**2.+(stat_y[temp_stat]-yoe)**2.+(stat_els[temp_stat]+10.)**2.)
answer_old=[toe,xoe,yoe]
print "Initial guess=",toe,ep_lon,ep_lat
# make G matrix and data
G=np.matrix(np.ones((NN,3))) # design matrix
d=np.matrix(np.ones((NN,1))) # data
D=np.matrix(np.ones((NN,1))) # distance

for i in range(0,NN):
	D[i]=np.sqrt((stat_x[i]-xoe)**2.+(stat_y[i]-yoe)**2.+(stat_els[i]+10.)**2.) # distance vector
	d[i]=stat_t[i]-toe-1./alpha*D[i] # data vector
	# define temporary constant which is same for each row
	a=1./alpha*1./D[i]
	G[i,0]=1.
	G[i,1]=a*(xoe-stat_x[i])
	G[i,2]=a*(yoe-stat_y[i])
# Compute least squares solution
GG=np.dot(np.transpose(G),G)
GGinv=np.linalg.inv(GG)
GGG=np.dot(GGinv,np.transpose(G))
mod=np.dot(GGG,d) # model data 
# correct initial guess for epicenter with model data
answer=np.asarray(mod).reshape(-1)+answer_old # converted to array from column vector
# convert answer into longitude and latitude for plotting
answer_lon=answer[1]/lonfac+lonc
answer_lat=answer[2]/latfac+latc
print "Solution after ",1," iterarion is ",answer[0],answer_lon,answer_lat

# make while loop
tol=0.1 # tolerance for exit loop
check=10.
count=0
numtries=50 # maximum number of iterations allowed
while check>tol and count<numtries:
	if (count==0): # create map from first iteration
		# make map with initial guess from previous iteration
		m = Basemap(llcrnrlon=-122.,llcrnrlat=32.,urcrnrlon=-108.,urcrnrlat=45.)
		# add stations
		x, y = m(stat_lons,stat_lats)
		m.scatter(x, y, marker='D',color='m')
		# add initial guess
		x1, y1 = m(ep_lon,ep_lat)
		m.scatter(x1,y1,marker='D',color='g')
		plt.annotate("Initial guess for epicenter", xy=(x1, y1),  xycoords='data',
			xytext=(-50, 5), textcoords='offset points',
			color='g')
		# add labels for stations
		x2, y2 = (-20,5)
		for i in range(0,NN):
			plt.annotate(stat_names[i], xy=(x[i], y[i]),  xycoords='data',
						xytext=(x2, y2), textcoords='offset points',
						color='m')
		# map cosmetics
		m.drawstates()
		m.drawcountries()
		m.drawmapboundary(fill_color='aqua')
		m.fillcontinents(color='coral',lake_color='aqua',zorder=0)
		m.drawcoastlines()
		# add least squares solution
		x3, y3 = m(answer_lon,answer_lat)
		m.scatter(x3,y3,marker='D',color='b')		
		mytitle='After '+str(count+1)+' iterations'
		plt.title(mytitle)
		plt.ion()
		plt.show()
		plt.pause(0.2)

	# save 1st solution for comparison later
	mod_old=mod
	toe=answer[0]
	xoe=answer[1]
	yoe=answer[2]
	answer_old=answer
	# make G matrix and data from new initial guess
	for i in range(0,NN):
		D[i]=np.sqrt((stat_x[i]-xoe)**2.+(stat_y[i]-yoe)**2.+(stat_els[i]+10.)**2.) # distance vector
		d[i]=stat_t[i]-toe-1./alpha*D[i] # data vector
		# define temporary constant which is same for each row
		a=1./alpha*1./D[i]
		G[i,0]=1.
		G[i,1]=a*(xoe-stat_x[i])
		G[i,2]=a*(yoe-stat_y[i])
	# Compute least squares solution
	GG=np.dot(np.transpose(G),G)
	GGG=np.dot(np.linalg.inv(GG),np.transpose(G))
	mod=np.dot(GGG,d)
	# controls for exiting loop
	check=np.max(mod_old-mod)
	# correct epicenter with model data
	answer=np.asarray(mod).reshape(-1)+answer_old # convert to array from column vector
	
	# overwrite old solution with red on map
	m.scatter(x3,y3,marker='D',color='r')
	# get epicenter in latitude and longitude for map
	answer_lon=answer[1]/lonfac+lonc
	answer_lat=answer[2]/latfac+latc
	print "Solution after ",count+2," iterarion is ",answer[0],answer_lon,answer_lat
	# add least squares solution to map
	x3, y3 = m(answer_lon,answer_lat)
	m.scatter(x3,y3,marker='D',color='b')
	if (check<tol): # this plot contains the final solution
		mytitle='Final answer after '+str(count+2)+' iterations'
		mylabel='t='+str(answer[0])
		plt.annotate(mylabel, xy=(x3, y3),  xycoords='data',
					xytext=(x2, y2), textcoords='offset points',
					color='b')
	else: # this plot contains an intermediate solution
		mytitle='After '+str(count+2)+' iterations'
	plt.title(mytitle)
	plt.show()
	if (check<tol):
		print "Solution converged, waiting to continue"
		raw_input()
		plt.close()
	else:
		plt.pause(1)
	count=count+1
print "The final answer is ",answer[0],answer_lon,answer_lat
if (count>=numtries):
	print "***did not converge*** \n"
# compute covariance matrix
covmat=0.8**2.*GGinv
print "The covariance matrix is"
print covmat
epmat=np.matrix(np.ones((2,2)))
# isolate lower right 2 by 2 sub-block
epmat[0,0]=covmat[1,1]
epmat[0,1]=covmat[1,2]
epmat[1,0]=covmat[2,1]
epmat[1,1]=covmat[2,2]
print "The elements of the covariance matrix associated with the epicenter are "
print epmat
print "time std=",np.sqrt(covmat[0,0])
eval, evecs=np.linalg.eig(epmat)
print "The eigenvalues associated withe errors in the epicenter are "
print eval
print "The concomitant eigenvectors are "
print evecs
# make another map with initial guess and final solution (no intermediate solutions)
pylab.ion()
map = Basemap(llcrnrlon=-122.,llcrnrlat=32.,urcrnrlon=-108.,urcrnrlat=45.)
# add stations
x, y = map(stat_lons,stat_lats)
map.scatter(x, y, marker='D',color='m')
# add initial guess
x1, y1 = map(ep_lon,ep_lat)
map.scatter(x1,y1,marker='D',color='g')
plt.annotate("Initial guess for epicenter", xy=(x1, y1),  xycoords='data',
			xytext=(-50, 5), textcoords='offset points',
			color='g')
# add least squares solution
x3, y3 = map(answer_lon,answer_lat)
map.scatter(x3,y3,marker='D',color='b',zorder=2)
mylabel="t="+str(answer[0])+r"$\pm$"+str(np.sqrt(covmat[0,0]))
x4,y4=(-80,5)
plt.annotate(mylabel, xy=(x3, y3),  xycoords='data',
			xytext=(x4, y4), textcoords='offset points',
			color='b')
# map cosmetics
map.drawstates()
map.drawcountries()
map.drawmapboundary(fill_color='aqua')
map.fillcontinents(color='coral',lake_color='aqua',zorder=0)
map.drawcoastlines()
# add labels for stations
x2, y2 = (-20,5)
for i in range(0,NN):
	plt.annotate(stat_names[i], xy=(x[i], y[i]),  xycoords='data',
				xytext=(x2, y2), textcoords='offset points',
				color='m')
mytitle='Final answer after '+str(count+1)+' iterations'
pylab.title(mytitle)
raw_input()