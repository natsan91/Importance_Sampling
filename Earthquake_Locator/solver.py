def ls_solver(stat_lons,stat_lats,stat_els,stat_t,stat_names,
           lonc,latc,lon_delt,lat_delt,lons,lats,answer,tol,numtries,plotq=1):
    import numpy as np
    from map_utilities import makemap
    import matplotlib.pyplot as plt 
#-------------------------------------------------------------------------
# The parameters of the solver function are : 
# stat_lons, stat_lats, stat_els, stat_t, stat_names are info about the stations
# lonc,latc,lon_delt,lat_delt are info about the region we're considering 
# (central coords and half-width and half-height of region)
# lons, lats are lists of the iterated guesses in degrees (to be filled in)
# answer is the initial guess (origin time and epicenter in km from center)
# tol, numtries are convergence tolerance and max number of iterations allowed
# have the option to turn off plotting intermediate steps (on by default)

# solver checks L1 norm of the LS correction
# solver returns :
# answer contains the origin time and epicenter in km from center
# count=# of iterations performed
# GGinv is the covariance matrix, for uncertainties in the answer  
#-------------------------------------------------------------------------
    # constants
    alpha=5.7 # speed km/s of earthquake waves (assumed constant)
    lonfac=85.18 # km to degree longitude in the region
    latfac=111.19 # km to degree latitude in the region
    # convert into km
    stat_x=lonfac*(stat_lons-lonc)
    stat_y=latfac*(stat_lats-latc)
    
    # make while loop
    NN=len(stat_lons)
    check=10. # initial value so that the loop executes at least once
    count=1
    #numtries=50 # maximum number of iterations allowed
    # make G matrix and data
    G=np.matrix(np.zeros((NN,3))) # design matrix
    D = np.matrix(np.zeros((NN,1))) # data
    d = np.matrix(np.zeros((NN,1))) # data
    while check>tol and count<numtries:
#        mod_old=mod
        # save 1st solution for comparison later
        toe=answer[0]
        xoe=answer[1]
        yoe=answer[2]
        answer_old=answer
        # make G matrix and data from new initial guess
        for i in range(0,NN):
            D[i]=np.sqrt((stat_x[i]-xoe)**2.+(stat_y[i]-yoe)**2.+(stat_els[i]+10.)**2.) # distance vector
            d[i]=stat_t[i]-toe-1./alpha*D[i] # data vector
            if D[i]==0:
                print(count,i)
            # define temporary constant which is same for each row
            a=1./alpha*1./D[i]
            G[i,0]=1.
            G[i,1]=a*(xoe-stat_x[i])
            G[i,2]=a*(yoe-stat_y[i])
        # Compute least squares solution
        GG=np.dot(np.transpose(G),G)
        GGinv=np.linalg.inv(GG)
        GGG=np.dot(GGinv,np.transpose(G))
        mod=np.dot(GGG,d)
        
        # correct epicenter with model data
        answer=np.asarray(mod).reshape(-1)+answer_old # convert to array from column vector
        # controls for exiting loop
        transmod=np.abs(mod) 
        transmod[1:2]=0.01*transmod[1:2]
        check=np.max(transmod)
    
    	# get epicenter in latitude and longitude for map
        answer_lon=answer[1]/lonfac+lonc
        answer_lat=answer[2]/latfac+latc
        lons[count]=answer_lon
        lats[count]=answer_lat
        if plotq:
            # print map for this iterate
            makemap(stat_lats,stat_lons,stat_names,lonc,latc,lon_delt,lat_delt,
                    True,count,lons,lats)
            plt.show()
            plt.pause(0.5)
        if (check<tol):
            print("Solution converged")
            plt.pause(1.)
            plt.close()
        count=count+1  
    return answer, GGinv, count