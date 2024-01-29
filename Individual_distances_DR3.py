import numpy as np
from astropy.io.votable import from_table, writeto,parse
import csv
from astropy.table import Table
import math as maths
import matplotlib.pyplot as plt
import emcee
from multiprocessing import Pool
import os
import corner
from scipy.stats import norm,cauchy
import matplotlib.mlab as mlab
from zero_point import zpt

os.environ["OMP_NUM_THREADS"] = "1"
def heal_pixel_finder(ra,dec):
    """finds the heal pixel given the ra and dec returns heal pixel"""
    if isinstance(ra, (float,np.floating)):
        closest_neighbour=(np.min((abs(360-ra+heal_pix_table[:,3]),abs(ra-heal_pix_table[:,3])),axis=0))**2+(np.min((abs(360-dec+heal_pix_table[:,4]),abs(dec-heal_pix_table[:,4])),axis=0))**2
        heal_pixel=heal_pix_table[np.where(closest_neighbour==min(closest_neighbour))][0]
    else:
        closest_neighbour_ra=[(np.min((abs(360-x+heal_pix_table[:,3]),abs(x-heal_pix_table[:,3])),axis=0))**2 for x in ra]
        closest_neighbour_dec=[(np.min((abs(360-x+heal_pix_table[:,4]),abs(x-heal_pix_table[:,4])),axis=0))**2 for x in dec]
        closest_neighbour=[x+y for (x,y) in zip(closest_neighbour_ra,closest_neighbour_dec)]
        heal_pixel=[heal_pix_table[np.where(x==min(x))][0] for x in closest_neighbour]
        heal_pixel=np.array(heal_pixel)
        return heal_pixel[:,0]

    return heal_pixel[0]
def l_alpha_beta_finder(heal_pixel):
    """insert a heal pixel and will get alpha beta and L in return"""
    if isinstance(heal_pixel, (float,np.floating)):
        parameters=summary[np.where(summary[:,0]==heal_pixel)]
        return parameters[0][5:8]
    else:
        parameters=[summary[np.where(summary[:,0]==x)][0][5:8] for x in heal_pixel]
        parameters=np.vstack(parameters)
        return parameters
def integrand(r,i):
    return posD(r)*rProbability(r)*clusterProbability(r,i)
def posD(r):
    if r<0:
        return 0
    return 1/(2*L**3)*r**2*np.exp(-r/L)
def rProbability(r):
    if r<0:
        return 0
    return np.exp(-1/(2*sigW**2)*(w+0.000029-1/r)**2)*1/(np.sqrt(2*np.pi*sigW**2))
def clusterProbability(r,i,j):
    return np.exp(-(r-rc[i])**2/(2*sc[i]**2))*1/(np.sqrt(2*np.pi*sc[i]**2))
def gammaSS(x,shape,scale):
    return x**(shape-1)*np.exp(-x/scale)/(scale**shape*maths.gamma(shape))

#LogIntergrand using 
def logInterR(r):
    probability=0.0
    for i in range( np.size(probDist,1)):
        # probability+=-(1/(2*sigW**2)*(w+0.000029-1/r)**2+r/L)+2*np.log(r/(2*np.pi*sc[i]*rc[i]*np.sqrt(L**3)))
        probability+=-(1/(2*sigW**2)*(w+0.000029-1/r)**2+(r-rc[i])**2/(2*sc[i]**2)+r/L)+2*np.log(r/(2*np.pi*sc[i]*rc[i]*np.sqrt(L**3)))
    return probability

#insert r distance and get the total log probability of that distance for the star, the function handles sc and rc in a vectorial way
def vectorProbability(r):
    if any(dist<0.0 for dist in r)==True:
        return -np.inf
    probability1=(w+0.000029-1/r)**2/(2*sigW**2)+r/L
    probability1*=len(sc)
    probabilityNormal=2*(len(sc)*np.log(r/(2**(1/2)*np.pi*L**(3/2))))-probConstant
    rArray=np.ones(len(sc))*r
    probability2=np.dot((rArray-rc)**2,1/(sc)**2)/2
    return -probability1+probabilityNormal-probability2
def vectorProbability2(r,w,sigW,L,alpha,beta,zero_point):
    if r<0.0:
        return -np.inf
    exp1=-(w-zero_point-1/r)**2/(2*sigW**2)-(r/L)**alpha+beta*np.log(r)
    
    rArray=np.ones(len(sc))*r
    exp2=-(rArray-rc)**2*1/(sc)**2/2
    exp2=np.sum(np.exp(exp2))
    exp2=np.log(exp2)
    
    
    total_probability=exp1+exp2
                                 
    return total_probability
def interR(r):
    probability=0.0
    for i in range( np.size(probDist,1)):
        for j in range( np.size(probDist,1)):
            probability+=integrand(r,i,j)
    probability/=np.size(probDist,1)
    return np.log(probability)
def oldLogProbability(r):
    if r[0]>0.0:
        probability=-((w+0.000029-1/r)**2/(2.0*sigW**2)+r/L)+2.0*(np.log(abs(r)/(2**(1/2)*sigW*np.pi*L**(3/2))))
        return probability[0]     
    else:
        probability= -np.inf
        return probability     
    return probability[0]


# def lohInterR(r,w,sigW,L,ProbDist):
#     logProbability=1
#         for i in range(len(ProbDist)):
summary=np.array(list(csv.reader(open('prior_summary.csv', 'rt'), delimiter=',')))  
summary=np.array(summary[1:],dtype=float) 
heal_pix_table=np.array(list(csv.reader(open('HEALpixel_level5_radec_longlat_coordinates.csv', 'rt'), delimiter=',')))   
heal_pix_table=np.array(heal_pix_table[1:],dtype=float)
#shows the list of stars and gets its parallax and sigW
list_stars = list(csv.reader(open('new_criteria_gaia_distances2.txt', 'rt'), delimiter='\t'))
list_stars = list_stars[1::]
list_stars = np.array(list_stars)
print('Input which cluster you would like to see from the list bellow')
print(list_stars[:,0])

def main_loop(star):
    #picks a star in the data set
    nwalkers=12
    ndim=1

    maxSteps=2.5e+2   
    rwalkers=np.random.normal(star['r_med_geo'],(star['r_hi_geo']-star['r_lo_geo'])/2,nwalkers)

    maxSteps=1e+4             
    nwalkers=16
    w=star['parallax']/1000
    sigW=star['parallax_error']/1000
    L=float(star['l'])
    alpha=float(star['alpha'])
    beta=float(star['beta'])
    zero_point=float(star['zero_point'])
    ndim=1
    rwalkersStarter=rwalkers
    probabilityDifference=1
    sigmaOld=0
    muOld=0
    sigma=20
    mu=0.1
    fails=1
    while(5*sigma/mu>1 and maxSteps<6.0e+5 and fails < 6):
        np.random.seed()
        rwalkers=[]
        while(len(rwalkers)!=nwalkers):
            rTemp=np.random.normal(star['r_med_geo'],(star['r_hi_geo']-star['r_lo_geo'])/2,1)
            if rTemp>0 and vectorProbability2(rTemp,w,sigW,L,alpha,beta,zero_point)!=-np.inf:
                rwalkers=np.append(rwalkers,rTemp)
        rwalkers=[[x] for x in rwalkers]
        sampler=emcee.EnsembleSampler(nwalkers, ndim, vectorProbability2,args=[w,sigW,L,alpha,beta,zero_point])
        autocorr=[]
        index=0
        oldTau=np.inf
        for sample in sampler.sample(rwalkers,iterations=int(maxSteps), progress=False):
                if sampler.iteration % 100:
                    continue
                tau=sampler.get_autocorr_time(tol=0)
                # print(tau)
                if np.isnan(tau):
                    break
                autocorr.append(tau)
                converged = np.all(tau * 100*fails < sampler.iteration)
                converged &= np.all(np.abs(oldTau - tau) / tau < 0.01)
                if converged:
                    break
                oldTau = tau
        sample1=sampler.get_chain(flat=True)
        if sampler.iteration==maxSteps:
            maxSteps*=2
        sample1=sorted(sample1)
        for number in range(len(sample1)):
            sample1[number]=sample1[number][0]
        (mu, sigma) = norm.fit(sample1)
        print(5*sigma/mu)
        
        # plt.figure()
        # plt.show()
        print(probabilityDifference,"probability Difference")
        if 5*sigma/mu>.3:
            fails+=1


    lowDistance1=sample1[int(np.round(len(sample1)*15.8/100))]
    meanDistance1=sample1[int(np.round(len(sample1)/2))]
    highDistance1=sample1[int(np.round(len(sample1)*84.2/100))] 
    return mu,sigma
    # newDistance=np.append(newDistance,mu)
    # newSigDistance=np.append(newSigDistance,sigma)
ncpu=64 
my_targets=np.array(list(csv.reader(open('my_targets.txt', 'rt'), delimiter=',')))  
for name in my_targets:
    name=name[0]
    print(name)
    position=int(np.where(list_stars==name)[0])
    
    #opens the data for the specific star you have chosen
    votable = parse(name+"_dr3.xml")
    data=votable.get_first_table().to_table(use_names_over_ids=True)
    # data=data[:10]
    # distanceData=list(csv.reader(open(name+'r_len.txt','rt'),delimiter=' '))
    # distanceData=np.array(distanceData)
    # data['r_len']=distanceData[:,2]
     
    #shorterns the data to the data around a specified pmra and pmdec
    data_short=Table(data[0],names=(data.colnames),meta={'name':'shorter star data'},masked=True)
    del data_short[0]
    for star in range(len(data)):
        if data['pmra'][star] >float(list_stars[position,3])-float(list_stars[position,4]) and data['pmra'][star] <float(list_stars[position,3])+float(list_stars[position,4]):
            if data['pmdec'][star] >float(list_stars[position,5])-float(list_stars[position,6]) and data['pmdec'][star] <float(list_stars[position,5])+float(list_stars[position,6]):
                data_short.add_row(data[star])
    data_short2=Table(data[0],names=(data.colnames),meta={'name':'shorter2 star data'},masked=True)
    
    l_alpha_beta=l_alpha_beta_finder(heal_pixel_finder(data_short['ra'],data_short['dec']))
    data_short['L']=l_alpha_beta[:,0]
    data_short['alpha']=l_alpha_beta[:,1]
    data_short['beta']=l_alpha_beta[:,2]
    phot_g_mean_mag=data_short['phot_g_mean_mag']
    nu_eff_used_in_astrometry=data_short['nu_eff_used_in_astrometry']
    pseudocolour=data_short['pseudocolour']
    ecl_lat=data_short['ecl_lat']
    astrometric_params_solved=data_short['astrometric_params_solved']
    zpt.load_tables()
    
    
    zero_point=zpt.get_zpt(phot_g_mean_mag, nu_eff_used_in_astrometry, pseudocolour, ecl_lat, astrometric_params_solved)/1000
    data_short['zero_point']=zero_point
    
    #imports SC and RC data
    cluster_distance=emcee.backends.HDFBackend(name+"_try_3.h5") 
    cluster_distance=cluster_distance.get_chain(flat=True)[200:]
    rc=cluster_distance[:,0]
    sc=cluster_distance[:,1]     
    #Constant to use in vectorial integral, it stays the same so I compute it first
    probConstant=2*np.sum(np.log(sc)+np.log(rc))
    probDist=np.array([sc,rc])
    
    
                
    newDistance=[]
    newSigDistance=[]
    current=0
    #individual star distance pdf with it being in a cluster
    print("there are", len(data_short)," many stars")
    badStars=0
    averageProbabilityDifference=0
    
    pool=Pool(processes=ncpu)
    results=pool.map(main_loop,data_short)
    results=np.array(results)
    newDistance=results[:,0]
    newSigDistance=results[:,1]
    data_short['r_est_cluster']=newDistance
    data_short['r_sigma_cluster']=newSigDistance
    
    votable = from_table(data_short)
    
    writeto(votable,name+"_dr3_distances3.xml")
# print(name,badStars)

# print(main_part(data[0]))

# binOld=np.linspace(min(sample1),max(sample1),num=1000)
# for number in range(len(sample1)):
#     sample1[number]=sample1[number][0]
# inds=np.digitize(np.transpose(sample1),np.transpose(binOld)[0])
# fig=plt.figure()
# ax = fig.add_subplot(111)
# inds=inds*max(sample1)/max(inds)
# n, bins, patches =ax.hist(sample1,density=True,bins=1000,histtype='step',label='PDF star using Cluster distance')

# locs, labels = plt.xticks()
# labels=[float(item)*(max(sample)/max(inds)) for item in locs]
# plt.xticks(locs, labels)

# data_short['r_est_cluster']=newDistance
# data_short['r_sigma_cluster']=newSigDistance

#votable = from_table(data_short)

#writeto(votable,name+"_dr3_distances2.xml")
#print(name,badStars)

#individual star distance without it being in a cluster
# nwalkers=8
# rwalkers=np.random.normal(1/w,1/sigW,nwalkers)
# rwalkers=[]
# while len(rwalkers)!=nwalkers:
#     test=np.random.normal(1/w,1/sigW,1)
#     if test>0.0:
#         rwalkers=np.append(rwalkers,test)
# with Pool() as pool:
#     sampler=emcee.EnsembleSampler(nwalkers,ndim,oldLogProbability,pool=pool)
#     state=sampler.run_mcmc(np.transpose([rwalkers]), 100, progress=True)
#     sampler.reset()
#     sampler.run_mcmc(state, maxSteps, progress=True)
# sample=sampler.get_chain(flat=True)
# sample=sorted(sample)
# lowDistance=sample[int(np.round(len(sample)*15.8/100))]
# meanDistance=sample[int(np.round(len(sample)/2))]
# highDistance=sample[int(np.round(len(sample)*84.2/100))]
# binOld=np.linspace(min(sample),max(sample),num=1000)
# for number in range(len(sample)):
#     sample[number]=sample[number][0]
# inds=np.digitize(np.transpose(sample),np.transpose(binOld)[0])
# # fig=plt.figure()
# # ax = fig.add_subplot(111)
# inds=inds*max(sample)/max(inds)
# n, bins, patches =ax.hist(sample,density=True,bins=1000,histtype='step',label='PDF of star without using cluster distance')
# locs, labels = plt.xticks()
# labels=[float(item)*(max(sample)/max(inds)) for item in locs]
# plt.xticks(locs, labels)



# plt.axvline(x=np.median(rc),label='median Rc',color='r')
# ax.set_yscale('log')
# plt.legend(loc='best')
# plt.show()
# coeff=[1/L,-2.0,(w+0.000029)/sigW**2,-1/sigW**2]
# roots=np.roots(coeff)


# fig=plt.figure()
# ax = fig.add_subplot(111)
# n, bins, patches =ax.hist(sample1,density=True,bins=1000,histtype='step',label='PDF star using Cluster distance')
# toSave=np.array([n, bins])
# (mu, sigma) = norm.fit(sample1)
# y=1/(sigma*np.sqrt(2*np.pi))*np.exp(-(bins-mu)**2/(2*sigma**2))
# ax.plot(bins,y)


# np.save(str(star['source_id'])+' distance',toSave)
