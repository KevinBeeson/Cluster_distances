
# -*- coding: utf-8 -*-
import numpy as np
from astropy.io.votable import parse
import csv
from astropy.table import Table,vstack
import mpmath as mp
import math
import os
import emcee
from multiprocessing import Pool
from zero_point import zpt
import scipy
from functools import partial
mp.dps=20
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
def integrand(r,w,sigW,L,rc,sc):
    return prior_distribution(r,L)*likelihood(r,w,sigW,L)*cluster_prior(r,rc,sc)
def prior_distribution(r,L,alpha,beta):
    """
    geometric distance prior
    Parameters
    ----------
    r : float
        distance from earth
    L,alpha,beta : float
        3 types of length scales.

    Returns
    -------
    TYPE
        unormalized prior distribution.

    """
    #un realisitic distance
    if r<0:
        return 0
    return r**beta*mp.exp(-(r/L)**alpha)
def likelihood(r,w,sigW,zero_point):
    """
    returns the likelihood for a certain star to be at point r given L

    Parameters
    ----------
    r : float
        distance.
    w : float
        parallax.
    sigW : float
        error in parallax.
    zero_point : float
        zero_point_error.

    Returns
    -------
    float
        likelihood.

    """
    if r<0:
        return 0
    return mp.exp(-1/(2*sigW**2)*(w-zero_point-1/r)**2)/(mp.sqrt(2*mp.pi*sigW**2))
def cluster_prior(r,rc,sc):
    """
    cluster prior without the gamma distribution

    Parameters
    ----------
    r : float
        distance.
    rc : float
        distance to the centre of the cluster.
    sc : float
        cluster spread.

    Returns
    -------
    float
        unormalized probability.

    """
    return mp.exp(-(r-rc)**2/(2*sc**2))/(mp.sqrt(2*mp.pi*sc**2))
def gammaSS(x,shape,scale):
    """
    gamma distribution for cluster 

    Parameters
    ----------
    x : float
        cluster spread.
    shape ,scale: float
        gamma distribution shapes.
    -------
    TYPE
        probability of the cluster paramaters being that.

    """
    return x**(shape-1)*np.exp(-x/scale)/(scale**shape*math.gamma(shape))
def exponent_negative(r):
    """
    exponent part of the conditional probability times -1 used to get the maximum of the function

    Parameters
    ----------
    r : distance
        DESCRIPTION.

    Returns
    -------
    float
        e.

    """
    return (-(r/L)**alpha -1/(2*sigw**2)*(w-zero_point-1/r)**2-1/(2*sc**2)*(r-rc)**2)*-1

def logVectorProb(r):
    if any(dist<0.0 for dist in r)==True:
        return -np.inf
    probability1=(w+0.000029-1/r)**2/(2*sigW**2)+r/L
    probability1=np.sum(probability1)
    probabilityNormal=2*(np.sum(np.log(r/(2*np.pi*L**(3/2))))-len(r)*np.log(sc*rc))
    rcArray=np.ones(len(r))*rc
    scArray=np.ones(len(r))*sc
    
    probability2=np.dot((r-rcArray)**2,1/(scArray)**2)/2
    return -probability1+probabilityNormal-probability2

def new_splitter(func,rc,sc,bins):
    r_lim=np.linspace(rc-sc*15,rc+sc*15,num=50)
    probabilities=[func(x) for x in r_lim]
    sum_probabilities=np.sum(probabilities)
    index_max=probabilities.index(max(probabilities))
    if abs(index_max-25)>0:
        r_lim=np.linspace(r_lim[index_max]-sc*15,r_lim[index_max]+sc*15,num=50)
        probabilities=[func(x) for x in r_lim]
        sum_probabilities=np.sum(probabilities)
        index_max=probabilities.index(max(probabilities))        
    
    test_probabilities=0
    index_shift=0
    while test_probabilities/sum_probabilities<0.90:
        test_probabilities+=(probabilities[index_max-index_shift]+probabilities[index_max+index_shift]+probabilities[index_max-1-index_shift]+probabilities[index_max+1+index_shift])*0.5
        index_shift+=1
    
    limits=np.linspace(r_lim[index_max-index_shift+1],r_lim[index_max+index_shift-1],num=bins-1)
    limits=np.append(-mp.inf,limits)
    limits=np.append(limits,mp.inf)
    return limits 
#gives the bins where the integral will integrate 
def splitter(func, bins,maxPoint):
    bins-=2
    if bins%2==0:
        return 0
    
    diff2=mp.diff(func,maxPoint,n=2)
    diff1=mp.diff(func,maxPoint,n=1)
    maxProb=func(maxPoint)
    probChange=maxProb/bins

    roots=np.roots([+1/2*diff2,diff1,probChange])
    max_change=roots[0]
    if max_change<0:
        max_change=roots[1]
        
    
    
    
    lowPoints= [0 for x in range(int((bins+1)/2))]  
    lowPoints[0]=maxPoint-max_change   
    for x in range(int((bins-1)/2)):
        diff1=mp.diff(func,lowPoints[x],n=1)
        change=probChange/diff1
        lowPoints[1+x]=lowPoints[x]-change
    lowPoints=lowPoints[::-1]    
        
    highPoints= [0 for x in range(int((bins+1)/2))]  
    highPoints[0]=maxPoint+max_change   
    for x in range(int((bins-1)/2)):
        diff1=mp.diff(func,highPoints[x],n=1)
        change=probChange/diff1
        highPoints[1+x]=highPoints[x]-change    
    
    points=np.hstack((-mp.inf,lowPoints,highPoints,+mp.inf))    
    return points
#Individual stars loop
def individual_star_loop(star,cluster):
        rc=cluster[0]
        sc=cluster[1]
        w=star[0]/1000
        sigW=star[1]/1000           
        L=float(star[3])
        alpha=float(star[4])
        beta=float(star[5])
        zero_point=float(star[2])
        f=lambda r: prior_distribution(r,L,alpha,beta)*likelihood(r,w,sigW,zero_point)*cluster_prior(r,rc,sc)
        bins=3        
        # #gets the bins of the integral
        # max_point=scipy.optimize.minimize_scalar(exponent_negative,bounds=[rc-sc*5,rc+sc*5],method='bounded')
        limits=new_splitter(f,rc,sc,bins)
        newIntegral=0.0
        for y in range(len(limits)-1):
            currentInt=mp.quad(f,[limits[y],limits[y+1]])
            newIntegral+=currentInt
        while newIntegral==0.0:
            if bins>50:
                return -np.inf
            newIntegral=0.0
            bins=bins*2+1
            limits=new_splitter(f,rc,sc,bins)
            for y in range(len(limits)-1):
                newIntegral+=mp.quad(f,[limits[y],limits[y+1]])
        return newIntegral
        # tempProbability+=mp.log(newIntegral)
#combining the probabilities
def totalProbability(walkers,ncpu=12):
    try:
        rc=walkers[0]
        sc=walkers[1]
        
        tempProbability=0.0
        if rc<=0.0 or sc<=0.0:
            return -np.inf
        with Pool(processes=ncpu ) as pool:
            probabilities=pool.map(partial(individual_star_loop,cluster=walkers),data_short_iteration)
        # for x in data_short:
        #     #gets the data from each individual stars
        #     w=x['parallax']/1000
        #     sigW=x['parallax_error']/1000           
        #     L=float(x['L'])
        #     alpha=float(x['alpha'])
        #     beta=float(x['beta'])
        #     zero_point=float(x['zero_point'])
        #     f=lambda r: prior_distribution(r,L,alpha,beta)*likelihood(r,w,sigW,zero_point)*cluster_prior(r,rc,sc)
        #     bins=3        
        #     # #gets the bins of the integral
        #     # max_point=scipy.optimize.minimize_scalar(exponent_negative,bounds=[rc-sc*5,rc+sc*5],method='bounded')
        #     limits=new_splitter(f,rc,sc,bins)
        #     newIntegral=0.0
        #     for y in range(len(limits)-1):
        #         currentInt=mp.quad(f,[limits[y],limits[y+1]])
        #         newIntegral+=currentInt
        #     while newIntegral==0.0:
        #         if bins>50:
        #             return -np.inf
        #         newIntegral=0.0
        #         bins=bins*2+1
        #         limits=new_splitter(f,rc,sc,bins)
        #         for y in range(len(limits)-1):
        #             newIntegral+=mp.quad(f,[limits[y],limits[y+1]])
            # tempProbability+=mp.log(newIntegral)
        probabilities=[float(x) for x in probabilities]
        tempProbability=np.sum(np.log(probabilities))
        tempProbability+=np.log(gammaSS(sc/1000,shape=2,scale=skew))
        return tempProbability
    except IndexError:
        return -np.inf
summary=np.array(list(csv.reader(open('prior_summary.csv', 'rt'), delimiter=',')))  
summary=np.array(summary[1:],dtype=float) 
heal_pix_table=np.array(list(csv.reader(open('HEALpixel_level5_radec_longlat_coordinates.csv', 'rt'), delimiter=',')))   
heal_pix_table=np.array(heal_pix_table[1:],dtype=float)

#shows the list of stars and gets its parallax and sigW
list_stars = list(csv.reader(open('new_criteria_gaia_distances2.txt', 'rt'), delimiter='\t'))
list_stars = list_stars[1::]
list_stars = np.array(list_stars)
# print('Input which cluster you would like to see from the list bellow')
target_list = list(csv.reader(open('my_targets.txt', 'rt'), delimiter='\t'))
# print(list_stars[:,0])
# name=input()

name='upper_Scorpius'
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

if name=='ASCC_16':
	data_short=vstack([x for x in data_short if abs(x['r_med_geo']-430)<2*(x['r_hi_geo']-x['r_lo_geo'])])
	
# dataShort2=Table(data[0],names=(data.colnames),meta={'name':'shorter2 star data'},masked=True)
# del dataShort2[0]
# for star in range(len(data_short)):
#     if not( data_short['r_est'][star]+2*(data_short['r_hi'][star]-data_short['r_est'][star]) <float(list_stars[position,-1])*0.8-20 or data_short['r_est'][star]-2*(data_short['r_est'][star]-data_short['r_lo'][star])>float(list_stars[position,-1])*1.2+20) :
#         dataShort2.add_row(data_short[star])
# data_short=dataShort2
# priorR=3600
# priorSigR=100
shape=2
skew=0.01

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


star=data_short[0]
L=star['L']
alpha=star['alpha']
beta=star['beta']
sigw=star['parallax_error']/1000
w=star['parallax']/1000
zero_point=star['zero_point']/1000
data_short_iteration=[[x['parallax'],x['parallax_error'],x['zero_point'],x['L'],x['alpha'],x['beta']] for x in data_short]
priorRc=np.mean(data_short['r_med_geo'])
priorSc=np.mean(data_short['r_hi_geo']-data_short['r_lo_geo'])/2
totalProbability([110.71825131,1.21249576])
# nwalkers=4
# ndim=2

# # #initializes the walkers and makes sure none of the rc or sc is negative
# scWalkers=[]
# rcWalkers=[]
# while(len(scWalkers)<nwalkers):
#     rcTemp=np.random.normal(priorRc,10.0,1)
#     scTemp=np.random.normal(5,1,1)
#     if rcTemp>0 and scTemp>0:
#         rcWalkers=np.append(rcWalkers,rcTemp)
#         scWalkers=np.append(scWalkers,scTemp)
# walkers=np.column_stack((rcWalkers,scWalkers))


# filename = name+"_smaller_scale.h5"
# backend = emcee.backends.HDFBackend(filename)
# backend.reset(nwalkers, ndim)

# # totalProbability([239.386071 , 6.04066001],12
# ncpu=64
# step_iteration=5
# sampler=emcee.EnsembleSampler(nwalkers, ndim, totalProbability,backend=backend,args=[ncpu])
# autocorr=[]
# oldTau=np.inf
# labels=['rc','sc']
# # totalProbability([305,10])
# for sample in sampler.sample(walkers,iterations=1875, progress=True):
#     if sampler.iteration % step_iteration:
#             continue
#     tau=sampler.get_autocorr_time(tol=0)
#     # sampler_mean=np.mean(sampler.get_chain(flat=True)[-int(sampler.iteration*0.6*nwalkers):],axis=1)
#     # shift_temp=shift_maker(sampler_mean)
#     # print("Autocorelation Time is :",tau)
#     # if sampler.iteration>np.nanmax(tau)*7 and sampler.iteration>min_steps:
#     #     temp_data=sampler.get_chain(flat=True)[-int(sampler.iteration*0.6*nwalkers):]
#     #     print('gaussian and moving difference '+ str(how_good_gaussian_fit(temp_data)))
#     #     if how_good_gaussian_fit(temp_data)<0.005:
#     #         print('converged_stoped')
#     #         break
#     if not np.any(np.isnan(tau)):
#         autocorr.append(tau)
#         converged = np.all(tau * 17 < sampler.iteration)
        
#         print('sampler autocorelation time rc,sc ', sampler.iteration/tau)

#         converged &= np.all(np.abs(oldTau - tau) / tau < 0.03)
#         print('change in sampler autocorelation time rc,sc', np.abs(oldTau - tau) / tau )
#         print('current mean',np.mean(sampler.get_chain(flat=True),axis=0))
#         if converged and sampler.iteration>100:
#             break
#         oldTau = tau

# toSave=sampler.get_chain(flat=True)

# np.save("Cluster distance"+name,toSave)

# L=200
# w=2
# sigW=1
# normalization=integrate.quad(probability,0,+np.inf,args=(w,sigW,L))
# r =np.linspace(-100,100,num=1e+5)
# posterior=[]
# prob=[]
# for x in r:
#     posterior.append(prior_distribution(x,L))
#     prob.append(probability(x,w,sigW,L))
# posterior=np.array(posterior)/normalization[0]
# prob=np.array(prob)/normalization[0]
# plt.figure()
# plt.plot(r,np.array(posterior))
# plt.figure()
# plt.plot(r,np.array(prob))
