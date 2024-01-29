#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep  9 14:30:59 2022

@author: kevin
"""
import matplotlib.pyplot as plt
from astropy.io.votable import parse,from_table,writeto
from astropy.table import Table, vstack
from functools import partial
import  numpy as np 
import requests
import re
import csv
from os.path import exists
from pathlib import Path
from matplotlib.collections import LineCollection
from scipy.interpolate import interp1d
from multiprocessing import Pool
from scipy.optimize import curve_fit
import scipy.stats as sps
import copy
import random
from time import sleep
from numba import jit

@jit(nopython=True,cache=True)
def normal(x,mu,sigma):
    return 1/(sigma*(2*3.141592653589793)**0.5)*np.exp(-1/2*((x-mu)/sigma)**2)
def get_isochrone(low_age,high_age,low_metalicty,high_metalicty,a_v,age_spacing=0.01,metalicty_spacing=0.01,iso_type='gaia'):
    """
    asks the padava server for an isochrone or a list of isochrones

    Parameters
    ----------
    low_age : float
        DESCRIPTION.
    high_age : float
        DESCRIPTION.
    low_metalicty : float
        DESCRIPTION.
    high_metalicty : float
        DESCRIPTION.
    a_v : float
        extinction.
    age_spacing : float, optional
        DESCRIPTION. The default is 0.01.
    metalicty_spacing : float, optional
        DESCRIPTION. The default is 0.01.

    an array of the isochrone
    -------
    None.

    """
    #if metallicity is given as a log ([M/H]) convert it here:


#    mass=[]
#    label=[]
#    mags=[]
#    imf=[]
#    teff=[]
#    logg=[]

    #parameters other than default
    d={
    'track_parsec': 'parsec_CAF09_v1.2S',
    'track_colibri':'parsec_CAF09_v1.2S_S35',
    'track_postagb':'no',
    'n_inTPC': 10,
    'eta_reimers': 0.2,
     # 'photsys_file': 'tab_mag_odfnew/tab_mag_gaiaDR2weiler.dat',
     'photsys_file':'tab_mag_odfnew/tab_mag_gaiaEDR3.dat',
   #   'photsys_file':'tab_mag_odfnew/tab_mag_gaia.dat',
  # 'photsys_file':'tab_mag_odfnew/tab_mag_2mass_spitzer_wise.dat',
    # 'photsys_file':'tab_mag_odfnew/tab_mag_panstarrs1.dat',
     #'photsys_file':'tab_mag_odfnew/tab_mag_gaiaDR2.dat',
    'photsys_version': 'OBC',
    'dust_sourceM': 'nodustM',
    'dust_sourceC': 'nodustC',
    'extinction_av': a_v,
    'extinction_coeff':'constant',
    'extinction_curve':'cardelli',   
    'imf_file': 'tab_imf/imf_kroupa_orig.dat',
    'isoc_isagelog':'1',
    'isoc_lagelow':low_age,
    'isoc_lageupp':high_age,
    'isoc_dlage':age_spacing, #steps ages
    'isoc_ismetlog':'1',
    'isoc_metlow':low_metalicty,
    'isoc_metupp':high_metalicty,
    'isoc_dmet':metalicty_spacing, #steps M/H
    'output_kind': 0,
    'submit_form': 'Submit'}
    if iso_type=='gaia':
        d['photsys_file']='tab_mag_odfnew/tab_mag_gaiaEDR3.dat'
    elif iso_type=='panstarrs':
        d['photsys_file']='tab_mag_odfnew/tab_mag_panstarrs1.dat'
    elif iso_type=='allwise':
        d['photsys_file']='tab_mag_odfnew/tab_mag_2mass_spitzer_wise.dat'
    #Check if we already downloaded this isochrone.
    #Isochrones are saved as txt files and the filename is the hash of the dictionary values.
    webserver = 'http://stev.oapd.inaf.it'
    c = requests.get(webserver + '/cgi-bin/cmd_3.7', params=d).text
#    print(c)
    aa = re.compile('output\d+')
    fname = aa.findall(c)
    if len(fname) > 0:
        url = '{0}/tmp/{1}.dat'.format(webserver, fname[0])
        #print url
        r = requests.get(url).text
        print('successfully gotten the ischrone')
    return r
    
def get_and_save_ishochrones(name,low_age,high_age,low_metalicty,high_metalicty,a_v,age_spacing=0.01,metalicty_spacing=0.01,special=False,iso_type='gaia'):
    """
    You can ask the function for a isochrone it will check if we already have it saved and if not it will ask from the padava servers

    Parameters
    ----------
    name : str
        name of the cluster.
    low_age : TYPE
        DESCRIPTION.
    high_age : TYPE
        DESCRIPTION.
    low_metalicty : TYPE
        DESCRIPTION.
    high_metalicty : TYPE
        DESCRIPTION.
    a_v : TYPE
        DESCRIPTION.
    age_spacing : TYPE, optional
        DESCRIPTION. The default is 0.01.
    metalicty_spacing : TYPE, optional
        DESCRIPTION. The default is 0.01.
    special : TYPE, optional
        DESCRIPTION. The default is False.

    Returns
    -------
    None.

    """
    Path('ALL_ISO').mkdir(parents=True,exist_ok=True)
    #gets the path name of the isochrone and checks if its already been downloaded

    if (low_age==high_age or high_age==None) and (low_metalicty==high_metalicty or high_metalicty==None):
        if special:
            path_name='ALL_ISO/'+name+'_'+special+'_age_'+str(low_age)+'_metalicty_'+str(low_metalicty)+'_extinction_'+str(a_v)
        else:
            path_name='ALL_ISO/age_'+str(low_age)+'_metalicty_'+str(low_metalicty)+'_extinction_'+str(a_v)
    elif low_metalicty==high_metalicty or high_metalicty==None:
        path_name='ALL_ISO/age_'+str(low_age)+'_to_'+str(high_age)+'_metalicty_'+str(low_metalicty)+'_extinction_'+str(a_v)

    elif low_age==high_age or high_age==None:
        path_name='ALL_ISO/age_'+str(low_age)+'_metalicty_'+str(low_metalicty)+'_to_'+str(high_metalicty)+'_extinction_'+str(a_v)


    else:
        
        if special:
            path_name='ALL_ISO/'+name+'_'+special+'_mixed'
        else:
            path_name='ALL_ISO/mixed'   
        if iso_type!='gaia':
            path_name+='_'+iso_type
        path_name+='.txt'
        #asks for the isochorne
        data=get_isochrone(low_age,high_age,low_metalicty,high_metalicty,a_v,age_spacing,metalicty_spacing,iso_type=iso_type)
        data=data.split('\n')
        header=data[11]
        data=data[0:len(data)-1]
        ammendedData=''
        for x in data:
            if x[0]!='#':
                numbers=np.array(x.split(' '))
                numbers=numbers[numbers!='']
                ammendedData+=' '.join(numbers)+'\n'
            else:
                numbers=np.array(x.split(' '))
                numbers=numbers[numbers!='']
                ammendedData+=' '.join(numbers)+'\n'
                
        ammendedData=ammendedData.split('\n')
        ammendedData=[i.split(' ') for i in ammendedData]
        ammendedData=ammendedData[0:len(ammendedData)-1]
        # ammendedData=np.array(ammendedData)
    
        toWrite=''
        for x in ammendedData:
            line=''
            for y in x:
                line+=str(y)+' '
            line=line[:-1]
            toWrite+=line+'\n'    
        file=open(path_name,'w')
        file.write(toWrite)
    if iso_type!='gaia':
        path_name+='_'+iso_type
    path_name+='.txt'
    if not exists(path_name):
        
        data=get_isochrone(low_age,high_age,low_metalicty,high_metalicty,a_v,iso_type=iso_type)
        data=data.split('\n')
        header=data[11]
        data=data[0:len(data)-1]
        ammendedData=''
        for x in data:
            if x[0]!='#':
                numbers=np.array(x.split(' '))
                numbers=numbers[numbers!='']
                ammendedData+=' '.join(numbers)+'\n'
            else:
                numbers=np.array(x.split(' '))
                numbers=numbers[numbers!='']
                ammendedData+=' '.join(numbers)+'\n'
                
        ammendedData=ammendedData.split('\n')
        ammendedData=[i.split(' ') for i in ammendedData]
        ammendedData=ammendedData[0:len(ammendedData)-1]
        # ammendedData=np.array(ammendedData)
    
        toWrite=''
        for x in ammendedData:
            line=''
            for y in x:
                line+=str(y)+' '
            line=line[:-1]
            toWrite+=line+'\n'    
        file=open(path_name,'w')
        file.write(toWrite)

def difference_iso(iso_1,iso_2):
    """
    Gets the difference in two iso's the iso needs to be in the cropped form. Will take into account number density

    Parameters
    ----------
    iso_1 : narray
        DESCRIPTION.
    iso_2 : TYPE
        DESCRIPTION.

    Returns
    -------
    diff : TYPE
        DESCRIPTION.

    """
    if isinstance(iso_2,np.ndarray):
        number_density_1=np.subtract(iso_1[1:,1],iso_1[:-1,1])
        number_density_2=np.subtract(iso_2[1:,1],iso_2[:-1,1])
        middle_point_1=np.add(iso_1[:-1,:],iso_1[1:,:])/2
        middle_point_2=np.add(iso_2[:-1,:],iso_2[1:,:])/2
        number_density_1_temp=[]
        middle_point_1_temp=[]
        for x,y in zip(number_density_1,middle_point_1):
            if x!=0.0:
                number_density_1_temp.append(x)
                middle_point_1_temp.append(y)
        middle_point_1_temp=np.vstack(middle_point_1_temp)
        number_density_2_temp=[]
        middle_point_2_temp=[]
        for x,y in zip(number_density_2,middle_point_2):
            if x!=0.0:
                number_density_2_temp.append(x)
                middle_point_2_temp.append(y)
        middle_point_2_temp=np.vstack(middle_point_2_temp)
        middle_point_2=(middle_point_2_temp.T*number_density_2_temp).T
        middle_point_1=(middle_point_1_temp.T*number_density_1_temp).T
    
        iso_matrix=(np.ones((len(middle_point_1_temp),len(middle_point_2_temp),3))*middle_point_2_temp[:,3:6]).T
        diff=np.sum([np.min(abs(x-y),axis=0) for (x,y) in zip(iso_matrix,middle_point_1_temp[:,3:6].T)])
        return diff
    else:
        diff=[]
        for array in iso_2:
            diff.append(difference_iso(iso_1,array))
        return diff
#interpolates linearly with mass (not log mass) as the basis.
def interpolate(oldIsochrone,scale=5):   
    mass=oldIsochrone[:,0]
    #changes log temperature to temp.
    oldIsochrone[:,2]=10**(oldIsochrone[:,2])

    
    newIso=[[None for y in range(len(oldIsochrone[0]))] for x in range((len(oldIsochrone)-1)*scale)]
    newIso=np.array(newIso)
    newIso[:,0]=np.hstack(np.transpose(np.linspace(mass[:-1],mass[1:],num=scale,endpoint=True)))
    newIso=np.array(newIso,dtype=float)
    for x in range(len(oldIsochrone[0])-1):
        f = interp1d(mass,oldIsochrone[:,x+1])
        newIso[:,x+1]=f(newIso[:,0]) 
    newIso[:,2]=np.log10(newIso[:,2])

    return np.array(newIso,dtype=float)
def iso_reader(name,low_age,low_metalicty,a_v,special=False,high_age=None, high_metalicty=None,iso_type='gaia'):
    if (low_age==high_age or high_age==None) and (low_metalicty==high_metalicty or high_metalicty==None):
        if special:
            path_name='ALL_ISO/'+name+'_'+special+'_age_'+str(low_age)+'_metalicty_'+str(low_metalicty)+'_extinction_'+str(a_v)
        else:
            path_name='ALL_ISO/age_'+str(low_age)+'_metalicty_'+str(low_metalicty)+'_extinction_'+str(a_v)
        if iso_type!='gaia':
            path_name+='_'+iso_type
        path_name+='.txt'
        iso = list(csv.reader(open(path_name, 'rt'), delimiter=' '))
        
        iso=iso[13:-1]
        iso=np.array(iso,dtype=float)
        #crops the isochrone down to initial ini mass, int_imf, log_T, G, G_Bp_bright, G_BP_faint, G_RP ,log_g , cur mass
        if iso_type=='gaia':
            iso=np.column_stack((iso[:,3],iso[:,4],iso[:,7],iso[:,25],iso[:,26],iso[:,27],iso[:,8],iso[:,5]))
            
                
            iso=[x for x in iso if x[3]<15]
        return [low_age,low_metalicty],np.vstack(iso)
    elif low_metalicty==high_metalicty or high_metalicty==None:
        path_name='ALL_ISO/age_'+str(low_age)+'_to_'+str(high_age)+'_metalicty_'+str(low_metalicty)+'_extinction_'+str(a_v)+'.txt'

    elif low_age==high_age or high_age==None:
        path_name='ALL_ISO/age_'+str(low_age)+'_metalicty_'+str(low_metalicty)+'_to_'+str(high_metalicty)+'_extinction_'+str(a_v)+'.txt'

    else:
        if special:
            path_name='ALL_ISO/'+name+'_'+special+'_mixed.txt'
        else:
            path_name='ALL_ISO/mixed.txt'
    iso = list(csv.reader(open(path_name, 'rt'), delimiter=' '))
    iso_full=[]
    iso_parameters=[]
    iso_temp=[]
    
    for x in iso:
        if x[0]!='#':
            iso_temp.append(x)
        else:
            if len(iso_temp):
                iso_temp=np.array(iso_temp,dtype=float)
                age_temp=iso_temp[0][2]
                metalicity_temp=iso_temp[0][1]

                iso_temp=np.column_stack((iso_temp[:,3],iso_temp[:,4],iso_temp[:,7],iso_temp[:,25],iso_temp[:,26],iso_temp[:,27],iso_temp[:,8],iso_temp[:,5]))
                iso_temp=[x for x in iso_temp if x[3]<15]

                iso_full.append(np.vstack(iso_temp))
                iso_parameters.append([age_temp,metalicity_temp])
                iso_temp=[]
    return iso_parameters,iso_full
class isochrone:
    def __init__(self,name,age,metalicity,extinction,special=False,interpolate_scale=10,high_age=None,high_metalicity=None,limits=None,iso_type='gaia'):
        self.name=name
        self.special=special
        self.extinction=extinction
        self.iso_type=iso_type
        get_and_save_ishochrones(name,age,high_age,metalicity,high_metalicity,extinction,special=special,iso_type=iso_type)
        self.parameters,self.isochrone=iso_reader(name, age, metalicity, extinction,special=special,high_age=high_age,high_metalicty=high_metalicity,iso_type=iso_type)
        if iso_type=='allwise':
            iso_temp=isochrone(name,age,metalicity,extinction,special,interpolate_scale,high_age,high_metalicity,limits,iso_type='gaia')
            self.iso_gaia=iso_temp.isochrone
        if limits:
            if len(limits)!=2:
                raise ValueError("wrong shape for limits")
            self.limit_g=limits
        if limits:
            if isinstance(self.isochrone, np.ndarray):
                self.isochrone=np.vstack([x for x in self.isochrone if x[3]-limits[0]>0 and x[3]-limits[1]<0])  

            else:
                iso_temp=[]
                for x in self.isochrone:
                    iso_temp.append(np.vstack([x for x in isochrone if x[3]-limits[0]>0 and x[3]-limits[1]<0]))
                self.isochrone=iso_temp
        if isinstance(self.isochrone, np.ndarray) and iso_type=='gaia' and interpolate_scale:
            self.isochrone=interpolate(self.isochrone,interpolate_scale)
        elif iso_type=='gaia'and interpolate_scale:
            iso_temp=[]
            for x in self.isochrone:
                iso_temp.append(interpolate(x,interpolate_scale))
            self.isochrone=iso_temp
    def plot(self):
        plt.figure()
        if isinstance(self.isochrone,np.ndarray):
            plt.plot(self.isochrone[:,5]-self.isochrone[:,4],self.isochrone[:,3])
        else:
            for (x,label) in zip(self.isochrone,self.parameters):
                temp_label='age '+ str(label[0]) +' metalicty '+str(label[1])
            
                plt.plot(x[:,4]-x[:,5],x[:,3],label=temp_label)
                plt.legend(loc='best')
        plt.xlabel('bp_rp')
        plt.ylabel('G')
        plt.gca().invert_yaxis()
def linear_fit(x, A, B): # this is your 'straight line' y=f(x)
    return A*x + B
class all_isochrones:
    def __init__(self,numbers,iso_bestfit,iso_low,iso_high,limits=None,gaussian=False):
        self.name=iso_bestfit.name
        self.iso_type=iso_bestfit.iso_type
        self.iso_bestfit=iso_bestfit
        self.iso_low=iso_low
        self.iso_high=iso_high
        self.extinction=iso_bestfit.extinction
        self.limits=limits
        if numbers!=0:
            self.difference_low=difference_iso(self.iso_low.isochrone,self.iso_bestfit.isochrone)
            self.difference_high=difference_iso(self.iso_high.isochrone,self.iso_bestfit.isochrone)
            if gaussian:
                metalicty_wanted=[]
                while len(metalicty_wanted)<numbers:
                    metalicty_temp=np.random.normal(iso_bestfit.parameters[1],abs(iso_bestfit.parameters[1]-iso_low.parameters[1]),1)
                    if metalicty_temp<iso_bestfit.parameters[1]:
                        metalicty_wanted.append(metalicty_temp)
                while len(metalicty_wanted)<numbers*2:
                    metalicty_temp=np.random.normal(iso_bestfit.parameters[1],abs(iso_bestfit.parameters[1]-iso_high.parameters[1]),1)
                    if metalicty_temp>iso_bestfit.parameters[1] and metalicty_temp<1.0:
                        metalicty_wanted.append(metalicty_temp)
                age_wanted=[]
                while len(age_wanted)<numbers:
                    age_temp=np.random.normal(iso_bestfit.parameters[0],abs(iso_bestfit.parameters[0]-iso_low.parameters[0]),1)
                    if age_temp<iso_bestfit.parameters[0]:
                        age_wanted.append(age_temp)
                while len(age_wanted)<numbers*2:
                    age_temp=np.random.normal(iso_bestfit.parameters[0],abs(iso_bestfit.parameters[0]-iso_high.parameters[0]),1)
                    if age_temp>iso_bestfit.parameters[0]:
                        age_wanted.append(age_temp)
                extinction_wanted=[]
                while len(extinction_wanted)<numbers:
                    extinction_temp=np.random.normal(iso_bestfit.extinction,abs(iso_bestfit.extinction-iso_low.extinction),1)
                    if extinction_temp<iso_bestfit.extinction:
                        extinction_wanted.append(extinction_temp)
                while len(extinction_wanted)<numbers*2:
                    extinction_temp=np.random.normal(iso_bestfit.extinction,abs(iso_bestfit.extinction-iso_high.extinction),1)
                    if extinction_temp>iso_bestfit.extinction:
                        extinction_wanted.append(extinction_temp)
                random.shuffle(extinction_wanted)
                random.shuffle(age_wanted)

                self.parameters=parameters=np.hstack((age_wanted,metalicty_wanted,extinction_wanted))
                isochrones_list=[]
                current_isochrone=0
                results=[]
                # for param in parameters:
                #     results.append(iso_pool(param))
                with Pool(processes=10) as pool:
                    results=pool.map(iso_pool,parameters)
                self.isochrone=results
    
            else:
            # self.isochrone_list=[[iso_low.metalicty,iso_low.age],[iso_bestfit.]]
                # define the normal distribution and PDF
                difference_wanted=[]
                dist = sps.norm(loc=0, scale=self.difference_low)
                percentile_pdf=np.linspace(0.5,1,num=numbers)
                for i in percentile_pdf:
                    difference_wanted.append(dist.ppf(i))
                self.points=difference_wanted
                parameters=[[self.iso_bestfit.parameters[0],self.iso_bestfit.parameters[1],0],[self.iso_low.parameters[0],self.iso_low.parameters[1],self.difference_low]]
                isochrones_list=[]
                iso_param_have=np.copy(parameters[0])
                isochrones_list.append(iso_bestfit.isochrone)
                difference_iso_have=[0,self.difference_low,np.inf]
                while len(difference_wanted)-2:
                    print('there are '+ str(len(difference_wanted)+numbers/2)+' isocrhones left')
                    print(str(len(isochrones_list))+' isocrhones are done')
                    digitized = np.digitize(difference_wanted, difference_iso_have)
                    bin_amounts = [len(digitized[digitized == i]) for i in range(1, len(difference_iso_have))]
                    bin_max_index=bin_amounts.index(max(bin_amounts))
                    point_max_index=int(len([x for x in digitized if x<bin_max_index]))
                    point_max_index+=int(np.round(bin_amounts[bin_max_index]/2))
                    point_max_wanted=difference_wanted[point_max_index]
                    
                    if bin_max_index==len(difference_iso_have)-2:
                        bin_low=min(5,len(parameters))
                        parameters_closest=np.array(parameters[-bin_low:])
                        
                        trying_loop=True
                        while trying_loop:
                            popt, pcov = curve_fit(linear_fit, parameters_closest[:,2], parameters_closest[:,0])
                            age_guess=linear_fit(point_max_wanted,popt[0] ,popt[1])
                            popt, pcov = curve_fit(linear_fit, parameters_closest[:,2], parameters_closest[:,1])
                            metalicty_guess=linear_fit(point_max_wanted,popt[0] ,popt[1])
        
                            if abs(age_guess-parameters_closest[1][0])>0.5 or abs(metalicty_guess-parameters_closest[1][1])>0.2:
                                diff=(point_max_wanted-parameters[-1][2])/2
                                point_max_wanted-=diff
                            elif [age_guess,metalicty_guess] in np.array(parameters)[:,:2]:
                                bin_max_index=np.where(np.array(parameters)[:,:2]==[age_guess,metalicty_guess])[0][0]
        
                                sign=np.sign(point_max_wanted-parameters[np.where(np.array(parameters)[:,:2]==[age_guess,metalicty_guess])[0][0]][2])
                                ratio=point_max_wanted/parameters[np.where(np.array(parameters)[:,:2]==[age_guess,metalicty_guess])[0][0]][2]
                                point_max_wanted*=ratio
                            else:
                                trying_loop=False
                    else:
                        loop=True
                        while loop:
                            bin_low=max(0,bin_max_index-2)
                            if bin_max_index+4<len(parameters):
                                bin_high=bin_max_index+4
                            else:
                                bin_high=len(parameters)
                            parameters_closest=np.array(parameters[bin_low:bin_high])
                            popt, pcov = curve_fit(linear_fit, parameters_closest[:,2], parameters_closest[:,0])
                            age_guess=linear_fit(point_max_wanted,popt[0] ,popt[1])
                            popt, pcov = curve_fit(linear_fit, parameters_closest[:,2], parameters_closest[:,1])
                            metalicty_guess=linear_fit(point_max_wanted,popt[0] ,popt[1])
                            if [age_guess,metalicty_guess] in np.array(parameters)[:,:2]:
                                bin_max_index=np.where(np.array(parameters)[:,:2]==[age_guess,metalicty_guess])[0][0]
                                ratio=point_max_wanted/parameters[np.where(np.array(parameters)[:,:2]==[age_guess,metalicty_guess])[0][0]][2]
                                point_max_wanted*=ratio
                            else:
                                loop=False
        
                        
                    if metalicty_guess>0.68:
                        metalicty_guess=0.68
                    elif metalicty_guess<-2.1:
                        metalicty_guess=-2.1
                    
                    if age_guess>10:
                        age_guess=10
                    if not ((metalicty_guess==0.68 or metalicty_guess==-2.1) and age_guess==10 and [age_guess,metalicty_guess] in np.array(parameters)[:,:2]):
                        iso_temp=isochrone('temp',age_guess,metalicty_guess,self.extinction,limits=limits,high_age=age_guess,high_metalicity=metalicty_guess)
                        difference_temp=difference_iso(iso_temp.isochrone,iso_bestfit.isochrone)
                        
                        difference_closest_index=np.where(abs(difference_wanted-difference_temp)==np.amin(abs(difference_wanted-difference_temp)))[0][0]
                        difference_closest=difference_wanted[difference_closest_index]
                        if abs(difference_closest-difference_temp)<self.difference_low/10:
                            del difference_wanted [difference_closest_index]
                            iso_param_have=np.vstack((iso_param_have,np.array([age_guess,metalicty_guess,difference_temp])))
                            isochrones_list.append(iso_temp.isochrone)
                            parameters.insert(bin_max_index+1,[age_guess,metalicty_guess,difference_temp])
                            difference_iso_have.insert(bin_max_index+1  ,difference_temp)
                        else:
                            parameters.insert(bin_max_index+1,[age_guess,metalicty_guess,difference_temp])
                            difference_iso_have.insert(bin_max_index+2,difference_temp)
                        parameters=sorted(parameters,key=lambda x:x[2])
                        difference_iso_have=sorted(difference_iso_have)
                    else:
                        del difference_wanted [point_max_index]
                    
                difference_wanted=[]
                dist = sps.norm(loc=0, scale=self.difference_high)
                percentile_pdf=np.linspace(0.5,1,num=numbers)
                for i in percentile_pdf:
                    difference_wanted.append(dist.ppf(i))
                self.points=difference_wanted
                
                parameters=[[self.iso_bestfit.parameters[0],self.iso_bestfit.parameters[1],0],[self.iso_high.parameters[0],self.iso_high.parameters[1],self.difference_high]]
        
                
                difference_iso_have=[0,self.difference_high,np.inf]
                while len(difference_wanted)-2:
                    print('there are '+ str(len(difference_wanted))+' isocrhones left')
                    print(str(len(isochrones_list))+' isocrhones are done')
                    digitized = np.digitize(difference_wanted, difference_iso_have)
                    bin_amounts = [len(digitized[digitized == i]) for i in range(1, len(difference_iso_have))]
                    bin_max_index=bin_amounts.index(max(bin_amounts))
                    point_max_index=int(len([x for x in digitized if x<bin_max_index]))
                    point_max_index+=int(np.round(bin_amounts[bin_max_index]/2))
                    point_max_wanted=difference_wanted[point_max_index]
                    
                    if bin_max_index==len(difference_iso_have)-2:
                        bin_low=min(5,len(parameters))
                        parameters_closest=np.array(parameters[-bin_low:])
                        
                        trying_loop=True
                        while trying_loop:
                            popt, pcov = curve_fit(linear_fit, parameters_closest[:,2], parameters_closest[:,0])
                            age_guess=linear_fit(point_max_wanted,popt[0] ,popt[1])
                            popt, pcov = curve_fit(linear_fit, parameters_closest[:,2], parameters_closest[:,1])
                            metalicty_guess=linear_fit(point_max_wanted,popt[0] ,popt[1])
        
                            if abs(age_guess-parameters_closest[1][0])>0.5 or abs(metalicty_guess-parameters_closest[1][1])>0.2:
                                diff=(point_max_wanted-parameters[-1][2])/2
                                point_max_wanted-=diff
                            elif [age_guess,metalicty_guess] in np.array(parameters)[:,:2]:
                                
                                bin_max_index=np.where(np.array(parameters)[:,:2]==[age_guess,metalicty_guess])[0][0]
        
                                sign=np.sign(point_max_wanted-parameters[np.where(np.array(parameters)[:,:2]==[age_guess,metalicty_guess])[0][0]][2])
                                ratio=point_max_wanted/parameters[np.where(np.array(parameters)[:,:2]==[age_guess,metalicty_guess])[0][0]][2]
                                point_max_wanted*=ratio
                            else:
                                trying_loop=False
                    else:
                        loop=True
                        while loop:
                            bin_low=max(0,bin_max_index-2)
                            if bin_max_index+4<len(parameters):
                                bin_high=bin_max_index+4
                            else:
                                bin_high=len(parameters)
                            parameters_closest=np.array(parameters[bin_low:bin_high])
                            popt, pcov = curve_fit(linear_fit, parameters_closest[:,2], parameters_closest[:,0])
                            age_guess=linear_fit(point_max_wanted,popt[0] ,popt[1])
                            popt, pcov = curve_fit(linear_fit, parameters_closest[:,2], parameters_closest[:,1])
                            metalicty_guess=linear_fit(point_max_wanted,popt[0] ,popt[1])
                            if [age_guess,metalicty_guess] in np.array(parameters)[:,:2]:
                                bin_max_index=np.where(np.array(parameters)[:,:2]==[age_guess,metalicty_guess])[0][0]
                                ratio=point_max_wanted/parameters[np.where(np.array(parameters)[:,:2]==[age_guess,metalicty_guess])[0][0]][2]
                                point_max_wanted*=ratio
                            else:
                                loop=False
                    if metalicty_guess>0.68:
                        metalicty_guess=0.68
                    elif metalicty_guess<-2.1:
                        metalicty_guess=-2.1
                    
                    if age_guess>10:
                        age_guess=10
                    elif age_guess<1.5:
                        age_guess=1.5
                    if not ((metalicty_guess==0.68 or metalicty_guess==-2.1) and (age_guess==10 or age_guess==1.5) and [age_guess,metalicty_guess] in np.array(parameters)[:,:2]):
                        iso_temp=isochrone('temp',age_guess,metalicty_guess,self.extinction,limits=limits,high_age=age_guess,high_metalicity=metalicty_guess)
                        difference_temp=difference_iso(iso_temp.isochrone,iso_bestfit.isochrone)
                        
                        difference_closest_index=np.where(abs(difference_wanted-difference_temp)==np.amin(abs(difference_wanted-difference_temp)))[0][0]
                        difference_closest=difference_wanted[difference_closest_index]
                        if abs(difference_closest-difference_temp)<self.difference_low/10:
                            del difference_wanted [difference_closest_index]
                            iso_param_have=np.vstack((iso_param_have,np.array([age_guess,metalicty_guess,difference_temp])))
                            isochrones_list.append(iso_temp.isochrone)
                            parameters.insert(bin_max_index+1,[age_guess,metalicty_guess,difference_temp])
                            difference_iso_have.insert(bin_max_index+1  ,difference_temp)
                        else:
                            parameters.insert(bin_max_index+1,[age_guess,metalicty_guess,difference_temp])
                            difference_iso_have.insert(bin_max_index+2,difference_temp)
                        parameters=sorted(parameters,key=lambda x:x[2])
                        difference_iso_have=sorted(difference_iso_have)
                    else:
                        del difference_wanted [point_max_index]
                self.isochrone=isochrones_list
                self.parameters=iso_param_have
        def plot(self):
            plt.figure()
            if isinstance(self.isochrone,np.ndarray):
                plt.plot(self.isochrone[:,5]-self.isochrone[:,4],self.isochrone[:,3])
            else:
                for (x,label) in zip(self.isochrone,self.parameters):
                    temp_label='age '+ str(label[0]) +' metalicty '+str(label[1])
                
                    plt.plot(x[:,4]-x[:,5],x[:,3],label=temp_label)
                if len(self.isochrone)<5:
                    plt.legend(loc='best')
            plt.xlabel('bp_rp')
            plt.ylabel('G')
            plt.gca().invert_yaxis()
            # for x in         
            # while x>
        # def splitting(self,)
def iso_pool(parameters_single):
    # print('asking for ')
    # print(parameters_single)
    iso_temp=isochrone('temp',parameters_single[0],parameters_single[1],parameters_single[2],high_age=parameters_single[0],high_metalicity=parameters_single[1])
    return iso_temp.isochrone

#interpolates linearly with mass (not log mass) as the basis.
def interpolateV(oldIsochrone):   
    iso_temp=copy.copy(oldIsochrone)
    mass=iso_temp[:,0]
    #changes log temperature to temp.
    iso_temp[:,2]=10**(iso_temp[:,2])
    #changes log g to g
#    oldIsochrone[:,7]=np.exp(oldIsochrone[:,7])
    scale=5
    
    #makes a new ISO with 5 times more slots 
    newIso=np.zeros(((len(iso_temp)-1)*scale,len(iso_temp[0])))
    
    
    newIso[:,0]=np.hstack(np.transpose(np.linspace(mass[:-1],mass[1:],num=5,endpoint=True)))
    
    
    newIso=np.array(newIso,dtype=float)

    for x in range(len(iso_temp[0])-1):
        f = interp1d(mass,iso_temp[:,x+1])
        newIso[:,x+1]=f(newIso[:,0]) 
    newIso[:,2]=np.log10(newIso[:,2])
#    newIso[:,7]=np.log(newIso[:,7])
    return np.array(newIso,dtype=float)
def parameter_error_calculator(stars,isochroneTemp,done):
    shift=stars['shift']
    point=np.array([stars['abs_phot_g_mean_mag'],stars['abs_phot_bp_mean_mag']-stars['abs_phot_rp_mean_mag']])
    sigma=[2.5/np.log(10)/stars['phot_g_mean_flux_over_error'],np.sqrt((2.5/np.log(10)/stars['phot_bp_mean_flux_over_error'])**2+(2.5/np.log(10)/stars['phot_rp_mean_flux_over_error'])**2)]
    sigmaD=5/np.log(10)*(stars['r_sigma_cluster'])/stars['r_est_cluster']
    sigma=np.array(sigma)
   
    sigma[0,:]=np.sqrt(sigma[0,:]**2+sigmaD**2)
    middlePoint=np.add(isochroneTemp[:-1,:],isochroneTemp[1:,:])/2
    numberDensity=np.subtract(isochroneTemp[1:,1],isochroneTemp[:-1,1])
    middlePointPhotometry=np.transpose([middlePoint[:,3],middlePoint[:,4]-middlePoint[:,5]])
    middleTemperature=np.add(10**isochroneTemp[1:,2],10**isochroneTemp[:-1,2])/2 
    isochroneTempPoint=np.array([isochroneTemp[:,3],isochroneTemp[:,4]-isochroneTemp[:,5]])
    distance=np.sqrt(np.sum((isochroneTempPoint[:,1:]-isochroneTempPoint[:,:-1])**2,axis=0))
    
#    middlelog_g=np.add(10**isochroneTemp[1:,6],10**isochroneTemp[:-1,6])/2
    middlelog_g=np.add(isochroneTemp[1:,6],isochroneTemp[:-1,6])/2

    middleMass=np.add(isochroneTemp[1:,-1],isochroneTemp[:-1,-1])/2
    
    
    e_mass=[]
    e_logg=[]
    e_teff=[]
    meanTemperature=stars['teff']
    meanLog_g=stars['logg']
    meanMass=stars['mass']
    for x in range(len(stars)):
        if np.all(done[x]==0):
            probabilityStar=-np.sum((point[:,x]-middlePointPhotometry)**2/(2*sigma[:,x]**2),axis=1)+shift[x]
            normalization=np.sum(numberDensity*np.exp(probabilityStar)*distance)
            upper=numberDensity*np.exp(probabilityStar)*distance
            e_teff.append(np.sqrt(np.sum((middleTemperature-meanTemperature[x])**2*upper)/normalization))
            e_logg.append(np.sqrt(np.sum((meanLog_g[x]-middlelog_g)**2*upper)/normalization))
            if e_logg[-1]<1e-5:
                e_logg[-1]=1e-5
            e_mass.append(np.sqrt(np.sum((meanMass[x]-middleMass)**2*upper)/normalization))
        else:
            e_teff.append(0)
            e_logg.append(0)
            e_mass.append(0)
    return[ e_logg,e_mass,e_teff]
def parameter_calculator(stars,isochroneTemp,done,parameter):
    shift=stars['shift']
    param_calculated=[]         
    point=np.array([stars['abs_phot_g_mean_mag'],stars['abs_phot_bp_mean_mag']-stars['abs_phot_rp_mean_mag']])
    sigma=[2.5/np.log(10)/stars['phot_g_mean_flux_over_error'],np.sqrt((2.5/np.log(10)/stars['phot_bp_mean_flux_over_error'])**2+(2.5/np.log(10)/stars['phot_rp_mean_flux_over_error'])**2)]
    sigmaD=5/np.log(10)*(stars['r_sigma_cluster'])/stars['r_est_cluster']
    sigma=np.array(sigma)
    sigma[0,:]=np.sqrt(sigma[0,:]**2+sigmaD**2)
    middlePoint=np.add(isochroneTemp[:-1,:],isochroneTemp[1:,:])/2
    numberDensity=np.subtract(isochroneTemp[1:,1],isochroneTemp[:-1,1])
    middlePointPhotometry=np.transpose([middlePoint[:,3],middlePoint[:,4]-middlePoint[:,5]])
    if parameter=='logg':
        middle_param=np.add(isochroneTemp[1:,6],isochroneTemp[:-1,6])/2
    if parameter=='teff':
        middle_param=np.add(10**isochroneTemp[1:,2],10**isochroneTemp[:-1,2])/2
    if parameter=='mass':
        middle_param=np.add(isochroneTemp[1:,-1],isochroneTemp[:-1,-1])/2

    isochroneTempPoint=np.array([isochroneTemp[:,3],isochroneTemp[:,4]-isochroneTemp[:,5]])

    
    distance=np.sqrt(np.sum((isochroneTempPoint[:,1:]-isochroneTempPoint[:,:-1])**2,axis=0))
    
    for x,y in zip(range(len(stars)),done):
        if y==0:
            temp_parameter=np.nan
            while np.isnan(temp_parameter):
                probabilityStar=-np.sum((point[:,x]-middlePointPhotometry)**2/(2*sigma[:,x]**2),axis=1)+shift[x]
                temp_parameter=np.sum(middle_param*numberDensity*np.exp(probabilityStar)*distance)/np.sum(numberDensity*np.exp(probabilityStar)*distance)
                if np.isnan(temp_parameter):
                    shift[x]=min(1+shift[x],700+min(np.sum((point[:,x]-middlePointPhotometry)**2/(2*sigma[:,x]**2),axis=1)))
                else:
                    param_calculated.append(temp_parameter)
        else:
            param_calculated.append(0)
    return param_calculated,shift


#Weights the log_g using Cluster distances 
def Wlog_gClusterV(stars,isochroneTemp,done):
    # crops the isochrone down to initial mass, int_imf, log_T, G, G_BP, G_RP ,log_g ,cur mass

    shift=stars['shift']
    logg=[]         
    point=np.array([stars['abs_phot_g_mean_mag'],stars['abs_phot_bp_mean_mag']-stars['abs_phot_rp_mean_mag']])
    sigma=[2.5/np.log(10)/stars['phot_g_mean_flux_over_error'],np.sqrt((2.5/np.log(10)/stars['phot_bp_mean_flux_over_error'])**2+(2.5/np.log(10)/stars['phot_rp_mean_flux_over_error'])**2)]
    sigmaD=5/np.log(10)*(stars['r_sigma_cluster'])/stars['r_est_cluster']
    sigma=np.array(sigma)
   
    sigma[0,:]=np.sqrt(sigma[0,:]**2+sigmaD**2)
    middlePoint=np.add(isochroneTemp[:-1,:],isochroneTemp[1:,:])/2
    numberDensity=np.subtract(isochroneTemp[1:,1],isochroneTemp[:-1,1])
    middlePointPhotometry=np.transpose([middlePoint[:,3],middlePoint[:,4]-middlePoint[:,5]])
    middlelogg=np.add(isochroneTemp[1:,6],isochroneTemp[:-1,6])/2
    isochroneTempPoint=np.array([isochroneTemp[:,3],isochroneTemp[:,4]-isochroneTemp[:,5]])

    
    distance=np.sqrt(np.sum((isochroneTempPoint[:,1:]-isochroneTempPoint[:,:-1])**2,axis=0))
    
    for x,y in zip(range(len(stars)),done):
        if y==0:
            temp_parameter=np.nan
            while np.isnan(temp_parameter):
                probabilityStar=-np.sum((point[:,x]-middlePointPhotometry)**2/(2*sigma[:,x]**2),axis=1)+shift[x]
                temp_parameter=np.sum(middlelogg*numberDensity*np.exp(probabilityStar)*distance)/np.sum(numberDensity*np.exp(probabilityStar)*distance)
                if np.isnan(temp_parameter):
                    shift[x]=min(20+shift[x],700+min(np.sum((point[:,x]-middlePointPhotometry)**2/(2*sigma[:,x]**2),axis=1)))
                else:
                    logg.append(temp_parameter)
        else:
            logg.append(0)
    return logg,shift

#Weights the log_g using Cluster distances 
def WmassClusterV(stars,isochroneTemp):
    shift=stars['shift']
    mass=[]         
    point=np.array([stars['abs_phot_g_mean_mag'],stars['abs_phot_bp_mean_mag']-stars['abs_phot_rp_mean_mag']])
    sigma=[2.5/np.log(10)/stars['phot_g_mean_flux_over_error'],np.sqrt((2.5/np.log(10)/stars['phot_bp_mean_flux_over_error'])**2+(2.5/np.log(10)/stars['phot_rp_mean_flux_over_error'])**2)]
    sigmaD=5/np.log(10)*(stars['r_sigma_cluster'])/stars['r_est_cluster']
    sigma=np.array(sigma)
   
    sigma[0,:]=np.sqrt(sigma[0,:]**2+sigmaD**2)
    middlePoint=np.add(isochroneTemp[:-1,:],isochroneTemp[1:,:])/2
    numberDensity=np.subtract(isochroneTemp[1:,1],isochroneTemp[:-1,1])
    middlePointPhotometry=np.transpose([middlePoint[:,3],middlePoint[:,4]-middlePoint[:,5]])
    middleMass=np.add(isochroneTemp[1:,-1],isochroneTemp[:-1,-1])/2
    isochroneTempPoint=np.array([isochroneTemp[:,3],isochroneTemp[:,4]-isochroneTemp[:,5]])

    
    distance=np.sqrt(np.sum((isochroneTempPoint[:,1:]-isochroneTempPoint[:,:-1])**2,axis=0))
    
    for x in range(len(stars)):
        temp_parameter=np.nan
        while np.isnan(temp_parameter):
            probabilityStar=-np.sum((point[:,x]-middlePointPhotometry)**2/(2*sigma[:,x]**2),axis=1)+shift[x]
            temp_parameter=np.sum(middleMass*numberDensity*np.exp(probabilityStar)*distance)/np.sum(numberDensity*np.exp(probabilityStar)*distance)
            if np.isnan(temp_parameter):
                shift[x]=min(20+shift[x],700+min(np.sum((point[:,x]-middlePointPhotometry)**2/(2*sigma[:,x]**2),axis=1)))
            else:
                mass.append(temp_parameter)
         
    return mass,shift


#Weights the temperature using Cluster distances 
def WtempVCluster(stars,isochroneTemp):
    shift=stars['shift']
    temperature=[]         
    point=np.array([stars['abs_phot_g_mean_mag'],stars['abs_phot_bp_mean_mag']-stars['abs_phot_rp_mean_mag']])
    sigma=[2.5/np.log(10)/stars['phot_g_mean_flux_over_error'],np.sqrt((2.5/np.log(10)/stars['phot_bp_mean_flux_over_error'])**2+(2.5/np.log(10)/stars['phot_rp_mean_flux_over_error'])**2)]
    sigmaD=5/np.log(10)*(stars['r_sigma_cluster'])/stars['r_est_cluster']
    sigma=np.array(sigma)
   
    sigma[0,:]=np.sqrt(sigma[0,:]**2+sigmaD**2)
    middlePoint=np.add(isochroneTemp[:-1,:],isochroneTemp[1:,:])/2
    numberDensity=np.subtract(isochroneTemp[1:,1],isochroneTemp[:-1,1])
    middlePointPhotometry=np.transpose([middlePoint[:,3],middlePoint[:,4]-middlePoint[:,5]])
    middleTemperature=np.add(10**isochroneTemp[1:,2],10**isochroneTemp[:-1,2])/2
    
    isochroneTempPoint=np.array([isochroneTemp[:,3],isochroneTemp[:,4]-isochroneTemp[:,5]])

    
    distance=np.sqrt(np.sum((isochroneTempPoint[:,1:]-isochroneTempPoint[:,:-1])**2,axis=0))
    
    for x in range(len(stars)):
        temp_parameter=np.nan
        while np.isnan(temp_parameter):
            probabilityStar=-np.sum((point[:,x]-middlePointPhotometry)**2/(2*sigma[:,x]**2),axis=1)+shift[x]
            temp_parameter=np.sum(middleTemperature*numberDensity*np.exp(probabilityStar)*distance)/np.sum(numberDensity*np.exp(probabilityStar)*distance)
            if np.isnan(temp_parameter):
                shift[x]=min(20+shift[x],700+min(np.sum((point[:,x]-middlePointPhotometry)**2/(2*sigma[:,x]**2),axis=1)))
            else:
                temperature.append(temp_parameter)
    return temperature,shift

#Weights the temperature using Cluster distances 
def UncertaintyV(stars,isochroneTemp):
    shift=stars['shift']
    point=np.array([stars['abs_phot_g_mean_mag'],stars['abs_phot_bp_mean_mag']-stars['abs_phot_rp_mean_mag']])
    sigma=[2.5/np.log(10)/stars['phot_g_mean_flux_over_error'],np.sqrt((2.5/np.log(10)/stars['phot_bp_mean_flux_over_error'])**2+(2.5/np.log(10)/stars['phot_rp_mean_flux_over_error'])**2)]
    sigmaD=5/np.log(10)*(stars['r_sigma_cluster'])/stars['r_est_cluster']
    sigma=np.array(sigma)
   
    sigma[0,:]=np.sqrt(sigma[0,:]**2+sigmaD**2)
    middlePoint=np.add(isochroneTemp[:-1,:],isochroneTemp[1:,:])/2
    numberDensity=np.subtract(isochroneTemp[1:,1],isochroneTemp[:-1,1])
    middlePointPhotometry=np.transpose([middlePoint[:,3],middlePoint[:,4]-middlePoint[:,5]])
    middleTemperature=np.add(10**isochroneTemp[1:,2],10**isochroneTemp[:-1,2])/2 
    isochroneTempPoint=np.array([isochroneTemp[:,3],isochroneTemp[:,4]-isochroneTemp[:,5]])
    distance=np.sqrt(np.sum((isochroneTempPoint[:,1:]-isochroneTempPoint[:,:-1])**2,axis=0))
    
#    middlelog_g=np.add(10**isochroneTemp[1:,6],10**isochroneTemp[:-1,6])/2
    middlelog_g=np.add(isochroneTemp[1:,6],isochroneTemp[:-1,6])/2

    middleMass=np.add(isochroneTemp[1:,-1],isochroneTemp[:-1,-1])/2
    
    
    sigmaMass=[]
    sigmaLog_g=[]
    sigmaTemperature=[]
    meanTemperature=stars['teff']
    meanLog_g=stars['logg']
    meanMass=stars['mass']
    for x in range(len(stars)):
        probabilityStar=-np.sum((point[:,x]-middlePointPhotometry)**2/(2*sigma[:,x]**2),axis=1)+shift[x]
        normalization=np.sum(numberDensity*np.exp(probabilityStar)*distance)
        upper=numberDensity*np.exp(probabilityStar)*distance
        sigmaTemperature.append(np.sqrt(np.sum((middleTemperature-meanTemperature[x])**2*upper)/normalization))
        sigmaLog_g.append(np.sqrt(np.sum((meanLog_g[x]-middlelog_g)**2*upper)/normalization))
        if sigmaLog_g[-1]<1e-5:
            sigmaLog_g[-1]=1e-5
        sigmaMass.append(np.sqrt(np.sum((meanMass[x]-middleMass)**2*upper)/normalization))     
    return sigmaTemperature,sigmaLog_g,sigmaMass
#see what is the power of the probability exponential
def powersClusterV(stars,isochroneTemp):
    power=[]
    point=np.array([stars['abs_phot_g_mean_mag'],stars['abs_phot_bp_mean_mag']-stars['abs_phot_rp_mean_mag']])
    sigma=[2.5/np.log(10)/stars['phot_g_mean_flux_over_error'],np.sqrt((2.5/np.log(10)/stars['phot_bp_mean_flux_over_error'])**2+(2.5/np.log(10)/stars['phot_rp_mean_flux_over_error'])**2)]
    sigmaD=5/np.log(10)*(stars['r_sigma_cluster'])/stars['r_est_cluster']
    sigma=np.array(sigma)
   
    sigma[0,:]=np.sqrt(sigma[0,:]**2+sigmaD**2)
    middlePointIso=np.add(isochroneTemp[:-1,:],isochroneTemp[1:,:])/2
    middlePointPhotometry=[middlePointIso[:,3],middlePointIso[:,4]-middlePointIso[:,5]]  
    for x in range(len(stars)):
        power.append(250+min(np.sum((point[:,x]-np.transpose(middlePointPhotometry))**2/(2*sigma[:,x]**2),axis=1)))
    return power

def powersClusterIndividual(stars,isochroneTemp):
    power=[]
    point=[stars['abs_phot_g_mean_mag'],stars['abs_phot_bp_mean_mag'],stars['abs_phot_rp_mean_mag']]
    point=np.array(point)
    sigma=[2.5/np.log(10)/stars['phot_g_mean_flux_over_error'],2.5/np.log(10)/stars['phot_bp_mean_flux_over_error'],2.5/np.log(10)/stars['phot_rp_mean_flux_over_error']]
    sigmaD=5/np.log(10)*(stars['r_sigma_cluster'])/stars['r_est_cluster']
    sigma=np.array(sigma)
    sigma=np.sqrt(sigma**2+sigmaD**2)  
    middlePointIso=np.add(isochroneTemp[:-1,:],isochroneTemp[1:,:])/2
    middlePointPhotometry=[middlePointIso[:,3],middlePointIso[:,4],middlePointIso[:,5]]  
    power=(200+min(np.sum((point-np.transpose(middlePointPhotometry))**2/(2*sigma**2),axis=1)))
    return power

def isoCropper(isochrone):
    #crops the isochrone down to initial ini mass, int_imf, log_T, G, G_Bp_bright, G_BP_faint, G_RP ,log_g , cur mass
    croppedIsochone=np.vstack((isochrone[:,3],isochrone[:,4],isochrone[:,7],isochrone[:,25],isochrone[:,26],isochrone[:,26],isochrone[:,27],isochrone[:,8],isochrone[:,5]))
    croppedIsochone=np.transpose(croppedIsochone)
    cIsoTemp=croppedIsochone[0,:]    
    dist=np.sqrt(np.sum((croppedIsochone[1:,[3,4,6]]-croppedIsochone[:-1,[3,4,6]])**2,axis=1))
    while np.any(dist>8.0):
        cIsoTemp=croppedIsochone[0,:]
        for i in range(len(dist)):
            if dist[i]<2.0:
                cIsoTemp=np.vstack((cIsoTemp,croppedIsochone[i+1,:]))
        croppedIsochone=cIsoTemp    
        dist=np.sqrt(np.sum((croppedIsochone[1:,[3,4,6]]-croppedIsochone[:-1,[3,4,6]])**2,axis=1))
    return croppedIsochone

def intergration_loop(data,cropped_iso,parameter,limit):
    #finished parameters goes into here
    parameters_loop=np.array(np.zeros(len(data)),dtype=object)
    
    #Does the first loop 
    parameters_old=np.zeros(len(data))    
    if parameter!='error':
        parameters_calculated,data['shift']=parameter_calculator(data,cropped_iso,parameters_loop,parameter)
    else:
        parameters_calculated=parameter_error_calculator(data, cropped_iso, parameters_loop)
    while np.any([np.all(x==0) for x in parameters_loop]):
    	#Interpolates spectra
        cropped_iso=interpolateV(cropped_iso)
        parameters_old=parameters_calculated
        if parameter!='error':
            parameters_calculated,data['shift']=parameter_calculator(data,cropped_iso,parameters_loop,parameter)
        else:
            parameters_calculated=parameter_error_calculator(data, cropped_iso, parameters_loop)
        if parameter!='error':
            if np.any((abs(np.array(parameters_old)-np.array(parameters_calculated))<limit) *( parameters_loop==0)) or len(cropped_iso)>1e7:
                pos_finished=np.where((abs(np.array(parameters_old)-np.array(parameters_calculated))<limit) *( parameters_loop==0))
                for x in pos_finished[0]:
                    parameters_loop[x]=parameters_calculated[x] 
                if len(cropped_iso)>1e7:
                    pos_finished=np.where([np.all(x==0) for x in parameters_loop])
                    for x in pos_finished[0]:
                        parameters_loop[x]=parameters_calculated[x] 

        else:
            if np.any(np.all((abs(np.array(parameters_old)-np.array(parameters_calculated)).T<limit),axis=1) *( [np.all(x==0) for x in parameters_loop]))or len(cropped_iso)>1e7:
                pos_finished=np.where(np.all((abs(np.array(parameters_old)-np.array(parameters_calculated)).T<limit),axis=1) *( [np.all(x==0) for x in parameters_loop]))
                for x in pos_finished[0]:
                    parameters_loop[x]=np.array(parameters_calculated).T[x] 
                if len(cropped_iso)>1e7:
                    pos_finished=np.where([np.all(x==0) for x in parameters_loop])
                    for x in pos_finished[0]:
                        parameters_loop[x]=np.array(parameters_calculated).T[x] 

    if parameter=='teff':
        data['teff']=parameters_loop
    elif parameter=='logg':
        data['logg']=parameters_loop
    elif parameter=='mass':
        data['mass']=parameters_loop
    elif parameter=='error':
        parameters_loop=np.vstack(parameters_loop)
        data['e_logg']=parameters_loop[:,0]
        data['e_mass']=parameters_loop[:,1]
        data['e_teff']=parameters_loop[:,2]

    return data,cropped_iso

    
#The main loop to interpolate iscorhones until the temperature doesnt change
def mainPart(dataLoop,isochrone):
    temperature_loop=[]
    # Zini MH logAge Mini int_IMF Mass logL logTe logg label McoreTP C_O period0 period1 pmode Mloss tau1m X Y Xc Xn Xo Cexcess Z 	 mbolmag Gmag G_BPmag G_RPmag
    # crops the isochrone down to initial mass, int_imf, log_T, G, G_BP, G_RP ,log_g ,cur mass
    cropped_iso=np.copy(isochrone)
    cIsoTemp=cropped_iso[0,:]
    
    dist=np.sqrt(np.sum((cropped_iso[1:,[3,4,6]]-cropped_iso[:-1,[3,4,6]])**2,axis=1))
    while np.any(dist>8.0):
        cIsoTemp=cropped_iso[0,:]
        for i in range(len(dist)):
            if dist[i]<2.0:
                cIsoTemp=np.vstack((cIsoTemp,cropped_iso[i+1,:]))
        cropped_iso=cIsoTemp    
        dist=np.sqrt(np.sum((cropped_iso[1:,[3,4,6]]-cropped_iso[:-1,[3,4,6]])**2,axis=1))
    
    original_isochrone=copy.deepcopy(cropped_iso)
    
    #splits the stars into dim and bright stars
    stars_loop=copy.copy(dataLoop)

    
    
    #calculates shift which is by how much we must shift the exponent so it will show up when integrating
    shift=powersClusterV(dataLoop,cropped_iso)
    stars_loop['shift']=shift
    

    
    print("calculating log_g")
    cropped_iso=copy.deepcopy(original_isochrone)
    stars_loop,cropped_iso=intergration_loop(stars_loop,cropped_iso,'logg',0.002)
    
    
    print("caluclating mass")
    cropped_iso=copy.deepcopy(original_isochrone)
    stars_loop,cropped_iso=intergration_loop(stars_loop,cropped_iso,'mass',0.002)
    print("caluclating teff")
    cropped_iso=copy.deepcopy(original_isochrone)
    stars_loop,cropped_iso=intergration_loop(stars_loop,cropped_iso,'teff',0.5)
    
    print("calculating uncertainty")
    cropped_iso=copy.deepcopy(original_isochrone)
    #logg mass teff
    stars_loop,cropped_iso=intergration_loop(stars_loop,cropped_iso,'error',[0.001,0.001,0.25])
    print("finsihed everything for the loop")
    # sigmaTemperatureBright,sigmaLog_gBright,sigmaMassBright=UncertaintyV(stars_loop,cropped_iso)
    
    # stars_loop['e_teff']=sigmaTemperatureBright


    # stars_loop['e_logg']=sigmaLog_gBright
    
    # stars_loop['e_mass']=sigmaMassBright


    stars_loop.sort('source_id')

    return stars_loop

def temperature_getter(stars_xml,isochrome):
    stars_temp=copy.deepcopy(stars_xml)
    data_temp=mainPart(stars_temp,isochrome)
    data_temp.sort('source_id')
    return data_temp['teff'],data_temp['e_teff'],data_temp['logg'],data_temp['e_logg'],data_temp['mass'],data_temp['e_mass']
def plot_hr(isochrones,data,save=False,paper=True,name=False,title=False):
    if not name:
        name=isochrones.name
    name=re.sub('_',' ', name)
    plt.rc('font', size=15)
    iso_bestfit=isochrones.iso_bestfit.isochrone
    iso_low=isochrones.iso_low.isochrone
    iso_high=isochrones.iso_high.isochrone
    fig=plt.figure(figsize=(5.0,6.0))
    ax=fig.gca()
    ax.tick_params(direction='in')
    if isochrones.iso_type=='gaia':
        
        plt.plot(iso_bestfit[:,4]-iso_bestfit[:,5],iso_bestfit[:,3],label=r'Best fit',c='red',linewidth=1.5,alpha=0.7)
        plt.plot(iso_low[:,4]-iso_low[:,5],iso_low[:,3],label=r'$\pm1 \sigma$',c='green',zorder=2,linewidth=1.5,alpha=0.7)
        plt.plot(iso_high[:,4]-iso_high[:,5],iso_high[:,3],c='green',zorder=2,linewidth=1.5,alpha=0.7)
        sigma=[2.5/np.log(10)/data['phot_g_mean_flux_over_error'],2.5/np.log(10)/data['phot_bp_mean_flux_over_error'],2.5/np.log(10)/data['phot_rp_mean_flux_over_error']]
        sigmaDCluster=5/np.log(10)*(data['r_sigma_cluster'])/data['r_est_cluster']
        # 
        errorYCluster=np.sqrt(sigmaDCluster**2+sigma[0]**2)
        errorX=np.sqrt(sigma[2]**2+sigma[1]**2)

        limitsX=[0.53,2.5]
        limitsY=[-1,10]
        plt.scatter(data['bp_rp'],data['abs_phot_g_mean_mag'],s=0.3,zorder=0,c='black' )
        # plt.errorbar(data['bp_rp'],data['abs_phot_g_mean_mag'],xerr=errorX, yerr=errorYCluster,fmt='o',ms=0.5,elinewidth=0.4,zorder=3,c='blue' )
        plt.xlabel(r'$G_{\rm{BP}}-G_{\rm{RP}}$',)
        plt.ylabel(r'$M_{\rm{G}}$')
        
        plt.xlim(limitsX)
        plt.ylim(limitsY)
        name_save=name
        if title:
            plt.title(title)
        if not paper:
            isochrone_fill=isochrones.isochrone
            for x in isochrone_fill:
                plt.plot(x[:,4]-x[:,5],x[:,3],alpha=0.1,c='blue',zorder=0)
    elif isochrones.iso_type=='allwise':
        iso_bestfit_gaia=isochrones.iso_bestfit.iso_gaia
        iso_low_gaia=isochrones.iso_low.iso_gaia
        iso_high_gaia=isochrones.iso_high.iso_gaia

        plt.plot(iso_bestfit_gaia[:,4]-iso_bestfit[:,-4],iso_bestfit_gaia[:,3],label=r'Best fit',c='red',linewidth=0.8)
        plt.plot(iso_low_gaia[:,4]-iso_low[:,-4],iso_low_gaia[:,3],label=r'$\pm1 \sigma$',c='green',zorder=2,linewidth=0.8)
        plt.plot(iso_high_gaia[:,4]-iso_high[:,-4],iso_high_gaia[:,3],c='green',zorder=2,linewidth=0.8)
    
        sigma=[2.5/np.log(10)/data['phot_g_mean_flux_over_error'],2.5/np.log(10)/data['phot_bp_mean_flux_over_error'],data['w1mpro_error']]
        sigmaDCluster=5/np.log(10)*(data['r_sigma_cluster'])/data['r_est_cluster']
        # 
        errorYCluster=np.sqrt(sigmaDCluster**2+sigma[0]**2)
        errorX=np.sqrt(sigma[2]**2+sigma[1]**2)

        limitsX=[np.nanmin(data['phot_bp_mean_mag']-data['w1mpro'])-0.1,max(data['phot_bp_mean_mag']-data['w1mpro'])+0.1]
        limitsY=[min(data['abs_phot_g_mean_mag'])-1,max(data['abs_phot_g_mean_mag'])+1]    
    
        plt.errorbar(data['phot_bp_mean_mag']-data['w1mpro'],data['abs_phot_g_mean_mag'],xerr=errorX, yerr=errorYCluster,fmt='o',ms=0.5,elinewidth=0.4,zorder=3,c='blue' )
        plt.xlabel(r'$G_{\rm{BP}}-W1$')
        plt.ylabel(r'$M_{\rm{G}}$')
        plt.xlim(limitsX)
        plt.ylim(limitsY)
        
        name_save=name+'_allwise'
    
        if not paper:
            isochrone_fill=isochrones.isochrone
            for x in isochrone_fill:
                plt.plot(x[:,4]-x[:,5],x[:,3],alpha=0.1,c='blue',zorder=0)
    elif isochrones.iso_type=='panstarrs':

        plt.plot(iso_bestfit[:,-6]-iso_bestfit[:,-2],iso_bestfit[:,-4],label=r'Best fit',c='red',linewidth=0.8)
        plt.plot(iso_low[:,-6]-iso_low[:,-2],iso_low[:,-4],label=r'$\pm1 \sigma$',c='green',zorder=2,linewidth=0.8)
        plt.plot(iso_high[:,-6]-iso_high[:,-2],iso_high[:,-4],c='green',zorder=2,linewidth=0.8)
    
        sigma=[data['i_mean_psf_mag_error'],data['g_mean_psf_mag_error'],data['y_mean_psf_mag_error']]
        sigmaDCluster=5/np.log(10)*(data['r_sigma_cluster'])/data['r_est_cluster']
        # 
        errorYCluster=np.sqrt(sigmaDCluster**2+sigma[0]**2)
        errorX=np.sqrt(sigma[2]**2+sigma[1]**2)
        
        Absmag=data['i_mean_psf_mag']-5*(np.log10(data['r_est_cluster'])-1)
        data['abs_i_mean_psf_mag']=Absmag

        limitsX=[np.nanmin(data['g_mean_psf_mag']-data['y_mean_psf_mag'])-0.1,max(data['g_mean_psf_mag']-data['y_mean_psf_mag'])+0.1]
        limitsY=[min(data['abs_i_mean_psf_mag'])-1,max(data['abs_i_mean_psf_mag'])+1]    
    
        plt.errorbar(data['g_mean_psf_mag']-data['y_mean_psf_mag'],data['abs_i_mean_psf_mag'],xerr=errorX, yerr=errorYCluster,fmt='o',ms=0.5,elinewidth=0.4,zorder=3,c='blue' )
        plt.xlabel(r'$G-Y1$')
        plt.ylabel(r'$M_{\rm{I}}$')
        plt.xlim(limitsX)
        plt.ylim(limitsY)
        
        name_save=name+'_allwise'
    
        if not paper:
            isochrone_fill=isochrones.isochrone
            for x in isochrone_fill:
                plt.plot(x[:,4]-x[:,5],x[:,3],alpha=0.1,c='blue',zorder=0)

    plt.legend(loc='best')
    if not name:
        name=isochrones.name
    plt.gca().invert_yaxis()
    plt.tight_layout()
    if save:
        plt.savefig(name_save+'.pdf',dpi=100,format='pdf')
def plot_keil(isochrones,data,save=False,paper=True,name=False):
    if not name:
        name=isochrones.name
    name=re.sub('_',' ', name)

    iso_bestfit=isochrones.iso_bestfit.isochrone
    iso_low=isochrones.iso_low.isochrone
    iso_high=isochrones.iso_high.isochrone
    fig=plt.figure(figsize=(4.0,5.0))
    ax=fig.gca()
    ax.tick_params(direction='in')
    plt.plot(10**iso_bestfit[:,2],iso_bestfit[:,-2],label=r'Best fit',c='red',linewidth=0.8)
    plt.plot(10**iso_low[:,2],iso_low[:,-2],label=r'$\pm1 \sigma$',c='green',zorder=2,linewidth=0.8)
    plt.plot(10**iso_high[:,2],iso_high[:,-2],c='green',zorder=2,linewidth=0.8)

    if not paper:
        isochrone_fill=isochrones.isochrone
        for x in isochrone_fill:
            plt.plot(10**x[:,2],x[:,-2],alpha=0.1,c='blue',zorder=0)

    plt.errorbar(data['teff'],data['logg'],xerr=data['e_teff'], yerr=data['e_logg'],fmt='o',ms=0.5,elinewidth=0.4,zorder=3,c='blue' )
    plt.xlabel(r'$G_{\rm{BP}}-G_{\rm{RP}}$')
    plt.ylabel(r'$M_{\rm{G}}$')


    plt.legend(loc='best')
    if not name:
        name=isochrones.name
    plt.gca().invert_yaxis()
    plt.tight_layout()
    if save:
        plt.savefig(name+'.pdf',dpi=100,format='pdf')
def plot_hr_teff(isochrones,data):
    iso_temp=isochrones.isochrone
    # iso_low=isochrones.iso_low.isochrone
    # iso_high=isochrones.iso_high.isochrone
    fig=plt.figure()
    ax=fig.gca()
    numberDensity=np.subtract(iso_temp[1:,1],iso_temp[:-1,1])
    np.append(numberDensity,numberDensity[-1])
    numberDensity*=1e3
    points = np.array([iso_temp[:,4]-iso_temp[:,5],iso_temp[:,3]]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    lc = LineCollection(segments, cmap='viridis', array=10**iso_temp[:,2],zorder=3)
    # Set the values used for colormapping
    lc.set_array(10**iso_temp[:-2,2])
     
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    # lc = LineCollection(segments, cmap='viridis', norm=norm)
    line = ax.add_collection(lc)
    plt.colorbar(line, label='teff')
    # plt.scatter(data['bp_rp'],data['abs_phot_g_mean_mag'],color='b',s=2)
    plt.scatter(data['bp_rp'],data['abs_phot_g_mean_mag'],c=data['teff'])
    plt.clim(min(10**iso_temp[:-2,2]),max(10**iso_temp[:-2,2]))
    plt.gca().invert_yaxis()
    plt.tight_layout()
def plot_hr_teff_2(isochrones,data):
    iso_temp=isochrones.isochrone
    # iso_low=isochrones.iso_low.isochrone
    # iso_high=isochrones.iso_high.isochrone
    fig=plt.figure()
    ax=fig.gca()
    plt.scatter(data['abs_phot_bp_mean_mag'],data['abs_phot_g_mean_mag'],c=data['teff'])
    numberDensity=np.subtract(iso_temp[1:,1],iso_temp[:-1,1])
    np.append(numberDensity,numberDensity[-1])
    numberDensity*=1e2
    points = np.array([iso_temp[:,4],iso_temp[:,3]]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    lc = LineCollection(segments, cmap='viridis', array=10**iso_temp[:,2],zorder=3)
    # Set the values used for colormapping
    lc.set_array(10**iso_temp[:-2,2])
     
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    # lc = LineCollection(segments, cmap='viridis', norm=norm)
    line = ax.add_collection(lc)
    plt.colorbar(line, label='teff')
    # plt.scatter(data['bp_rp'],data['abs_phot_g_mean_mag'],color='b',s=2)
    plt.gca().invert_yaxis()
    plt.tight_layout()
def plot_hr_teff_3(isochrones,data):
    iso_temp=isochrones.isochrone
    # iso_low=isochrones.iso_low.isochrone
    # iso_high=isochrones.iso_high.isochrone
    fig=plt.figure()
    ax=fig.gca()
    plt.scatter(data['abs_phot_rp_mean_mag'],data['abs_phot_g_mean_mag'],c=data['teff'])
    numberDensity=np.subtract(iso_temp[1:,1],iso_temp[:-1,1])
    np.append(numberDensity,numberDensity[-1])
    numberDensity*=1e4
    points = np.array([iso_temp[:,5],iso_temp[:,3]]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    lc = LineCollection(segments, cmap='viridis', array=10**iso_temp[:,2],zorder=3,linewidths=numberDensity)
    # Set the values used for colormapping
    lc.set_array(10**iso_temp[:-2,2])
     
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    # lc = LineCollection(segments, cmap='viridis', norm=norm)
    line = ax.add_collection(lc)
    plt.colorbar(line, label='teff')
    # plt.scatter(data['bp_rp'],data['abs_phot_g_mean_mag'],color='b',s=2)
    plt.gca().invert_yaxis()
    plt.tight_layout()
def number_density_lines(isochrones,data):
    iso_bestfit=isochrones.iso_bestfit.isochrone
    iso_low=isochrones.iso_low.isochrone
    iso_high=isochrones.iso_high.isochrone
    fig=plt.figure()
    ax=fig.gca()
    plt.scatter(data['bp_rp'],data['abs_phot_g_mean_mag'],c=data['teff'])
    numberDensity=np.subtract(iso_bestfit[1:,1],iso_bestfit[:-1,1])
    np.append(numberDensity,numberDensity[-1])
    plt.plot(data['bp_rp'],data['abs_phot_g_mean_mag'],linewidths=numberDensity)
    # plt.scatter(data['bp_rp'],data['abs_phot_g_mean_mag'],color='b',s=2)
    plt.gca().invert_yaxis()
    plt.tight_layout()
def curating_results(results,data):
    print('calculating things')
    teff_mean=[]
    teff_sigma=[]
    logg_mean=[]
    logg_sigma=[]
    teff_array=[]
    logg_array=[]
    probability_grid=[]
    teff_raw=[]
    logg_raw=[]
    e_logg_raw=[]
    e_teff_raw=[]
    mass_mean=[]
    mass_sigma=[]
    mass_raw=[]
    e_mass_raw=[]
    for star_index in range(len(data)): 
        temperature_star=[x[0][star_index] for x in results ]
        logg_star=[x[2][star_index] for x in results ]
        error_temperature=[x[1][star_index] for x in results ]
        error_logg=[x[3][star_index] for x in results ]
        mass_star=[x[4][star_index] for x in results ]
        error_mass=[x[5][star_index] for x in results ]
        
        mass_lin=np.linspace(min(mass_star)-2*max(error_mass),max(mass_star)+2*max(error_mass),num=100)
        teff_lin=np.linspace(min(temperature_star)-2*max(error_temperature),max(temperature_star)+2*max(error_temperature),num=100)
        teff_array.append(teff_lin)
        logg_lin=np.linspace(min(logg_star)-2*max(error_logg),max(logg_star)+2*max(error_logg),num=100)
        logg_array.append(logg_lin)
        xx,yy=np.meshgrid(teff_lin,logg_lin)
        zz=np.zeros((len(teff_lin),len(logg_lin)))
        for teff,e_teff,logg,e_logg in zip(temperature_star,error_temperature,logg_star,error_logg):
            for i in range(len(teff_lin)):
                for j in range(len(logg_lin)):
                    zz[i,j]+=normal(xx[j,j],teff,e_teff)*normal(yy[i,j],logg,e_logg)
                    
        zz/=len(temperature_star)
        probability_grid.append(zz)
        teff_zero=np.zeros(len(teff_lin))
        for teff,e_teff in zip(temperature_star,error_temperature):
            teff_zero+=[normal(x,teff,e_teff) for x in teff_lin]
        teff_mean.append(teff_lin[np.where(teff_zero==np.max(teff_zero))][0])
        teff_sigma.append(np.sqrt(sum([(teff_mean[star_index]-x)**2*y for (x,y) in zip(teff_lin,teff_zero)])/sum(teff_zero)))
        
        logg_zero=np.zeros(len(logg_lin))
        for logg,e_logg in zip(logg_star,error_logg):
            logg_zero+=[normal(x,logg,e_logg) for x in logg_lin]
        logg_mean.append(logg_lin[np.where(logg_zero==np.max(logg_zero))][0])
        logg_sigma.append(np.sqrt(sum([(logg_mean[star_index]-x)**2*y for (x,y) in zip(logg_lin,logg_zero)])/sum(logg_zero)))
        
        mass_zero=np.zeros(len(mass_lin))
        for mass,e_mass in zip(mass_star,error_mass):
            mass_zero+=[normal(x,mass,e_mass) for x in mass_lin]
        mass_mean.append(mass_lin[np.where(mass_zero==np.max(mass_zero))][0])
        mass_sigma.append(np.sqrt(sum([(mass_mean[star_index]-x)**2*y for (x,y) in zip(mass_lin,mass_zero)])/sum(mass_zero)))
    
        teff_raw.append(temperature_star)
        e_teff_raw.append(error_temperature)
        logg_raw.append(logg_star)
        e_logg_raw.append(error_logg)
        mass_raw.append(mass_star)
        e_mass_raw.append(error_mass)
    data['teff']=teff_mean
    data['e_teff']=teff_sigma
    data['logg']=logg_mean
    data['e_logg']=logg_sigma
    data['teff_array']=teff_array
    data['logg_array']=logg_array
    data['probability_grid']=probability_grid
    data['teff_raw']=teff_raw
    data['e_teff_raw']=e_teff_raw
    data['logg_raw']=logg_raw
    data['e_logg_raw']=e_logg_raw
    data['mass']=mass_mean
    data['e_mass']=mass_sigma
    data['mass_raw']=mass_raw
    data['e_mass_raw']=e_mass_raw
    return data
def curating_loop(loop):
        temperature_star=loop[0]
        logg_star=loop[2]
        error_temperature=loop[1]
        error_logg=loop[3] 
        mass_star=loop[4] 
        error_mass=loop[5] 
        
        mass_lin=np.linspace(min(mass_star)-2*max(error_mass),max(mass_star)+2*max(error_mass),num=100)
        teff_lin=np.linspace(min(temperature_star)-2*max(error_temperature),max(temperature_star)+2*max(error_temperature),num=100)
        logg_lin=np.linspace(min(logg_star)-2*max(error_logg),max(logg_star)+2*max(error_logg),num=100)
        xx,yy=np.meshgrid(teff_lin,logg_lin)
        zz=np.zeros((len(teff_lin),len(logg_lin)))
        for teff,e_teff,logg,e_logg in zip(temperature_star,error_temperature,logg_star,error_logg):
            for i in range(len(teff_lin)):
                for j in range(len(logg_lin)):
                    zz[i,j]+=normal(xx[j,j],teff,e_teff)*normal(yy[i,j],logg,e_logg)
                    
        zz/=len(temperature_star)
        teff_zero=np.zeros(len(teff_lin))
        for teff,e_teff in zip(temperature_star,error_temperature):
            teff_zero+=[normal(x,teff,e_teff) for x in teff_lin]
        
        logg_zero=np.zeros(len(logg_lin))
        for logg,e_logg in zip(logg_star,error_logg):
            logg_zero+=[normal(x,logg,e_logg) for x in logg_lin]
        
        mass_zero=np.zeros(len(mass_lin))
        for mass,e_mass in zip(mass_star,error_mass):
            mass_zero+=[normal(x,mass,e_mass) for x in mass_lin]
        mass_mean=mass_lin[np.where(mass_zero==np.max(mass_zero))][0]
        mass_sigma=np.sqrt(sum([(mass_mean-x)**2*y for (x,y) in zip(mass_lin,mass_zero)])/sum(mass_zero))
        teff_mean=teff_lin[np.where(teff_zero==np.max(teff_zero))][0]
        teff_sigma=np.sqrt(sum([(teff_mean-x)**2*y for (x,y) in zip(teff_lin,teff_zero)])/sum(teff_zero))
        logg_mean=logg_lin[np.where(logg_zero==np.max(logg_zero))][0]
        logg_sigma=np.sqrt(sum([(logg_mean-x)**2*y for (x,y) in zip(logg_lin,logg_zero)])/sum(logg_zero))
        probability_grid=zz
        teff_raw=temperature_star
        e_teff_raw=error_temperature
        logg_raw=logg_star
        e_logg_raw=error_logg
        mass_raw=mass_star
        e_mass_raw=error_mass
        teff_array=teff_lin
        logg_array=logg_lin
        # pbar.update(1)
        return teff_mean,teff_sigma,logg_mean,logg_sigma,mass_mean,mass_sigma,teff_raw,e_teff_raw,logg_raw,\
                e_logg_raw,mass_raw,e_mass_raw,teff_array,logg_array,probability_grid
# isochrone('temp',7.790753776549865,-0.02616414896675308,0.05,high_age=7.790753776549865,high_metalicity=-0.02616414896675308)
cluster_details_all = list(csv.reader(open('dr3_clusters.txt', 'rt'), delimiter=','))
cluster_details_all = np.array(cluster_details_all)
print(cluster_details_all[:,0])
# name=input('what isochrone do you want?')
name='NGC_2682'
cluster_details=[x for x in cluster_details_all if x[0]==name][0]
iso_type='gaia'

low_age=float(cluster_details[1])
low_metalicty=float(cluster_details[2])
low_extinction=float(cluster_details[3])
best_age=float(cluster_details[4])
best_metalicty=float(cluster_details[5])
best_extinction=float(cluster_details[6])
high_age=float(cluster_details[7])
high_metalicty=float(cluster_details[8])
high_extinction=float(cluster_details[9])


np.random.seed(0)
random.seed(0)
temperatureLimit=0.1
interpolate_scale=20

#opens the data for the specific star you have chosen
if iso_type=='gaia':
    votable = parse(name+"_dr3_distances3.xml")
elif iso_type=='allwise':
    votable = parse(name+"_allwise.xml")
elif iso_type=='panstarrs':
    votable = parse(name+"_panstarrs.xml")
data=votable.get_first_table().to_table(use_names_over_ids=True)

data=vstack(data)
data.sort('source_id')


print(name)
print(len(data))
# data.sort('phot_g_mean_mag')
# data_turning=[x for x in data if x['bp_rp']>0.95 and x['bp_rp']<1.2 and x['phot_g_mean_mag']>2.8+9.6 and x['phot_g_mean_mag']<3.3+9.6]
# data=data[:100]
# data=vstack([data,data_turning[2]])
# data=data[1:3]
    
#Changes the magnitudes of the stars to absolute magnitudes
# data=vstack([x for x in data if abs(x['r_est_cluster']-np.mean(data['r_est_cluster']))/x['r_sigma_cluster']<3])
print(len(data))
Absmag=data['phot_g_mean_mag']-5*(np.log10(data['r_est_cluster'])-1)
data['abs_phot_g_mean_mag']=Absmag
Absmag=data['phot_rp_mean_mag']-5*(np.log10(data['r_est_cluster'])-1)
data['abs_phot_rp_mean_mag']=Absmag
Absmag=data['phot_bp_mean_mag']-5*(np.log10(data['r_est_cluster'])-1)
data['abs_phot_bp_mean_mag']=Absmag

sigma=[2.5/np.log(10)/data['phot_g_mean_flux_over_error'],2.5/np.log(10)/data['phot_bp_mean_flux_over_error'],2.5/np.log(10)/data['phot_rp_mean_flux_over_error']]
sigmaDCluster=5/np.log(10)*(data['r_sigma_cluster'])/data['r_est_cluster']
# 
errorYCluster=np.sqrt(sigmaDCluster**2+sigma[0]**2)
errorX=np.sqrt(sigma[2]**2+sigma[1]**2)

limitsX=[min(data['bp_rp'])-0.1,max(data['bp_rp'])+0.1]
limitsY=[min(data['abs_phot_g_mean_mag'])-1,max(data['abs_phot_g_mean_mag'])+1]
# plot_hr(iso_blanco, data) 
limit_g=None
iso_bestfit=isochrone(name,best_age,best_metalicty,best_extinction,special='Best_fit',limits=limit_g,high_age=best_age,high_metalicity=best_metalicty,iso_type=iso_type,interpolate_scale=interpolate_scale)
iso_low_age=isochrone(name, low_age, best_metalicty, best_extinction,special='low_fit_age',limits=limit_g,high_age=low_age,high_metalicity=best_metalicty,iso_type=iso_type,interpolate_scale=interpolate_scale)
iso_high_age=isochrone(name,high_age,best_metalicty,best_extinction,special='high_fit_age',limits=limit_g,high_age=high_age,high_metalicity=best_metalicty,iso_type=iso_type,interpolate_scale=interpolate_scale)
iso_all_age=all_isochrones(0,iso_bestfit , iso_low_age, iso_high_age,limits=limit_g)

iso_low_metalicty=isochrone(name, best_age, low_metalicty, best_extinction,special='low_fit_metalicty',limits=limit_g,high_age=best_age,high_metalicity=low_metalicty,iso_type=iso_type,interpolate_scale=interpolate_scale)
iso_high_metalicty=isochrone(name,best_age,high_metalicty,best_extinction,special='high_fit_metalicty',limits=limit_g,high_age=best_age,high_metalicity=high_metalicty,iso_type=iso_type,interpolate_scale=interpolate_scale)
iso_all_metalicty=all_isochrones(0,iso_bestfit , iso_low_metalicty, iso_high_metalicty,limits=limit_g)

iso_low_extinction=isochrone(name, best_age, best_metalicty, low_extinction,special='low_fit_extinction',limits=limit_g,high_age=best_age,high_metalicity=best_metalicty,iso_type=iso_type,interpolate_scale=interpolate_scale)
iso_high_extinction=isochrone(name, best_age, best_metalicty, high_extinction,special='high_fit_extinction',limits=limit_g,high_age=best_age,high_metalicity=best_metalicty,iso_type=iso_type,interpolate_scale=interpolate_scale)
iso_all_extinction=all_isochrones(0,iso_bestfit , iso_low_extinction, iso_high_extinction,limits=limit_g)

iso_low=isochrone(name, low_age, low_metalicty, low_extinction,special='low_fit',limits=limit_g,high_age=low_age,high_metalicity=low_metalicty)
iso_high=isochrone(name, high_age, high_metalicty, high_extinction,special='high_fit',limits=limit_g,high_age=high_age,high_metalicity=high_metalicty)

iso_gaussian=all_isochrones(5000,iso_bestfit , iso_low, iso_high,limits=limit_g,gaussian=True)#

if name=='NGC_2682':
    #delete row 542
    data.remove_row(542)
temperature_getter(data,iso_gaussian.isochrone[5])
# data=mainPart(data, iso_gaussian.isochrone[115])
# step_size=1
# for x in np.arange(0,len(iso_gaussian.isochrone),step_size):
# 	print('starting a loop at '+str(x))
# 	with Pool(processes=step_size) as pool:
# 		results_gaussian=pool.map(partial(temperature_getter,data), iso_gaussian.isochrone[x:x+step_size])

# np.save('temp gaussian save.npy',results_gaussian)

# print('curating results')

# # data.sort('source_id')
# data_gaussian=copy.deepcopy(data)


# teff=[]
# logg=[]
# mass=[]
# e_teff=[]
# e_logg=[]
# e_mass=[]
# for value in range(len(results_gaussian[1][0])):
#     teff.append([x[0][value] for x in results_gaussian ])
#     e_teff.append([x[1][value] for x in results_gaussian ])
#     logg.append([x[2][value] for x in results_gaussian ])
#     e_logg.append([x[3][value] for x in results_gaussian ])
#     mass.append([x[4][value] for x in results_gaussian ])
#     e_mass.append([x[5][value] for x in results_gaussian ])
# new_data_loop=[np.vstack((x1,x2,x3,x4,x5,x6)) for (x1,x2,x3,x4,x5,x6) in zip(teff,e_teff,logg,e_logg,mass,e_mass)]

# def weighted_parameter(parameter_raw,e_parameter):
#     weighted_mean=np.sum(parameter_raw*(1/e_parameter**2))/np.sum((1/e_parameter**2))
#     var_weighted=np.sum((weighted_mean-parameter_raw)**2*(1/e_parameter**2))/np.sum((1/e_parameter**2))
#     total_error=np.sqrt(var_weighted+np.mean(e_parameter)**2)
#     return weighted_mean,total_error

# def currating_weighted(loop):
#     teff_raw=loop[0]
#     logg_raw=loop[2]
#     e_teff_raw=loop[1]
#     e_logg_raw=loop[3] 
#     mass_raw=loop[4] 
#     e_mass_raw=loop[5] 
#     weighted_teff,e_weighted_teff=weighted_parameter(teff_raw,e_teff_raw)
#     weighted_logg,e_weighted_logg=weighted_parameter(logg_raw,e_logg_raw)
#     weighted_mass,e_weighted_mass=weighted_parameter(mass_raw,e_mass_raw)
#     return weighted_teff,e_weighted_teff,weighted_logg,e_weighted_logg, weighted_mass,e_weighted_mass,teff_raw,e_teff_raw,logg_raw,\
#             e_logg_raw,mass_raw,e_mass_raw,
# with Pool(processes=144) as pool:
#     results_gaussian=pool.map(currating_weighted, new_data_loop)


# parameters=['teff','e_teff','logg','e_logg','mass','e_mass','teff_raw','e_teff_raw','logg_raw','e_logg_raw','mass_raw','e_mass_raw']
# for value,param in enumerate(parameters):
#     data_gaussian[param]=[x[value] for x in results_gaussian]
    



# votable=from_table(data_gaussian)
# writeto(votable,name+"_bp_rp_run_2.xml")


# votable=from_table(data_gaussian)
# writeto(votable,name+"_bp_rp_run.xml")

# data['teff']=results[0]
# data['logg']=results[2]
# # data_temp=curating_results(results,data)
# data=data[9:11]
# a=temperature_getter(data, iso_gaussian.isochrone[289])
# for value,iso in enumerate(iso_gaussian.isochrone,288):
#     # print(str(value*8)+' until '+ str((value+1)*8))
#     print(value)
#     temperature_getter(data,iso_gaussian.isochrone[value])
    # iso=iso_gaussian.isochrone[value*8:(value+1)*8]
    # with Pool(processes=8) as pool:
    #     results_gaussian=pool.map(partial(temperature_getter,data),iso )
# data=data[:10]
# data_new=mainPart(data,iso_bestfit.isochrone)
# # data=mainPart(data, iso_bestfit.isochrone)
# # # plt.figure()
# # plt.scatter(iso_all_age.parameters[:,0],iso_all_age.parameters[:,2])
# # plt.xlabel('age')
# # plt.ylabel('differences')


# # plt.figure()
# # plt.scatter(iso_all_metalicty.parameters[:,1],iso_all_metalicty.parameters[:,2])
# # plt.xlabel('metalicty')
# # plt.ylabel('differences')

# # plt.figure()
# # plt.hist(iso_all_metalicty.parameters[:,2])

# # plot_hr(iso_all_metalicty,data,name=name+' metalicity range',save=True,title=r'Metalicty')
# # plot_hr(iso_all_age,data,name=name+' age range',save=True,title=r'Age')
# # plot_hr(iso_all_extinction,data,name=name+' extinction range',save=True,title=r'Extinction')


# # print('calculating temperatures')
# # with Pool(processes=10) as pool:
# # #     results_age=pool.map(partial(temperature_getter,data), iso_all_age.isochrone)
# # data=data[1:3]
# # temperature_getter(data,iso_gaussian.isochrone[1])
# with Pool(processes=30) as pool:
#     results_gaussian=pool.map(partial(temperature_getter,data), iso_gaussian.isochrone)


# # temperature_getter(iso_blanco.isochrone[0],data)
# # temperature_getter(iso_blanco.isochrone[0],data)
# print('calculating temperatures')
# with Pool(processes=10) as pool:
#     results_metalicty=pool.map(partial(temperature_getter,data), iso_all_metalicty.isochrone)
# with Pool(processes=10) as pool:
#     results_metalicty=pool.map(partial(temperature_getter,data), iso_gaussian.isochrone)

# #only uses the stars in the area defined in the file  


# # votable = from_table(data)
# data_age=copy.deepcopy(data)
# data_metalicty=copy.deepcopy(data)
# data_age=curating_results(results_age,data_age)
# # data_metalicty=curating_results(results_metalicty,data_metalicty)
# data.sort('source_id')
# data_gaussian=copy.deepcopy(data)
# data_gaussian=curating_results(results_gaussian,data_gaussian)
# votable=from_table(data_gaussian)
# writeto(votable,name+"_bp_rp_run.xml")

# plot_hr_teff(iso_bestfit, data_gaussian)

# number=2
def plot_pdf(number,data):
    data=data[number]
    plt.figure()
    temperature_star=data['teff_raw']
    logg_star=data['logg_raw']
    teff_error=data['e_teff_raw']
    logg_error=data['e_logg_raw']
    # plt.scatter(temperature_star,logg_star)
    plt.errorbar(temperature_star,logg_star,xerr=teff_error, yerr=logg_error,fmt='o',ms=0.2,elinewidth=0.4,zorder=0)
    xx=data['teff_array']
    yy=data['logg_array']
    zz=data['probability_grid']
    plt.contour(xx,yy,zz)
    plt.xlabel('teff')
    plt.ylabel('logg')   
    plt.gca().invert_yaxis()
    plt.gca().invert_xaxis()


# # teff1=3312
# # logg1=4.7169
# # teff2=3259.4
# # logg2=4.6682

# # m=(logg1-logg2)/(teff1-teff2)
# # c=logg2-m*teff2
def inner_plotter(axes,star,corners,fig):
    plt.rc('font',size=10)

    left, bottom, width, height = axes
    ax2 = fig.add_axes([left, bottom, width, height])
    
    con = ConnectionPatch(xyA=(star['bp_rp'],star['abs_phot_g_mean_mag']), coordsA=ax.transData,
                          xyB=(corners[0],corners[1]), coordsB=ax2.transAxes)
    fig.add_artist(con)
    con = ConnectionPatch(xyA=(star['bp_rp'],star['abs_phot_g_mean_mag']), coordsA=ax.transData,
                          xyB=(corners[2],corners[3]), coordsB=ax2.transAxes)
    fig.add_artist(con)
    
    temperature_star=star['teff_raw']
    logg_star=star['logg_raw']
    teff_error=star['e_teff_raw']
    logg_error=star['e_logg_raw']
    ax2.scatter(temperature_star,logg_star,s=0.1)
    # ax2.errorbar(temperature_star,logg_star,xerr=teff_error, yerr=logg_error,fmt='o',ms=0.2,elinewidth=0.4,zorder=0)
    xx=star['teff_array']
    yy=star['logg_array']
    zz=star['probability_grid']
    ax2.contour(xx,yy,zz)
    ax2.set_xlabel(r'$T_{\rm{eff}}$')
    ax2.set_ylabel(r'$\rm{log}(g)$') 
    
# from matplotlib.patches import ConnectionPatch
# iso_bestfit=iso_bestfit.isochrone
# name='NGC_2682_final_run.xml'
# votable = parse(name)
# data_small_run=votable.get_first_table().to_table(use_names_over_ids=True)
# # plot_hr_teff(iso_all_age,data)
# name=re.sub('_',' ', name)
# data_small_run.sort(['abs_phot_g_mean_mag'])

# fig=plt.figure(figsize=(16.0,9.0))
# plt.rc('font',size=15)
# plt.title(r'NGC 2682')
# ax=fig.gca()
# ax.tick_params(direction='in')
# plt.plot(iso_bestfit[:,4]-iso_bestfit[:,5],iso_bestfit[:,3],label=r'Best fit',c='Blue',linewidth=1.5)
# sigma=[2.5/np.log(10)/data['phot_g_mean_flux_over_error'],2.5/np.log(10)/data['phot_bp_mean_flux_over_error'],2.5/np.log(10)/data['phot_rp_mean_flux_over_error']]
# sigmaDCluster=5/np.log(10)*(data['r_sigma_cluster'])/data['r_est_cluster']
# # 
# errorYCluster=np.sqrt(sigmaDCluster**2+sigma[0]**2)
# errorX=np.sqrt(sigma[2]**2+sigma[1]**2)

# limitsX=[0.53,max(data_small_run['bp_rp'])+0.1]
# limitsY=[min(data_small_run['abs_phot_g_mean_mag']),max(data_small_run['abs_phot_g_mean_mag'])]


# # plt.errorbar(data['bp_rp'],data['abs_phot_g_mean_mag'],xerr=errorX, yerr=errorYCluster,fmt='o',ms=1,elinewidth=1,zorder=3,c='black' )
# plt.scatter(data_small_run['bp_rp'],data_small_run['abs_phot_g_mean_mag'],s=0.3,c='black' )

# plt.xlabel(r'$G_{\rm{BP}}-G_{\rm{RP}}$')
# plt.ylabel(r'$M_{\rm{G}}$')
# data

# plt.xlim(limitsX)
# plt.ylim(limitsY)
# plt.gca().invert_yaxis()
# plt.tight_layout()
# inner_plotter([0.15, 0.70, 0.2, 0.25],data_small_run[len(data_small_run)//30],[0,0,0,1],fig)

# inner_plotter([0.05, 0.10, 0.2, 0.25],data_small_run[len(data_small_run)//10],[0,1,1,1],fig)

# # inner_plotter([0.1, 0.70, 0.2, 0.25],data_small_run[len(data_small_run)//10+1],[0,0,1,0],fig)
# # inner_plotter([0.10, 0.170, 0.17, 0.21],data_small_run[len(data_small_run)*3//10],[1,0,1,1],fig)
# # inner_plotter([0.80, 0.270, 0.17, 0.21],data_small_run[len(data_small_run)*4//10],[1,0,0,0],fig)
# inner_plotter([0.40, 0.400, 0.2, 0.25],data_small_run[len(data_small_run)*5//10],[0,0,1,0],fig)

# plt.savefig(name+' pdf examples.pdf',dpi=100,format='pdf')
# # inner_plotter( [0.25, 0.65,0.15, 0.2],data[len(data)//10],[1,0,0,0],fig)

# # inner_plotter( [0.35, 0.55,0.15, 0.2],data_small_run[len(data_small_run)*6//10],[1,0,0,0],fig)

# # inner_plotter( [0.15, 0.55,0.15, 0.2],data[2],[1,0,0,0])

# # inner_plotter( [0.45, 0.20, 0.15, 0.2],data_small_run[len(data_small_run)*7//10],[0,1,1,1],fig)
# inner_plotter( [0.65, 0.35, 0.2, 0.25],data_small_run[len(data_small_run)*8//10],[1,0,0,0],fig)
# # inner_plotter( [0.65, 0.25, 0.15, 0.2],data[len(data)*5//10],[1,0,0,0],fig)
# # inner_plotter( [0.75, 0.20, 0.15, 0.2],data[len(data)*6//10],[0,1,1,1],fig)
# # inner_plotter( [0.85, 0.35, 0.15, 0.2],data[len(data)*7//10],[1,0,0,0],fig)
# # inner_plotter( [0.95, 0.25, 0.15, 0.2],data[len(data)*8//10],[1,0,0,0],fig)





# fig.tight_layout()
