#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 16 15:33:39 2022

@author: kevin
"""
from astropy.io.votable import parse,from_table,writeto

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import kde
from pathlib import Path
import requests
import re
from os.path import exists
import csv
from scipy.interpolate import interp1d
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
        if isinstance(self.isochrone, np.ndarray) and iso_type=='gaia' and interpolate_scale!=1:
            self.isochrone=interpolate(self.isochrone,interpolate_scale)
        elif iso_type=='gaia'and interpolate_scale!=1:
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
def plot_hr(isochrones,data,axis,save=False,paper=True,name=False,title=False):
    if not name:
        name=isochrones.name
    name=re.sub('_',' ', name)
    plt.rc('font', size=15)
    iso=isochrones.isochrone
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
        
plt.rcParams['font.family']='Times New Roman'
plt.rcParams['xtick.direction']='in'
plt.rcParams['ytick.direction']='in'
plt.rcParams['text.usetex'] = True
cluster_details_all = list(csv.reader(open('dr3_clusters.txt', 'rt'), delimiter=','))
cluster_details_all = np.array(cluster_details_all)
# name=input('cluster')
# cluster_details=[x for x in cluster_details_all if x[0]==name][0]

# low_age=float(cluster_details[1])
# low_metalicty=float(cluster_details[2])
# low_extinction=float(cluster_details[3])
# best_age=float(cluster_details[4])
# best_metalicty=float(cluster_details[5])
# best_extinction=float(cluster_details[6])
# high_age=float(cluster_details[7])
# high_metalicty=float(cluster_details[8])
# high_extinction=float(cluster_details[9])

iso_type='gaia'
interpolate_scale=1
# #opens the data for the specific star you have chosen
# if iso_type=='gaia':
#     votable = parse(name+"_dr3_distances3.xml")
# elif iso_type=='allwise':
#     votable = parse(name+"_allwise.xml")
# elif iso_type=='panstarrs':
#     votable = parse(name+"_panstarrs.xml")
# data=votable.get_first_table().to_table(use_names_over_ids=True)
# iso_bestfit=isochrone(name,best_age,best_metalicty,best_extinction,special='Best_fit',high_age=best_age,high_metalicity=best_metalicty,iso_type=iso_type,interpolate_scale=interpolate_scale)
# plot_hr(iso_bestfit, data)

shape=(3,3)
ax=[]
fig = plt.figure(figsize=shape)
for row in range(shape[0]):
    for col in range(shape[1]):
        ax.append(plt.subplot2grid(shape=shape,loc=(row,col)))
        
for value,(axis,name) in enumerate(zip(ax[:len(cluster_details_all)],cluster_details_all[:,0])):
    print(name)
    cluster_details=[x for x in cluster_details_all if x[0]==name][0]
    name_title=re.sub('_',' ', name)

    low_age=float(cluster_details[1])
    low_metalicty=float(cluster_details[2])
    low_extinction=float(cluster_details[3])
    best_age=float(cluster_details[4])
    best_metalicty=float(cluster_details[5])
    best_extinction=float(cluster_details[6])
    high_age=float(cluster_details[7])
    high_metalicty=float(cluster_details[8])
    high_extinction=float(cluster_details[9])

    iso_type='gaia'
    #opens the data for the specific star you have chosen
    if iso_type=='gaia':
        votable = parse(name+"_dr3_distances3.xml")
    elif iso_type=='allwise':
        votable = parse(name+"_allwise.xml")
    elif iso_type=='panstarrs':
        votable = parse(name+"_panstarrs.xml")
    data=votable.get_first_table().to_table(use_names_over_ids=True)
    Absmag=data['phot_g_mean_mag']-5*(np.log10(data['r_est_cluster'])-1)
    data['abs_phot_g_mean_mag']=Absmag

    iso_bestfit=isochrone(name,best_age,best_metalicty,best_extinction,special='Best_fit',high_age=best_age,high_metalicity=best_metalicty,iso_type=iso_type,interpolate_scale=interpolate_scale)
    iso=iso_bestfit.isochrone
    axis.plot(iso[:,4]-iso[:,5],iso[:,3],label=r'Best fit',c='blue',linewidth=1.5,alpha=0.7)
    axis.scatter(data['bp_rp'],data['abs_phot_g_mean_mag'],s=0.3,zorder=0,c='black' ,alpha=0.9)
    axis.set_xlim((max(min(data['bp_rp']),min(iso[:,4]-iso[:,5]))-0.1,max(data['bp_rp'])))
    axis.set_ylim((min(data['abs_phot_g_mean_mag']),max(data['abs_phot_g_mean_mag'])))    
    name_title = " ".join(word[0].upper()+word[1:] for word in name_title.split(" "))
    if name_title=='NGC 2682':
        name_title='M67'
    axis.set_title(name_title)
    axis.invert_yaxis()
    if value==0 :
        axis.set_ylabel(r'$M_{\rm{G}}$')
    if value==6:
        axis.set_xlabel(r'$\rm{G}_{\rm{BP}}-\rm{G}_{\rm{RP}}$')
# ax[-1].xaxis.set_visible(False)
# ax[-1].yaxis.set_visible(False)

# for spine in ['top', 'right', 'left', 'bottom']:
#    ax[-1].spines[spine].set_visible(False)
fig.set_size_inches(3.32088003321*2,6)
plt.tight_layout(pad=0.5)

plt.savefig('test_test.pdf')
