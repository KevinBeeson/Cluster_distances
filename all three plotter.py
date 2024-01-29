from astropy.table import Table
import csv
import matplotlib.pyplot as plt
from astropy.io.votable import parse
import numpy as np
import os

#gets the names and the parameters of the stars
list_stars = list(csv.reader(open('new_criteria_gaia_distances2.txt', 'rt'), delimiter='\t'))
list_stars = list_stars[1::]
list_stars = np.array(list_stars)
print('Input which cluster you would like to see from the list bellow')
print(list_stars[:,0])
# for name in list_stars[:,0]:
#     print(name)
#     try:
name=input()
path=os.path.abspath(__file__)
#iso = iso[1::]
#iso = np.array(iso,dtype=float)
#
##fig = plt.figure()
##ax = fig.add_subplot(111, projection='3d')
#
##crops the isochrone down to initial mass, int_imf, log_T, G, G_Bp_bright, G_BP_faint, G_RP
#cIso=np.vstack((iso[:,3],iso[:,4],iso[:,7],iso[:,25],iso[:,26],iso[:,27],iso[:,28]))
#cIso=np.transpose(cIso)

opacity=0.2
thick=1.5

position=int(np.where(list_stars==name)[0])

#opens the data for the specific star you have chosen
votable = parse(name+"Final.xml")
data=votable.get_first_table().to_table(use_names_over_ids=True)

#Changes the magnitudes of the stars to absolute magnitudes
Absmag=data['phot_g_mean_mag']-5*(np.log10(data['r_estNew'])-1)
data['abs_phot_g_mean_mag']=Absmag
Absmag=data['phot_rp_mean_mag']-5*(np.log10(data['r_estNew'])-1)
data['abs_phot_rp_mean_mag']=Absmag
Absmag=data['phot_bp_mean_mag']-5*(np.log10(data['r_estNew'])-1)
data['abs_phot_bp_mean_mag']=Absmag


sigma=[2.5/np.log(10)/data['phot_g_mean_flux_over_error'],2.5/np.log(10)/data['phot_bp_mean_flux_over_error'],2.5/np.log(10)/data['phot_rp_mean_flux_over_error']]
sigmaDCluster=5/np.log(10)*(data['r_sigma'])/data['r_estNew']
errorYCluster=np.sqrt(sigmaDCluster**2+sigma[0]**2)
errorX=np.sqrt(sigma[2]**2+sigma[1]**2)


votable=parse(name+"_panstar.xml")
dataPanstar=votable.get_first_table().to_table(use_names_over_ids=True)


Absmag=dataPanstar['g_mean_psf_mag']-5*(np.log10(dataPanstar['r_estNew'])-1)
dataPanstar['abs_g_mean_psf_mag']=Absmag
Absmag=dataPanstar['r_mean_psf_mag']-5*(np.log10(dataPanstar['r_estNew'])-1)
dataPanstar['abs_r_mean_psf_mag']=Absmag
Absmag=dataPanstar['i_mean_psf_mag']-5*(np.log10(dataPanstar['r_estNew'])-1)
dataPanstar['abs_i_mean_psf_mag']=Absmag
Absmag=dataPanstar['z_mean_psf_mag']-5*(np.log10(dataPanstar['r_estNew'])-1)
dataPanstar['abs_z_mean_psf_mag']=Absmag
Absmag=dataPanstar['y_mean_psf_mag']-5*(np.log10(dataPanstar['r_estNew'])-1)
dataPanstar['abs_y_mean_psf_mag']=Absmag

sigmaPanstar=[dataPanstar[ 'g_mean_psf_mag_error'],dataPanstar[ 'r_mean_psf_mag_error'],dataPanstar[ 'i_mean_psf_mag_error'],dataPanstar[ 'z_mean_psf_mag_error'],dataPanstar[ 'y_mean_psf_mag_error']]
sigmaDCluster=5/np.log(10)*(dataPanstar['r_sigma'])/dataPanstar['r_estNew']
x=0
dataPanstarNew=Table(dataPanstar[0],names=(dataPanstar.colnames),meta={'name':'shorter star data'},masked=True)
del dataPanstarNew[0]
sigmaPanstarNew=[]
sigmaDClusterNew=[]
for y in np.transpose(sigmaPanstar):
    if not np.any(np.isnan(y)):
        sigmaPanstarNew.append(y)
        sigmaDClusterNew.append(sigmaDCluster[x])
        dataPanstarNew.add_row(dataPanstar[x])
    x+=1
sigmaDCluster=np.array(sigmaDClusterNew)  
sigmaPanstar=np.transpose(sigmaPanstarNew)
errorYPanstar=np.array(np.sqrt(sigmaPanstar**2+sigmaDCluster**2))
errorXPanstar=np.array(np.sqrt(sigmaPanstar[0]**2+sigmaPanstar[4]**2))
sigmaPanstarNew=[dataPanstarNew[ 'g_mean_psf_mag_error']/dataPanstarNew['g_mean_psf_mag'],dataPanstarNew[ 'r_mean_psf_mag_error']/dataPanstarNew['r_mean_psf_mag'],dataPanstarNew[ 'i_mean_psf_mag_error']/dataPanstarNew['i_mean_psf_mag'],dataPanstarNew[ 'z_mean_psf_mag_error']/dataPanstarNew['z_mean_psf_mag'],dataPanstarNew[ 'y_mean_psf_mag_error']/dataPanstarNew['y_mean_psf_mag']]



votable=parse(name+"_wise.xml")
dataWise=votable.get_first_table().to_table(use_names_over_ids=True)



#Changes the magnitudes of the stars to absolute magnitudes
Absmag=dataWise['phot_g_mean_mag']-5*(np.log10(dataWise['r_estNew'])-1)
dataWise['abs_phot_g_mean_mag']=Absmag
Absmag=dataWise['phot_rp_mean_mag']-5*(np.log10(dataWise['r_estNew'])-1)
dataWise['abs_phot_rp_mean_mag']=Absmag
Absmag=dataWise['phot_bp_mean_mag']-5*(np.log10(dataWise['r_estNew'])-1)
dataWise['abs_phot_bp_mean_mag']=Absmag
Absmag=dataWise['w1mpro']-5*(np.log10(dataWise['r_estNew'])-1)
dataWise['abs_w1mpro']=Absmag
Absmag=dataWise['w2mpro']-5*(np.log10(dataWise['r_estNew'])-1)
dataWise['abs_w2mpro']=Absmag
Absmag=dataWise['w3mpro']-5*(np.log10(dataWise['r_estNew'])-1)
dataWise['abs_w3mpro']=Absmag
Absmag=dataWise['w4mpro']-5*(np.log10(dataWise['r_estNew'])-1)
dataWise['abs_w4mpro']=Absmag
Absmag=dataWise['w1gmag']-5*(np.log10(dataWise['r_estNew'])-1)
dataWise['abs_w1gmag']=Absmag
Absmag=dataWise['w2gmag']-5*(np.log10(dataWise['r_estNew'])-1)
dataWise['abs_w2gmag']=Absmag
Absmag=dataWise['w3gmag']-5*(np.log10(dataWise['r_estNew'])-1)
dataWise['abs_w3gmag']=Absmag
Absmag=dataWise['w4gmag']-5*(np.log10(dataWise['r_estNew'])-1)
dataWise['abs_w4gmag']=Absmag


sigmaWise=[dataWise[ 'w1mpro_error'],dataWise[ 'w2mpro_error'],dataWise[ 'w3mpro_error'],dataWise[ 'w4mpro_error']]
sigmaDCluster=5/np.log(10)*(dataWise['r_sigma'])/dataWise['r_estNew']
sigmaWiseGaia=[2.5/np.log(10)/dataWise['phot_g_mean_flux_over_error'],2.5/np.log(10)/dataWise['phot_bp_mean_flux_over_error'],2.5/np.log(10)/dataWise['phot_rp_mean_flux_over_error']]

x=0
dataWiseNew=Table(dataWise[0],names=(dataWise.colnames),meta={'name':'shorter star data'},masked=True)
del dataWiseNew[0]
sigmaWiseNew=[]
sigmaDClusterNew=[]
sigmaDClusterNewGaia=[]
for y in np.transpose(sigmaWise):
    if not np.any(np.isnan(y)):
        sigmaWiseNew.append(y)
        sigmaDClusterNew.append(sigmaDCluster[x])
        sigmaDClusterNewGaia.append(np.transpose(sigmaWiseGaia)[x])
        dataWiseNew.add_row(dataWise[x])
    x+=1
sigmaDCluster=np.array(sigmaDClusterNew)  
sigmaWise=np.transpose(sigmaWiseNew)
sigmaDClusterNewGaia=np.transpose(np.array(sigmaDClusterNewGaia))

errorYWise=np.array(np.sqrt(sigmaWise**2+sigmaDCluster**2))
errorYWiseGaia=np.array(np.sqrt(sigmaDClusterNewGaia**2+sigmaDCluster**2))
errorXWise=np.array(np.sqrt(sigmaWise[0]**2+sigmaDClusterNewGaia[1]**2))



plotter=np.array([0.0])
age=0.0
clusterDistance=np.mean(data['r_estNew'])
# for x in iso:
#     if x[:,25]+5*(np.log10(clusterDistance)-1)>11.99:
#         faintIso=np.append()




limitsX=[min(data['bp_rp'])-0.02,max(data['bp_rp'])+0.02]
limitsY=[min(data['abs_phot_g_mean_mag'])-0.1,max(data['abs_phot_g_mean_mag'])+0.1]

limitsXWise=[min(dataWise['phot_bp_mean_mag']-dataWise['w1mpro'])-0.02,max(dataWise['phot_bp_mean_mag']-dataWise['w1mpro'])+0.02]
limitsYWise=[min(dataWise['abs_phot_g_mean_mag'])-0.1,max(dataWise['abs_phot_g_mean_mag'])+0.1]

limitsXPanstar=[min(dataPanstar['g_mean_psf_mag']-dataPanstar['y_mean_psf_mag'])-0.02,max(dataPanstar['g_mean_psf_mag']-dataPanstar['y_mean_psf_mag'])+0.02]
limitsYPanstar=[min(dataPanstar['abs_i_mean_psf_mag'])-0.1,max(dataPanstar['abs_i_mean_psf_mag'])+0.1]

extintion='0.0'
iso = list(csv.reader(open(name+'Errors.txt', 'rt'), delimiter=' '))
isoWise = list(csv.reader(open(name+'WiseErrors.txt', 'rt'), delimiter=' '))
isoPanstar = list(csv.reader(open(name+'PanstarErrors.txt', 'rt'), delimiter=' '))

isoMid = list(csv.reader(open(name+'ISO.txt', 'rt'), delimiter=' '))
isoMidWise = list(csv.reader(open(name+'Wise ISO.txt', 'rt'), delimiter=' '))
isoMidPanstar = list(csv.reader(open(name+'Panstar ISO.txt', 'rt'), delimiter=' '))


for x in iso:
    if x[1]=='Attention:':
        extintion=x[6][3:-1]
        break
row=0 
fig, (ax1, ax2, ax3) = plt.subplots(1,3)

for x in iso:
    #Reads all the data first before plotting it
    if x[0]!='#' and x[0]!='#isochrone':
        age=x[2]
        metal=x[1]
        if float(x[25])<25:
            if plotter[0]==0.0:
                plotter=np.array(x)
                plotterWise=np.array(isoWise[row])
                plotterPanstar=np.array(isoPanstar[row])
            else:
                while isoWise[row][0][0]=='#':
                    row+=1
                plotter=np.vstack([plotter,np.array(x)])
                plotterWise=np.vstack([plotterWise,isoWise[row]])
                plotterPanstar=np.vstack([plotterPanstar,isoPanstar[row]])
    elif age!=0.0:
        plotter=plotter.astype('float64')
        plotterWise=plotterWise.astype('float64')
        plotterPanstar=plotterPanstar.astype('float64')
        
        faintIso=np.array([i for i in plotter if (i[25]+5*(np.log10(clusterDistance)-1)>9.99 and i[27]-i[28]<4.0)])
        
        brightIso=np.array([i for i in plotter if i[25]+5*(np.log10(clusterDistance)-1)<11.99])
        brightWise=[]
        faintWise=[]
        for i in range(len(plotter)):
            if (plotter[i][25]+5*(np.log10(clusterDistance)-1)>9.99 and plotter[i][27]-plotter[i][28]<4.0):
                faintWise.append(plotterWise[i])
            if plotter[i][25]+5*(np.log10(clusterDistance)-1)<11.99:
                brightWise.append(plotterWise[i])
                
        brightWise=np.array(brightWise)
        faintWise=np.array(faintWise)
        # [iso[i]so[i] for i in iso if i[25]+5*(np.log10(clusterDistance)-1)>11.99]
        #     faint=np.append(faint,i)
        # [i for i in iso if i[25]+5*(np.log10(clusterDistance)-1)>11.99]
        #     faint=np.append(faint,i)
        if len(brightIso)!=0:
            ax1.plot(brightIso[:,26]-brightIso[:,28],brightIso[:,25],color='r',alpha=opacity,linewidth=thick)
            # plt.scatter(brightIso[:,26]-brightIso[:,28],brightIso[:,25],color='r',s=10)
    
            ax1.plot(brightIso[:,26]-brightIso[:,28],brightIso[:,25]-0.752,color='k',linestyle='dashed',alpha=opacity,linewidth=thick)
        ax1.plot(faintIso[:,27]-faintIso[:,28],faintIso[:,25],color='g',alpha=opacity,linewidth=thick)
        # plt.scatter(faintIso[:,27]-faintIso[:,28],faintIso[:,25],color='g',s=10)

        ax1.plot(faintIso[:,27]-faintIso[:,28],faintIso[:,25]-0.752,color='k',linestyle='dashed',alpha=opacity,linewidth=thick)
#        ax1.scatter(data['bp_rp'],data['abs_phot_g_mean_mag'],color='b',s=4)
        # ax1.errorbar(data['bp_rp'],data['abs_phot_g_mean_mag'],xerr=errorX, yerr=errorYCluster,fmt='o',ms=1,elinewidth=2)
        
#        ax2.scatter(dataWise['w1gmag']-dataWise['w4gmag'],dataWise['abs_w3gmag'],color='b',s=2)
        # ax2.scatter(dataWise['phot_bp_mean_mag']-dataWise['w1mpro'],dataWise['abs_phot_g_mean_mag'],label='without errors',color='g',s=4)
        ax2.plot(faintIso[:,27]-faintWise[:,35],faintIso[:,25],color='g',alpha=opacity,linewidth=thick)
        ax2.plot(faintIso[:,27]-faintWise[:,35],faintIso[:,25]-0.752,color='k',linestyle='dashed',alpha=opacity,linewidth=thick)
        if len(brightWise)!=0:
            ax2.plot(brightIso[:,26]-brightWise[:,35],brightIso[:,25],color='g',alpha=opacity,linewidth=thick)
            ax2.plot(brightIso[:,26]-brightWise[:,35],brightIso[:,25]-0.752,color='k',linestyle='dashed',alpha=opacity,linewidth=thick)

        # ax2.errorbar(dataWiseNew['phot_bp_mean_mag']-dataWiseNew['w1mpro'],dataWiseNew['abs_phot_g_mean_mag'],color='b',label='with errors',xerr=errorXWise, yerr=errorYWiseGaia[0,:],fmt='o',ms=1,elinewidth=2)
        # ax2.legend(loc='best',fontsize=20)

        
        # ax3.scatter(dataPanstar['g_mean_psf_mag']-dataPanstar['y_mean_psf_mag'],dataPanstar['abs_i_mean_psf_mag'],label='without errors',color='g',s=4)
        # ax3.errorbar(dataPanstarNew['g_mean_psf_mag']-dataPanstarNew['y_mean_psf_mag'],dataPanstarNew['abs_i_mean_psf_mag'],color='b',label='with errors',xerr=errorXPanstar, yerr=errorYPanstar[2,:],fmt='o',ms=1,elinewidth=2)
        ax3.plot(plotterPanstar[:,25]-plotterPanstar[:,29],plotterPanstar[:,27],color='g',alpha=opacity,linewidth=thick)
        ax3.plot(plotterPanstar[:,25]-plotterPanstar[:,29],plotterPanstar[:,27]-0.752,color='k',linestyle='dashed',alpha=opacity,linewidth=thick)
        # plt.scatter(faintIso[8,27]-faintIso[8,28],faintIso[8,25],color='k',s=2000,marker='+')

        # plt.scatter(data['bp_rp'],data['abs_phot_g_mean_mag'],color='b',s=2)
        # plt.errorbar(data['bp_rp'],data['abs_phot_g_mean_mag'],xerr=errorX, yerr=errorYCluster,fmt='o',ms=0.1,elinewidth=1.5,color='b')
#        plt.errorbar(data['bp_rp'],data['abs_phot_g_mean_mag'],xerr=errorX, yerr=errorYCluster,fmt='o',ms=1.5,elinewidth=2.5)
#        plt.scatter(dataHermes['bp_rp'],dataHermes['abs_phot_g_mean_mag'],s=5,c='blue',linewidth=15,marker='+',label='Observed by GALAH also')
#        plt.scatter(notHermesStars['bp_rp'],notHermesStars['abs_phot_g_mean_mag'],s=5,c='black',label='Observed by Gaia only')
        # if np.mean(data['phot_g_mean_mag'])<10.99:
        #     plt.plot(plotter[:,26]-plotter[:,28],plotter[:,25],color='r',label='bright')
        #     plt.plot(plotter[:,26]-plotter[:,28],plotter[:,25]-0.752,color='k',label='Binary',linestyle='dashed')

        # else:
        #     plt.plot(plotter[:,27]-plotter[:,28],plotter[:,25],color='g',label='faint')
        #     plt.plot(plotter[:,27]-plotter[:,28],plotter[:,25]-0.752,color='k',label='Binary',linestyle='dashed')


        plt.xlim(limitsX)
        plt.ylim(limitsY)
        plt.title(name+' HR diagram',fontsize=25)
        plt.xlabel('$G_{BP}- G_{RP}$',fontsize=20)
        plt.gca().invert_yaxis()
        plt.ylabel('Absolute G mag',fontsize=20)
        fig.set_size_inches(26,13)
        plt.xticks(fontsize=15)
        plt.yticks(fontsize=15)
        plt.tight_layout(pad=2.0)

        # plt.gca().set_aspect('equal', adjustable='box')
        # plt.savefig('D:\Work\Temperature weighter\\' +name +'Iso\\'+name+' age='+str(age)+' metalicity='+str(metal)+' extinction Av='+str(extintion)+'.png')
        # plt.close(fig)
        age=0.0
        plotter=np.array([0.0])
    row+=1

for x in isoMid:
    if x[1]=='Attention:':
        extintion=x[6][3:-1]
        break                
row=0
for x in isoMid:
    if x[0]!='#' and x[0]!='#isochrone':
        age=x[2]
        metal=x[1]
        if float(x[25])<25:
            if plotter[0]==0.0:
                plotter=np.array(x)
                plotterWise=np.array(isoMidWise[row])
                plotterPanstar=np.array(isoMidPanstar[row])
            else:

                while row<len(isoMidWise) and isoMidWise[row][0][0]=='#':
                    row+=1
                if  row<len(isoMidWise) and isoMidWise[row][1]!='terminated':

                    plotterWise=np.vstack([plotterWise,isoMidWise[row]])
                    plotterPanstar=np.vstack([plotterPanstar,isoMidPanstar[row]])
    
                plotter=np.vstack([plotter,np.array(x)])
    
    elif age!=0.0:
            
        plotter=plotter.astype('float64')
        plotterWise=plotterWise.astype('float64')
        plotterPanstar=plotterPanstar.astype('float64')
        
        

        # [iso[i]so[i] for i in iso if i[25]+5*(np.log10(clusterDistance)-1)>11.99]
        #     faint=np.append(faint,i)
        # fig, (ax1, ax2, ax3) = plt.subplots(1,3)
        if len(brightIso)!=0:
            ax1.plot(brightIso[:,26]-brightIso[:,28],brightIso[:,25],color='r',label='Bright isochrone')
            # plt.scatter(brightIso[:,26]-brightIso[:,28],brightIso[:,25],color='r',s=10)
    
            ax1.plot(brightIso[:,26]-brightIso[:,28],brightIso[:,25]-0.752,color='k',linestyle='dashed')
        ax1.plot(faintIso[:,27]-faintIso[:,28],faintIso[:,25],color='g',label='Faint isochrone')
        # plt.scatter(faintIso[:,27]-faintIso[:,28],faintIso[:,25],color='g',s=10)

        ax1.plot(faintIso[:,27]-faintIso[:,28],faintIso[:,25]-0.752,color='k',label='Binary',linestyle='dashed')
#        ax1.scatter(data['bp_rp'],data['abs_phot_g_mean_mag'],color='b',s=4)
        ax1.errorbar(data['bp_rp'],data['abs_phot_g_mean_mag'],xerr=errorX, yerr=errorYCluster,fmt='o',ms=1,elinewidth=2)
        
        
        
        if len(plotter)!=len(plotterWise):
            plotterNew=[]
            plotterWiseNew=[]
            limit=0.08
            start=0
            
            if len(plotter)>len(plotterWise):
                for count,x in enumerate(plotter):
                    pos=start

                    while pos<len(plotterWise) and (abs(x[3]-plotterWise[pos][3])>limit and abs(x[3]-plotterWise[pos][3])<0.1)  :
                        pos+=1                 
                    if pos<len(plotterWise) and abs(x[3]-plotterWise[pos][3])<limit:
                        plotterNew.append(x)
                        plotterWiseNew.append(plotterWise[pos])
                        start+=1
            elif len(plotterWise)>len(plotter):
                for count,x in enumerate(plotterWise):
                    pos=0
                    while pos<len(plotterWise) and abs(x[3]-plotter[pos][3])>limit and abs(x[3]-plotter[pos][3])<0.1:
                        pos+=1                    
                    if pos<len(plotterWise) and abs(x[3]-plotter[pos][3])<limit:
                        plotterWiseNew.append(x)
                        plotterNew.append(plotter[pos])
                        pos+=1

            plotter=plotterNew
            plotterWise=plotterWiseNew
        faintIso=np.array([i for i in plotter if (i[25]+5*(np.log10(clusterDistance)-1)>9.99 and i[27]-i[28]<4.0)])
        brightIso=np.array([i for i in plotter if i[25]+5*(np.log10(clusterDistance)-1)<11.99])

        brightWise=[]
        faintWise=[]
        for i in range(len(plotter)):
            if (plotter[i][25]+5*(np.log10(clusterDistance)-1)>9.99 and plotter[i][27]-plotter[i][28]<4.0):
                faintWise.append(plotterWise[i])
            if plotter[i][25]+5*(np.log10(clusterDistance)-1)<11.99:
                brightWise.append(plotterWise[i])

        brightWise=np.array(brightWise)
        faintWise=np.array(faintWise)      
        
        
        
#        ax2.scatter(dataWise['w1gmag']-dataWise['w4gmag'],dataWise['abs_w3gmag'],color='b',s=2)
        ax2.scatter(dataWise['phot_bp_mean_mag']-dataWise['w1mpro'],dataWise['abs_phot_g_mean_mag'],label='Without errors',color='orange',s=4)
        ax2.plot(faintIso[:,27]-faintWise[:,35],faintIso[:,25],color='g')
        ax2.plot(faintIso[:,27]-faintWise[:,35],faintIso[:,25]-0.752,color='k',linestyle='dashed')
        if len(brightWise)!=0:
            ax2.plot(brightIso[:,26]-brightWise[:,35],brightIso[:,25],color='g')
            ax2.plot(brightIso[:,26]-brightWise[:,35],brightIso[:,25]-0.752,color='k',label='Binary',linestyle='dashed')

        ax2.errorbar(dataWiseNew['phot_bp_mean_mag']-dataWiseNew['w1mpro'],dataWiseNew['abs_phot_g_mean_mag'],color='b',label='With errors',xerr=errorXWise, yerr=errorYWiseGaia[0,:],fmt='o',ms=1,elinewidth=2)

        
        ax3.scatter(dataPanstar['g_mean_psf_mag']-dataPanstar['y_mean_psf_mag'],dataPanstar['abs_i_mean_psf_mag'],label='Without errors',color='orange',s=4)
        ax3.errorbar(dataPanstarNew['g_mean_psf_mag']-dataPanstarNew['y_mean_psf_mag'],dataPanstarNew['abs_i_mean_psf_mag'],label='With errors',xerr=errorXPanstar, yerr=errorYPanstar[2,:],fmt='o',ms=1,elinewidth=2)
        ax3.plot(plotterPanstar[:,25]-plotterPanstar[:,29],plotterPanstar[:,27],color='g')
        ax3.plot(plotterPanstar[:,25]-plotterPanstar[:,29],plotterPanstar[:,27]-0.752,color='k',label='Binary',linestyle='dashed')
#        plt.errorbar(data['bp_rp'],data['abs_phot_g_mean_mag'],xerr=errorX, yerr=errorYCluster,fmt='o',ms=1.5,elinewidth=2.5)
#        plt.scatter(dataHermes['bp_rp'],dataHermes['abs_phot_g_mean_mag'],s=5,c='blue',linewidth=15,marker='+',label='Observed by GALAH also')
#        plt.scatter(notHermesStars['bp_rp'],notHermesStars['abs_phot_g_mean_mag'],s=5,c='black',label='Observed by Gaia only')
        # if np.mean(data['phot_g_mean_mag'])<10.99:
        #     plt.plot(plotter[:,26]-plotter[:,28],plotter[:,25],color='r',label='bright')
        #     plt.plot(plotter[:,26]-plotter[:,28],plotter[:,25]-0.752,color='k',label='Binary',linestyle='dashed')

        # else:
        #     plt.plot(plotter[:,27]-plotter[:,28],plotter[:,25],color='g',label='faint')
        #     plt.plot(plotter[:,27]-plotter[:,28],plotter[:,25]-0.752,color='k',label='Binary',linestyle='dashed')


        #fig.suptitle(name+' HR diagrams',fontsize=20)

        
        if name=='NGC_2204':
            limitsX=[limitsX[0],2.2]
            limitsY=[limitsY[1]]
            limitsXWise=[limitsXWise[0],3.5]
            limitsXPanstar=[-1,2.2]
        ax1.set_xlim(limitsX)
        ax1.set_ylim(limitsY)
        ax1.invert_yaxis()
        ax1.set_title(' Gaia',fontsize=35)
        ax1.set_ylabel('Absolute G mag',fontsize=25)
        ax1.set_xlabel( '$G_{BP}- G_{RP}$',fontsize=25)
        # ax1.legend(loc='best',fontsize=25)       
        ax1.tick_params(axis='both', labelsize= 20)
        
        ax2.set_title('Allwise ',fontsize=35)
        ax2.set_ylabel('Absolute G mag',fontsize=25)
        ax2.set_xlabel( '$G_{BP}$-W1',fontsize=25)
        ax2.set_ylim(limitsYWise)
        ax2.set_xlim(limitsXWise)
        ax2.invert_yaxis()
        ax2.tick_params(axis='both', labelsize= 20)
        # ax2.legend(loc='best',fontsize=25)
        
        
        ax3.set_title('Panstarrs',fontsize=35)
        ax3.set_ylabel('Absolute I mag ',fontsize=25)
        ax3.set_ylim(limitsYPanstar)
        ax3.set_xlim(limitsXPanstar)
        ax3.set_xlabel( 'G-Y',fontsize=25)
        ax3.invert_yaxis()
        # ax3.legend(loc='best',fontsize=25)

        ax3.tick_params(axis='both', labelsize= 20)
        plt.tight_layout(pad=3.0)
        
        # plt.gca().set_aspect('equal', adjustable='box')
        age=0.0
        # plt.savefig(name+'HR with errors.png')

        plotter=np.array([0.0])      
        fig.savefig(name+'all three.pdf',format='pdf',dpi=3000)
    row+=1
    # except:
    #     print('no data')
plt.tight_layout()       
# for x in iso:
#     #Reads all the data first before plotting it
#     if x[0]!='#' and x[0]!='#isochrone':
#         age=x[2]
#         metal=x[1]
#         if float(x[25])<25:
#             if plotter[0]==0.0:
#                 plotter=np.array(x)
#                 plotterWise=np.array(isoWise[row])
#                 plotterPanstar=np.array(isoPanstar[row])
#             else:
#                 plotter=np.vstack([plotter,np.array(x)])
#                 plotterWise=np.vstack([plotterWise,isoWise[row]])
#                 plotterPanstar=np.vstack([plotterPanstar,isoPanstar[row]])
#     elif age!=0.0:
#         plotter=plotter.astype('float64')
#         plotterWise=plotterWise.astype('float64')
#         plotterPanstar=plotterPanstar.astype('float64')
        
#         faintIso=np.array([i for i in plotter if (i[25]+5*(np.log10(clusterDistance)-1)>9.99 and i[27]-i[28]<4.0)])
#         brightIso=np.array([i for i in plotter if i[25]+5*(np.log10(clusterDistance)-1)<11.99])


#         brightWise=[]
#         faintWise=[]
#         for i in range(len(plotter)):
#             if (plotter[i][25]+5*(np.log10(clusterDistance)-1)>9.99 and plotter[i][27]-plotter[i][28]<4.0):
#                 faintWise.append(plotterWise[i])
#             if plotter[i][25]+5*(np.log10(clusterDistance)-1)<11.99:
#                 brightWise.append(plotterWise[i])

#         brightWise=np.array(brightWise)
#         faintWise=np.array(faintWise)
#         # [iso[i]so[i] for i in iso if i[25]+5*(np.log10(clusterDistance)-1)>11.99]
#         #     faint=np.append(faint,i)
#         fig, (ax1, ax2, ax3) = plt.subplots(1,3)
#         if len(brightIso)!=0:
#             ax1.plot(brightIso[:,26]-brightIso[:,28],brightIso[:,25],color='r',label='bright')
#             # plt.scatter(brightIso[:,26]-brightIso[:,28],brightIso[:,25],color='r',s=10)
    
#             ax1.plot(brightIso[:,26]-brightIso[:,28],brightIso[:,25]-0.752,color='k',label='Bright Binary',linestyle='dashed')
#         ax1.plot(faintIso[:,27]-faintIso[:,28],faintIso[:,25],color='g',label='faint')
#         # plt.scatter(faintIso[:,27]-faintIso[:,28],faintIso[:,25],color='g',s=10)

#         ax1.plot(faintIso[:,27]-faintIso[:,28],faintIso[:,25]-0.752,color='k',label='Faint Binary',linestyle='dashed')
# #        ax1.scatter(data['bp_rp'],data['abs_phot_g_mean_mag'],color='b',s=4)
#         ax1.errorbar(data['bp_rp'],data['abs_phot_g_mean_mag'],xerr=errorX, yerr=errorYCluster,fmt='o',ms=1,elinewidth=2)
        
# #        ax2.scatter(dataWise['w1gmag']-dataWise['w4gmag'],dataWise['abs_w3gmag'],color='b',s=2)
#         ax2.scatter(dataWise['phot_bp_mean_mag']-dataWise['w1mpro'],dataWise['abs_phot_g_mean_mag'],label='without errors',color='g',s=4)
#         ax2.plot(faintIso[:,27]-faintWise[:,35],faintIso[:,25],color='r')
#         ax2.plot(faintIso[:,27]-faintWise[:,35],faintIso[:,25]-0.752,color='k',linestyle='dashed')
#         ax2.plot(brightIso[:,26]-brightWise[:,35],brightIso[:,25],color='r')
#         ax2.plot(brightIso[:,26]-brightWise[:,35],brightIso[:,25]-0.752,color='k',label='Binary',linestyle='dashed')

#         ax2.errorbar(dataWiseNew['phot_bp_mean_mag']-dataWiseNew['w1mpro'],dataWiseNew['abs_phot_g_mean_mag'],color='b',label='with errors',xerr=errorXWise, yerr=errorYWiseGaia[0,:],fmt='o',ms=1,elinewidth=2)
#         ax2.legend(loc='best',fontsize=20)

        
#         ax3.scatter(dataPanstar['g_mean_psf_mag']-dataPanstar['y_mean_psf_mag'],dataPanstar['abs_i_mean_psf_mag'],label='without errors',color='g',s=4)
#         ax3.errorbar(dataPanstarNew['g_mean_psf_mag']-dataPanstarNew['y_mean_psf_mag'],dataPanstarNew['abs_i_mean_psf_mag'],color='b',label='with errors',xerr=errorXPanstar, yerr=errorYPanstar[2,:],fmt='o',ms=1,elinewidth=2)
#         ax3.plot(plotterPanstar[:,25]-plotterPanstar[:,29],plotterPanstar[:,27],color='r')
#         ax3.plot(plotterPanstar[:,25]-plotterPanstar[:,29],plotterPanstar[:,27]-0.752,color='k',label='Binary',linestyle='dashed')

#        # if np.mean(data['phot_g_mean_mag'])<10.99:
#         #     plt.plot(plotter[:,26]-plotter[:,28],plotter[:,25],color='r',label='bright')
#         #     plt.plot(plotter[:,26]-plotter[:,28],plotter[:,25]-0.752,color='k',label='Binary',linestyle='dashed')

#         # else:
#         #     plt.plot(plotter[:,27]-plotter[:,28],plotter[:,25],color='g',label='faint')
#         #     plt.plot(plotter[:,27]-plotter[:,28],plotter[:,25]-0.752,color='k',label='Binary',linestyle='dashed')
# #        plt.legend(loc='best',fontsize=25)
# #        plt.title(name+' HR diagram',fontsize=25)
# #        ax1.set_xlabel('$\G_{BP}-$\G_{RP}',fontsize=20)

        
        
# #        plt.tight_layout()
# #        # plt.gca().set_aspect('equal', adjustable='box')
# #        plt.savefig('/home/kevin/Documents/Temperature weighter/'+name+'Iso/'+name +' age='+age+' metalicity='+metal+' extinction Av='+extintion+'.png')
# #        plt.close(fig)
#         fig.set_size_inches(26,13)
#         fig.tight_layout(pad=3.0)
#         fig.subplots_adjust(top=0.88)

#         fig.show()

#         age=0.0
#        plotter=np.array([0.0])    
        
        
#plots after the crop
#plt.scatter(dataShort['pmra'],dataShort['pmdec'],s=1)
#plt.figure()
#plt.show()
#
#plt.scatter(dataShort['abs_phot_bp_mean_mag']-dataShort['abs_phot_rp_mean_mag'],dataShort['abs_phot_g_mean_mag'],s=1,color='r')
#plt.scatter(iso[:,26]-iso[:,28],iso[:,25],color='b',s=1)
#plt.gca().invert_yaxis()
#plt.figure()
#plt.show()
            
#fig = plt.figure()
#plt.scatter(dataShort['abs_phot_bp_mean_mag']-dataShort['abs_phot_rp_mean_mag'],dataShort['abs_phot_g_mean_mag'],s=1,color='r')
#plt.scatter(cIso[:,4]-cIso[:,6],cIso[:,3],color='b',s=1)
#plt.gca().invert_yaxis()
#
#plt.show()
#
#
#fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')
#ax.scatter(dataShort['abs_phot_g_mean_mag'],dataShort['abs_phot_bp_mean_mag'],dataShort['abs_phot_rp_mean_mag'],s=1,color='r')
#ax.scatter(cIso[:,3],cIso[:,4],cIso[:,6],color='b',s=1)
#plt.show()