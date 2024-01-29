from scipy import stats
import csv
from astropy import units as u
from astropy.coordinates import SkyCoord
from astroquery.gaia import Gaia
import numpy as np
from astropy.table import Table
from astropy.io.votable import from_table,writeto,parse
from astropy.table import join
import matplotlib.pyplot as plt
from astropy.io import fits
from astroquery.simbad import Simbad


#c = SkyCoord(ra=280, dec=-60, unit=(u.degree, u.degree), frame='icrs')
#radius = u.Quantity(2.73, u.deg)
#j = Gaia.cone_search_async(c, radius)
#r = j.get_results()
  #        dist.r_est, dist.r_lo, dist.r_hi, dist.r_len, dist.result_flag, dist.modality_flag,gaia_source.visibility_periods_used,gaia_source.parallax ,gaia_source.pmra,gaia_source.pmdec,gaia_source.source_id,gaia_source.ra,gaia_source.ra_error,gaia_source.dec,gaia_source.dec_error,gaia_source.parallax_error,gaia_source.phot_g_mean_mag,gaia_source.bp_rp,gaia_source.radial_velocity,gaia_source.radial_velocity_error,gaia_source.phot_variable_flag,gaia_source.teff_val,gaia_source.a_g_val\



textBlocks=["SELECT TOP 300000\
          gaia_source.phot_variable_flag,gaia_source.parallax,gaia_source.parallax_error,dist.r_est,dist.r_lo,dist.r_hi,dist.r_len,gaia_source.source_id,gaia_source.ra,gaia_source.dec,gaia_source.phot_g_mean_mag,gaia_source.phot_bp_mean_mag,gaia_source.phot_rp_mean_mag,gaia_source.phot_g_mean_flux_over_error,gaia_source.bp_rp,phot_rp_mean_flux_over_error,gaia_source.phot_bp_mean_flux_over_error,dist.r_est,gaia_source.pmra,gaia_source.pmdec,gaia_source.bp_rp\
          FROM gaiadr2.gaia_source\
          INNER JOIN external.gaiadr2_geometric_distance AS dist \
              ON gaiadr2.gaia_source.source_id=dist.source_id\
          WHERE CONTAINS(POINT('ICRS',gaiadr2.gaia_source.ra,gaiadr2.gaia_source.dec),CIRCLE('ICRS',","))=1\
          AND  (gaiadr2.gaia_source.visibility_periods_used>=7)\
          AND pmra IS NOT NULL\
          AND pmdec IS NOT NULL\
          AND bp_rp IS NOT NULL\
          AND phot_g_mean_flux_over_error>25\
          AND phot_rp_mean_flux_over_error>10\
          AND phot_bp_mean_flux_over_error>10\
          AND phot_bp_rp_excess_factor < 1.3+0.06*power(phot_bp_mean_mag-phot_rp_mean_mag,2)\
          AND phot_bp_rp_excess_factor > 1.0+0.015*power(phot_bp_mean_mag-phot_rp_mean_mag,2)\
          AND  astrometric_chi2_al/(astrometric_n_good_obs_al-5)<1.44*greatest(1,exp(-0.4*(phot_g_mean_mag-19.5)))\
          AND result_flag = 1\
          AND parallax>"," AND parallax <"]

          # INNER JOIN gaiaedr3.dr2_neighbourhood AS neighbourhood \
          #     ON gaiaedr3.source_id=neighbourhood.dr3_source_id\
          # INNER JOIN external.gaiadr2_geometric_distance AS dist \
          #     ON neighbourhood.dr2_neighbourhood=dist.source_id\

textBlocksDR3=["SELECT TOP 300000\
          gaia_source.dr2_radial_velocity,gaia_source.ra_error,gaia_source.dec_error,gaia_source.pmra_error,gaia_source.pmdec_error,gaia_source.dr2_radial_velocity_error,dist.r_est,dist.r_len,r_hi,r_lo,neighbourhood.dr2_source_id,gaia_source.pseudocolour,gaia_source.parallax,gaia_source.nu_eff_used_in_astrometry,gaia_source.ecl_lat,gaia_source.astrometric_params_solved,gaia_source.parallax_error,gaia_source.source_id,gaia_source.ra,gaia_source.dec,gaia_source.phot_g_mean_mag,gaia_source.phot_bp_mean_mag,gaia_source.phot_rp_mean_mag,gaia_source.phot_g_mean_flux_over_error,gaia_source.bp_rp,phot_rp_mean_flux_over_error,gaia_source.phot_bp_mean_flux_over_error,gaia_source.pmra,gaia_source.pmdec,gaia_source.bp_rp\
          FROM gaiaedr3.gaia_source\
          INNER JOIN gaiaedr3.dr2_neighbourhood AS neighbourhood \
              ON gaia_source.source_id=neighbourhood.dr3_source_id\
          INNER JOIN external.gaiadr2_geometric_distance AS dist \
              ON neighbourhood.dr2_source_id=dist.source_id\
          WHERE CONTAINS(POINT('ICRS',gaiaedr3.gaia_source.ra,gaiaedr3.gaia_source.dec),CIRCLE('ICRS',","))=1\
          AND  (gaiaedr3.gaia_source.visibility_periods_used>=7)\
          AND pmra IS NOT NULL\
          AND pmdec IS NOT NULL\
          AND bp_rp IS NOT NULL\
          AND phot_g_mean_flux_over_error>25\
          AND phot_rp_mean_flux_over_error>10\
          AND phot_bp_mean_flux_over_error>10\
          AND phot_bp_rp_excess_factor < 1.3+0.06*power(phot_bp_mean_mag-phot_rp_mean_mag,2)\
          AND phot_bp_rp_excess_factor > 1.0+0.015*power(phot_bp_mean_mag-phot_rp_mean_mag,2)\
          AND  astrometric_chi2_al/(astrometric_n_good_obs_al-5)<1.44*greatest(1,exp(-0.4*(phot_g_mean_mag-19.5)))\
          AND parallax>"," AND parallax <"]
textBlocksDR3_dist=["SELECT TOP 300000\
          dist.r_med_geo,dist.r_lo_geo,dist.r_hi_geo,dist.r_med_photogeo,dist.r_lo_photogeo,dist.r_hi_photogeo,gaia_source.source_id,gaia_source.ra,gaia_source.ra_error,gaia_source.dec,gaia_source.dec_error,gaia_source.parallax,gaia_source.parallax_error,gaia_source.parallax_over_error,gaia_source.pm,gaia_source.pmra,gaia_source.ecl_lat,gaia_source.pmra_error,gaia_source.pmdec,gaia_source.pmdec_error,gaia_source.astrometric_params_solved,gaia_source.nu_eff_used_in_astrometry,gaia_source.pseudocolour,gaia_source.pseudocolour_error,gaia_source.visibility_periods_used,gaia_source.phot_g_mean_flux,gaia_source.phot_g_mean_flux_error,gaia_source.phot_g_mean_flux_over_error,gaia_source.phot_g_mean_mag,gaia_source.phot_bp_mean_flux,gaia_source.phot_bp_mean_flux_error,gaia_source.phot_bp_mean_flux_over_error,gaia_source.phot_bp_mean_mag,gaia_source.phot_rp_mean_flux,gaia_source.phot_rp_mean_flux_error,gaia_source.phot_rp_mean_flux_over_error,gaia_source.phot_rp_mean_mag,gaia_source.phot_bp_rp_excess_factor,gaia_source.bp_rp,gaia_source.radial_velocity,gaia_source.radial_velocity_error,gaia_source.vbroad,gaia_source.phot_variable_flag,gaia_source.l,gaia_source.b,gaia_source.teff_gspphot,gaia_source.teff_gspphot_lower,gaia_source.teff_gspphot_upper,gaia_source.logg_gspphot,gaia_source.logg_gspphot_lower,gaia_source.logg_gspphot_upper,gaia_source.mh_gspphot,gaia_source.mh_gspphot_lower \
          FROM gaiadr3.gaia_source\
          INNER JOIN external.gaiaedr3_distance AS dist \
              ON gaia_source.source_id=dist.source_id\
          WHERE CONTAINS(POINT('ICRS',gaiadr3.gaia_source.ra,gaiadr3.gaia_source.dec),CIRCLE('ICRS',","))=1\
          AND  (gaiadr3.gaia_source.visibility_periods_used>=7)\
          AND pmra IS NOT NULL\
          AND pmdec IS NOT NULL\
          AND bp_rp IS NOT NULL\
          AND phot_g_mean_flux_over_error>25\
          AND phot_rp_mean_flux_over_error>10\
          AND phot_bp_mean_flux_over_error>10\
          AND phot_bp_rp_excess_factor < 1.3+0.06*power(phot_bp_mean_mag-phot_rp_mean_mag,2)\
          AND phot_bp_rp_excess_factor > 1.0+0.015*power(phot_bp_mean_mag-phot_rp_mean_mag,2)\
          AND  astrometric_chi2_al/(astrometric_n_good_obs_al-5)<1.44*greatest(1,exp(-0.4*(phot_g_mean_mag-19.5)))\
          AND parallax>"," AND parallax <"]
# textBlocks=["SELECT TOP 300000\
#           gaia_source
#           FROM gaiadr2.gaia_source\
#           INNER JOIN external.gaiadr2_geometric_distance AS dist \
#               ON gaiadr2.gaia_source.source_id=dist.source_id\
#           WHERE CONTAINS(POINT('ICRS',gaiadr2.gaia_source.ra,gaiadr2.gaia_source.dec),CIRCLE('ICRS',","))=1\
#           AND  (gaiadr2.gaia_source.visibility_periods_used>=7)\
#           AND pmra IS NOT NULL\
#           AND pmdec IS NOT NULL\
#           AND bp_rp IS NOT NULL\
#           AND phot_g_mean_flux_over_error>50\
#           AND phot_rp_mean_flux_over_error>20\
#           AND phot_bp_mean_flux_over_error>20\
#           AND phot_bp_rp_excess_factor < 1.3+0.06*power(phot_bp_mean_mag-phot_rp_mean_mag,2)\
#           AND phot_bp_rp_excess_factor > 1.0+0.015*power(phot_bp_mean_mag-phot_rp_mean_mag,2)\
#           AND  astrometric_chi2_al/(astrometric_n_good_obs_al-5)<1.44*greatest(1,exp(-0.4*(phot_g_mean_mag-19.5)))\
#           AND result_flag = 1\
#           AND parallax>"," AND parallax <"]
textBlocksShort=["SELECT TOP 20000\
          dist.r_est,dist.r_lo,dist.r_hi,gaia_source.source_id,gaia_source.phot_g_mean_mag,gaia_source.phot_bp_mean_mag,gaia_source.phot_rp_mean_mag,gaia_source.phot_g_mean_flux_over_error,gaia_source.bp_rp,phot_rp_mean_flux_over_error,gaia_source.phot_bp_mean_flux_over_error,dist.r_est,gaia_source.pmra,gaia_source.pmdec,gaia_source.bp_rp\
          FROM gaiadr2.gaia_source\
          INNER JOIN external.gaiadr2_geometric_distance AS dist \
              ON gaiadr2.gaia_source.source_id=dist.source_id\
          WHERE CONTAINS(POINT('ICRS',gaiadr2.gaia_source.ra,gaiadr2.gaia_source.dec),CIRCLE('ICRS',","))=1\
          AND pmra IS NOT NULL\
          AND pmdec IS NOT NULL\
          AND bp_rp IS NOT NULL\
          AND parallax>"," AND parallax <"]

list_stars = list(csv.reader(open('new_criteria_gaia_distances2.txt', 'rt'), delimiter='\t'))
#list_stars = list(csv.reader(open('Cluster_list_n_body.txt', 'rt'), delimiter=','))


list_stars = list_stars[1::]
list_stars = np.array(list_stars)

list_mine=list(csv.reader(open('my_targets.txt', 'rt'), delimiter='\t'))
print('Input which cluster you would like to see from the list bellow')
print(list_stars[:,0])
# name=input()
for name in list_stars[:,0]:
    # name=name[0]

    # votable = parse(name+"_dr3.xml")
    # data_old=votable.get_first_table().to_table(use_names_over_ids=True)

    results_table =Simbad.query_object(name,u.deg)
    print(results_table)
    # back_up_data=list(csv.reader(open('no_simbad_open_clusters.txt', 'rt'), delimiter='\t'))    
    # back_up_data = back_up_data[1::]
    # back_up_data = np.array(back_up_data)
    
    if results_table is None:
        position=np.where(back_up_data[:,0]==name)
        ra=back_up_data[position][0][1]
        dec=back_up_data[position][0][2]
        coord=SkyCoord(ra=ra,dec=dec,unit=(u.hourangle, u.deg),frame='icrs')    
    else:
        coord=SkyCoord(ra=results_table['RA'][0],dec=results_table['DEC'][0],unit=(u.hourangle, u.deg),frame='icrs')
    # ra='17.033'
    # dec='-16 17 00'
    # coord=SkyCoord(ra=ra,dec=dec,unit=(u.hourangle, u.deg),frame='icrs')
    
    # coord=SkyCoord(ra=results_table['RA'][0],dec=results_table['DEC'][0],unit=(u.hourangle, u.deg),frame='icrs')
    print(coord)
    
    
    position=int(np.where(list_stars==name)[0])
    request=textBlocksDR3_dist[0]+str(coord.ra.degree) +' , ' + str(coord.dec.degree) +','+ str(list_stars[position,9]) +textBlocksDR3_dist[1] + str(float(list_stars[position,7])-float(list_stars[position,8])) + textBlocksDR3_dist[2] + str(float(list_stars[position,7])+float(list_stars[position,8]))
    
    job = Gaia.launch_job_async(request,output_format="votable",output_file=(name+" dr3.xml"),dump_to_file=False)
    data=job.get_results()
    #
    # shorterns the data to the data around a specified pmra and pmdec
    dataShort=Table(data[0],names=(data.colnames),meta={'name':'shorter star data'},masked=True)
    del dataShort[0]
    for star in range(len(data)):
        if data['pmra'][star] >float(list_stars[position,3])-float(list_stars[position,4]) and data['pmra'][star] <float(list_stars[position,3])+float(list_stars[position,4]):
            if data['pmdec'][star] >float(list_stars[position,5])-float(list_stars[position,6]) and data['pmdec'][star] <float(list_stars[position,5])+float(list_stars[position,6]):
                dataShort.add_row(data[star])
    
    # votable = from_table(dataShort)
    # writeto(votable,name+"_dr3.xml")
    
    # plt.figure()
    # plt.scatter(data['pmra'],data['pmdec'],s=0.1)
    
    # plt.axvline(x=float(list_stars[position,3])-float(list_stars[position,4]),color='r')
    # plt.axvline(x=float(list_stars[position,4])+float(list_stars[position,3]),color='r')
    # plt.axhline(y=float(list_stars[position,5])-float(list_stars[position,6]),color='r')
    # plt.axhline(y=float(list_stars[position,5])+float(list_stars[position,6]),color='r')
    # plt.title(name)
    # plt.xlabel('pmra')
    # plt.ylabel('pmdec')
    # plt.ylim([float(list_stars[position,5])-float(list_stars[position,6])*2,float(list_stars[position,5])+float(list_stars[position,6])*2])

    # plt.xlim([float(list_stars[position,3])-float(list_stars[position,4])*2,float(list_stars[position,3])+float(list_stars[position,4])*2])
    # plt.savefig('/home/kevin/Documents/Temperature_weighter_2/plots/'+name+'_pmra_pmdec',dpi=300)
    
    # plt.show()
    # print(len(dataShort),' len corrected ')
    # print(len(data_old),' len old')
#import dr6.0.fits



# request2=textBlocksShort[0]+str(coord.ra.degree) +' , ' + str(coord.dec.degree) +','+ str(list_stars[position,9]) +textBlocksShort[1] + str(float(list_stars[position,7])-float(list_stars[position,8])) + textBlocksShort[2] + str(float(list_stars[position,7])+float(list_stars[position,8]))
# job2 = Gaia.launch_job_async(request2,output_format="votable",output_file=(name+"2UFages.xml"),dump_to_file=False)
# # data2=job2.get_results()
# print(len(dataShort))
# distance=''
# for x in data:
#     distance+=str(x['ra'])+'\t'+str(x['dec'])+'\n'

# file=open(name+"distance.txt",'w')
# file.write(distance)
#job = Gaia.launch_job_async(clusters[name][0],output_format="votable",output_file=(name+".xml"),dump_to_file=True)
#data=job.get_results()

#dataShort=Table(data[0],names=(data.colnames),meta={'name':'shorter star data'},masked=True)
#del dataShort[0]
#for star in range(len(data)):
#    if data['pmra'][star] >float(list_stars[position,3])-float(list_stars[position,4])/2. and data['pmra'][star] <float(list_stars[position,3])+float(list_stars[position,4])/2.:
#        if data['pmdec'][star] >float(list_stars[position,5])-float(list_stars[position,6])/2. and data['pmdec'][star] <float(list_stars[position,5])+float(list_stars[position,6])/2.:
#            dataShort.add_row(data[star])
#            
#dataShortUF=Table(data[0],names=(data.colnames),meta={'name':'shorter star data'},masked=True)
#del dataShortUF[0]
#for star in range(len(data2)):
#    if data2['pmra'][star] >float(list_stars[position,3])-float(list_stars[position,4])/2. and data2['pmra'][star] <float(list_stars[position,3])+float(list_stars[position,4])/2.:
#        if data2['pmdec'][star] >float(list_stars[position,5])-float(list_stars[position,6])/2. and data2['pmdec'][star] <float(list_stars[position,5])+float(list_stars[position,6])/2.:
#            dataShortUF.add_row(data2[star])


# print('Making Graphs')
# dataShort2=Table(data[0],names=(data.colnames),meta={'name':'shorter2 star data'},masked=True)
# del dataShort2[0]
# for star in range(len(dataShort)):
#     if not( dataShort['r_est'][star]+2*(dataShort['r_hi'][star]-dataShort['r_est'][star]) <float(list_stars[position,-1])*0.8-20 or dataShort['r_est'][star]-2*(dataShort['r_est'][star]-dataShort['r_lo'][star])>float(list_stars[position,-1])*1.2+20) :
#         dataShort2.add_row(dataShort[star])
# dataShort=dataShort2
#g1=np.array(dataShort['phot_g_mean_mag'])
#bp1=np.array(dataShort['phot_bp_mean_mag'])
#rp1=np.array(dataShort['phot_rp_mean_mag'])
#bp_rp1=np.array(dataShort['bp_rp'])
##photometric=np.vstack([g,bp,rp])
#
#photometric=np.vstack([bp1,g1,rp1])
#kernel=stats.gaussian_kde(photometric)
#
#g=np.array(dataShortUF['phot_g_mean_mag'])
#bp=np.array(dataShortUF['phot_bp_mean_mag'])
#rp=np.array(dataShortUF['phot_rp_mean_mag'])
#bp_rp=np.array(dataShortUF['bp_rp'])
#
##photometric=np.vstack([g,bp,rp])
#
#photometric=np.vstack([bp,g,rp])
#
#kernelUF=stats.gaussian_kde(photometric)


#density=kernelUF([g,bp,rp])/kernel([g,bp,rp])
#density=np.divide(kernelUF([bp1,g1,rp1]),kernel([bp1,g1,rp1]))
#dataShort['density']=density
#
#vot=from_table(dataShort)
#writeto(vot,name+"ages.xml")



#gMin=g.min()
#gMax=g.max()
#rpMin=rp.min()
#rpMax=rp.max()
#bpMin=bp.min()
#bpMax=bp.max()
#bp_rpMin=bp_rp.min()
#bp_rpMax=bp_rp.max()
#
#BP_RP,G = np.mgrid[ bp_rpMin:bp_rpMax+1:100j,gMin:gMax+1:100j]
#positions=np.vstack([BP_RP.ravel(),G.ravel()])
#Z1 = np.reshape(kernel(positions).T, BP_RP.shape) 
#Z2=np.reshape(kernelUF(positions).T, BP_RP.shape) 
#fig, ax = plt.subplots()
#ax.imshow(np.rot90(np.divide(Z1,Z2)), cmap=plt.cm.gist_earth_r,
#          extent=[ bp_rpMin, bp_rpMax+1,gMin, gMax+1])
#
#ax.scatter(dataShort['bp_rp'],dataShort['phot_g_mean_mag'],color='r',s=1)
#ax.set_xlim([bp_rpMin, bp_rpMax+1])
#ax.set_ylim([gMin, gMax+1])
##ax.scatter(dataShortUF['bp_rp'],dataShortUF['phot_g_mean_mag'],color='b',s=0.1)
#plt.show()

#G, BP , RP= np.mgrid[gMin:gMax:100j, bpMin:bpMax:100j, rpMin:rpMax:100j]
#
#positions=np.vstack([G.ravel(),BP.ravel(),RP.ravel()])
#
#
#Z1 = np.reshape(kernel(positions).T, G.shape) 
#Z2=np.reshape(kernelUF(positions).T, G.shape) 
#fig, ax = plt.subplots()
#ax.imshow(np.rot90(Z), cmap=plt.cm.gist_earth_r,
#          extent=[xmin, xmax, ymin, ymax])
#ax.plot( bp_rp,g, 'k.', markersize=2)
#ax.set_xlim([xmin, xmax])
#ax.set_ylim([ymin, ymax])
#plt.show()

#Absmag=data['phot_g_mean_mag']-5*(np.log10(data['r_est'])-1)
#data['abs_phot_g_mean_mag']=Absmag
#
#pl
# ##t.scatter(data['pmra'],data['pmdec'],s=1)



# plt.figure()
# plt.hist(data['r_est'],bins=100)
# plt.hist(dataShort['r_est'],bins=100)
# plt.show()
#
#plt.scatter(data['bp_rp'],data['abs_phot_g_mean_mag'],s=1)
#plt.gca().invert_yaxis()
#plt.figure()
#plt.show()
#
#
#dataShort=Table(data[0],names=(data.colnames),meta={'name':'shorter star data'},masked=True)
#del dataShort[0]
#for star in range(len(data)):
#    if data['pmra'][star] >clusters[name][1] and data['pmra'][star] <clusters[name][2]:
#        if data['pmdec'][star] >  clusters[name][3]and data['pmdec'][star]<clusters[name][4]:
#            dataShort.add_row(data[star])
#            
#            
#plt.scatter(dataShort['pmra'],dataShort['pmdec'],s=1)
#plt.figure()
#plt.show()
#
#plt.scatter(dataShort['bp_rp'],dataShort['abs_phot_g_mean_mag'],s=1)
#plt.gca().invert_yaxis()
#plt.figure()
#plt.show()