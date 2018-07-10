#!/usr/bin/env python
"""
------------------------------------------------------------------------

GRB Catalog cross referencing
------------------------------------------------------------------------
"""

import urllib
from astropy.io import fits,ascii
from astropy.table import Table
import numpy as np
from astropy.time import Time
from astropy.coordinates import SkyCoord
from astropy import units as u
import matplotlib.pylab as plot
from astropy.cosmology import WMAP9 as cosmo
from scipy.optimize import curve_fit
import scipy.odr as odr
from matplotlib.ticker import LogLocator,MultipleLocator
from scipy import interpolate

def logline(beta,x):

	m=beta[0]
	b=beta[1]
#	logx=np.log10(x)
#	f=10**b*x**m
	f=b*(x/1e52)**m

	return f

# def part_dev(eng,model,alpha=alpha,beta=beta,epeak=epeak,epiv=100):

# 	dfda=[] ; dfdb=[] ; dfdep=[]
# 	if model=='pl':
# 		dfda=-(eng/epiv)**(-alpha)*np.log(eng/epiv)
# 	if model=='cpl':
# 		dfda=1/epeak*(eng/epiv)**alpha*exp(-eng/epeak*(2+alpha))*(epeak*np.log(eng/epiv)-eng)
# 		dfdep=(2+alpha)*eng/epeak**2*(eng/epiv)**alpha*exp(-(2+alpha)*eng/epeak)

def dist2z(dist):
	zs=np.arange(0.01,10,0.01)
	dists=cosmo.luminosity_distance(zs).value
	f=interpolate.interp1d(dists,zs,bounds_error=False,fill_value="extrapolate",kind='linear')
	z=f(dist)

	return z

def kcorr(input_Emin,input_Emax,output_Emin,output_Emax,gamma,z):

	eng=np.logspace(np.log10(output_Emin),np.log10(output_Emax),100)/(1.+z)
	eng2=np.logspace(np.log10(input_Emin),np.log10(input_Emax),100)

	f1=pl(eng,gamma,epiv=1)
	f2=pl(eng2,gamma,epiv=1)

	k1=np.trapz(f1*eng,eng)
	k2=np.trapz(f2*eng2,eng2)

	k=k1/k2

	return k

def pl(eng,alpha,epiv=100):

	f=(eng/epiv)**alpha

	return f

def sbpl(eng,alpha,epeak,beta,epiv=100,delta=0.3):

	m=(beta-alpha)/2
	b=(alpha+beta)/2
	q=np.log10(eng/epeak)/delta
	qpiv=np.log10(epiv/epeak)/delta
	apiv=m*delta*np.log((np.exp(qpiv)+np.exp(-qpiv))/2)
	a=m*delta*np.log((np.exp(q)+np.exp(-q))/2)

	f=(eng/epiv)**b*10**(a-apiv)

	return f

def comp(eng,alpha,epeak,epiv=100):
  
  f=(eng/epiv)**alpha*np.exp(-eng*(2.+alpha)/epeak)

  return f

def band(eng,alpha,epeak,beta,enorm=100):

  w1=np.where(eng <= (alpha-beta)*epeak/(2+alpha))
  w2=np.where(eng >= (alpha-beta)*epeak/(2+alpha))
  f=np.zeros(len(eng))
  if len(w1[0]) > 0:
  	f[w1]=(eng[w1]/enorm)**alpha*np.exp(-(2+alpha)*eng[w1]/epeak)
  if len(w2[0]) > 0:
  	f[w2]=((alpha-beta)*epeak/enorm/(2+alpha))**(alpha-beta)*np.exp(beta-alpha)*(eng[w2]/enorm)**beta

  return f

def bbody(eng,kt):

	f=8.0525*eng**2/(kt**4*(np.exp(eng/kt)-1))

	return f

def match_catalogs_name(name1,name2):

	ind_dict = dict((k,i) for i,k in enumerate(name1))
	inter = set(ind_dict).intersection(name2)
	m1 = [ ind_dict[x] for x in inter ]

	ind_dict = dict((k,i) for i,k in enumerate(name2))
	inter = set(ind_dict).intersection(name1)
	m2 = [ ind_dict[x] for x in inter ]

	return m1,m2

def match_catalogs_time_coord(ra1,dec1,met1,ra2,dec2,met2):

	c=SkyCoord(ra=ra1*u.deg,dec=dec1*u.deg)
	d=SkyCoord(ra=ra2*u.deg,dec=dec2*u.deg)

	m1=[]
	m2=[]
	for i in range(len(c)):

		dist=c[i].separation(d)
		fsep=abs(met1[i]-met2)
		w=np.where((fsep < 60) & (dist<50.*u.deg))
		if len(w[0])>0:
			m1=np.append(m1,i)
			m2=np.append(m2,w[0])

	m1=m1.astype('int')
	m2=m2.astype('int')

	return m1,m2


def load_GBM(dir='/Users/jracusin/GBM/'):

	gbm=fits.open(dir+'gbmgrbcat.fits')
	gbm=gbm[1].data
	cols=gbm.columns

	met0=Time('2001-01-01 00:00:00',format='iso',scale='utc')
#	colname=[re.split(";|=|'",str(col))[2] for col in cols]
	
	#gbmcol=[gbm.field(c) for c in range(len(colname))]
	# doesn't work

	mo=[] ; alpha=[] ; alpha_neg_err=[] ; alpha_pos_err=[]
	epeak=[] ; epeak_neg_err=[] ; epeak_pos_err=[]
	beta=[] ; beta_neg_err=[] ; beta_pos_err=[]
	ra=[] ; dec=[] ; t90=[] ; t90_err=[] ; t90_start=[]
	fluence=[] ; fluence_err=[] 
	gbmname=[] ; trigtime=[] ; met=[]
	pflux_mo=[] ; pflux=[] ; pflux_err=[] ; pflux_alpha=[] 
	pflux_alpha_neg_err=[] ; pflux_alpha_pos_err=[]
	pflux_epeak=[] ; pflux_epeak_neg_err=[] ; pflux_epeak_pos_err=[]
	pflux_beta=[] ; pflux_beta_neg_err=[] ; pflux_beta_pos_err=[]
	pfluxph64=[] ; pfluxph256=[] ; pfluxph1024=[]
	i=0

	for j in range(len(gbm)):
		gbmname=np.append(gbmname,gbm.NAME[i]) 
		pflux_mo=np.append(pflux_mo,gbm.PFLX_BEST_FITTING_MODEL[i].strip())
		mo=np.append(mo,gbm.FLNC_BEST_FITTING_MODEL[i].strip())
		ra=np.append(ra,float(gbm.RA[i]))
		dec=np.append(dec,float(gbm.DEC[i]))

		if gbm.T90[i] != '        ':
			t90=np.append(t90,float(gbm.T90[i]))
			t90_err=np.append(t90_err,float(gbm.T90_ERROR[i]))
			t90_start=np.append(t90_start,float(gbm.T90_START[i]))
		else:
			t90=np.append(t90,np.nan)
			t90_err=np.append(t90_err,np.nan)
			t90_start=np.append(t90_start,np.nan)
		if gbm.FLUENCE[i] != '          ':
			fluence=np.append(fluence,float(gbm.FLUENCE[i]))
			fluence_err=np.append(fluence_err,float(gbm.FLUENCE_ERROR[i]))
#			flux=np.append(flux,float(gbm.FLUX_1024[i]))
#			flux_err=np.append(flux_err,float(gbm.FLUX_1024_ERROR[i]))
		else:
			fluence=np.append(fluence,0.)
			fluence_err=np.append(fluence_err,0)
#			flux=np.append(flux,0.)
#			flux_err=np.append(flux_err,0)
		if gbm.FLUX_64[i] != '         ':
			pfluxph64=np.append(pfluxph64,float(gbm.FLUX_64[i]))
		else: pfluxph64=np.append(pfluxph64,np.nan)
		if gbm.FLUX_256[i] != '         ':
			pfluxph256=np.append(pfluxph256,float(gbm.FLUX_256[i]))
		else: pfluxph256=np.append(pfluxph256,np.nan)
		if gbm.FLUX_1024[i] != '         ':
			pfluxph1024=np.append(pfluxph1024,float(gbm.FLUX_1024[i]))
		else: pfluxph1024=np.append(pfluxph1024,np.nan)

		if 'BAND' in gbm.FLNC_BEST_FITTING_MODEL[i]:
			alpha=np.append(alpha,gbm.FLNC_BAND_ALPHA[i])
			alpha_neg_err=np.append(alpha_neg_err,gbm.FLNC_BAND_ALPHA_NEG_ERR[i])
			alpha_pos_err=np.append(alpha_pos_err,gbm.FLNC_BAND_ALPHA_POS_ERR[i])
			epeak=np.append(epeak,gbm.FLNC_BAND_EPEAK[i])
			epeak_neg_err=np.append(epeak_neg_err,gbm.FLNC_BAND_EPEAK_NEG_ERR[i])
			epeak_pos_err=np.append(epeak_pos_err,gbm.FLNC_BAND_EPEAK_POS_ERR[i])
			beta=np.append(beta,gbm.FLNC_BAND_BETA[i])
			beta_neg_err=np.append(beta_neg_err,gbm.FLNC_BAND_BETA_NEG_ERR[i])
			beta_pos_err=np.append(beta_pos_err,gbm.FLNC_BAND_BETA_POS_ERR[i])

		elif 'COMP' in gbm.FLNC_BEST_FITTING_MODEL[i]:
			alpha=np.append(alpha,gbm.FLNC_COMP_INDEX[i])
			alpha_neg_err=np.append(alpha_neg_err,gbm.FLNC_COMP_INDEX_NEG_ERR[i])
			alpha_pos_err=np.append(alpha_pos_err,gbm.FLNC_COMP_INDEX_POS_ERR[i])
			epeak=np.append(epeak,gbm.FLNC_COMP_EPEAK[i])
			epeak_neg_err=np.append(epeak_neg_err,gbm.FLNC_COMP_EPEAK_NEG_ERR[i])
			epeak_pos_err=np.append(epeak_pos_err,gbm.FLNC_COMP_EPEAK_POS_ERR[i])
			beta=np.append(beta,0.)
			beta_neg_err=np.append(beta_neg_err,0.)
			beta_pos_err=np.append(beta_pos_err,0.)
		elif 'SBPL' in gbm.FLNC_BEST_FITTING_MODEL[i]:
			alpha=np.append(alpha,gbm.FLNC_SBPL_INDX1[i])
			alpha_neg_err=np.append(alpha_neg_err,gbm.FLNC_SBPL_INDX1_NEG_ERR[i])
			alpha_pos_err=np.append(alpha_pos_err,gbm.FLNC_SBPL_INDX1_POS_ERR[i])
			epeak=np.append(epeak,gbm.FLNC_SBPL_BRKEN[i])
			epeak_neg_err=np.append(epeak_neg_err,gbm.FLNC_SBPL_BRKEN_NEG_ERR[i])
			epeak_pos_err=np.append(epeak_pos_err,gbm.FLNC_SBPL_BRKEN_POS_ERR[i])
			beta=np.append(beta,gbm.FLNC_SBPL_INDX2[i])
			beta_neg_err=np.append(beta_neg_err,gbm.FLNC_SBPL_INDX2_NEG_ERR[i])
			beta_pos_err=np.append(beta_pos_err,gbm.FLNC_SBPL_INDX2_POS_ERR[i])
		elif 'PLAW' in gbm.FLNC_BEST_FITTING_MODEL[i]:
			alpha=np.append(alpha,gbm.FLNC_PLAW_INDEX[i])
			alpha_neg_err=np.append(alpha_neg_err,gbm.FLNC_PLAW_INDEX_NEG_ERR[i])
			alpha_pos_err=np.append(alpha_pos_err,gbm.FLNC_PLAW_INDEX_POS_ERR[i])
			epeak=np.append(epeak,0.)
			epeak_neg_err=np.append(epeak_neg_err,0.)
			epeak_pos_err=np.append(epeak_pos_err,0.)
			beta=np.append(beta,0.)
			beta_neg_err=np.append(beta_neg_err,0.)
			beta_pos_err=np.append(beta_pos_err,0.)
		else:
			alpha=np.append(alpha,0.)
			alpha_neg_err=np.append(alpha_neg_err,0.)
			alpha_pos_err=np.append(alpha_pos_err,0.)
			epeak=np.append(epeak,0.)
			epeak_neg_err=np.append(epeak_neg_err,0.)
			epeak_pos_err=np.append(epeak_pos_err,0.)
			beta=np.append(beta,0.)
			beta_neg_err=np.append(beta_neg_err,0.)
			beta_pos_err=np.append(beta_pos_err,0.)

		if 'BAND' in gbm.PFLX_BEST_FITTING_MODEL[i]:
			pflux=np.append(pflux,float(gbm.PFLX_BAND_ERGFLUX[i]))
			pflux_err=np.append(pflux_err,float(gbm.PFLX_BAND_ERGFLUX_ERROR[i]))
			pflux_alpha=np.append(pflux_alpha,gbm.PFLX_BAND_ALPHA[i])
			pflux_alpha_neg_err=np.append(pflux_alpha_neg_err,gbm.PFLX_BAND_ALPHA_NEG_ERR[i])
			pflux_alpha_pos_err=np.append(pflux_alpha_pos_err,gbm.PFLX_BAND_ALPHA_POS_ERR[i])
			pflux_epeak=np.append(pflux_epeak,gbm.PFLX_BAND_EPEAK[i])
			pflux_epeak_neg_err=np.append(pflux_epeak_neg_err,gbm.PFLX_BAND_EPEAK_NEG_ERR[i])
			pflux_epeak_pos_err=np.append(pflux_epeak_pos_err,gbm.PFLX_BAND_EPEAK_POS_ERR[i])
			pflux_beta=np.append(pflux_beta,gbm.PFLX_BAND_BETA[i])
			pflux_beta_neg_err=np.append(pflux_beta_neg_err,gbm.PFLX_BAND_BETA_NEG_ERR[i])
			pflux_beta_pos_err=np.append(pflux_beta_pos_err,gbm.PFLX_BAND_BETA_POS_ERR[i])

		elif 'COMP' in gbm.PFLX_BEST_FITTING_MODEL[i]:
			pflux=np.append(pflux,float(gbm.PFLX_COMP_ERGFLUX[i]))
			pflux_err=np.append(pflux_err,float(gbm.PFLX_COMP_ERGFLUX_ERROR[i]))
			pflux_alpha=np.append(pflux_alpha,gbm.PFLX_COMP_INDEX[i])
			pflux_alpha_neg_err=np.append(pflux_alpha_neg_err,gbm.PFLX_COMP_INDEX_NEG_ERR[i])
			pflux_alpha_pos_err=np.append(pflux_alpha_pos_err,gbm.PFLX_COMP_INDEX_POS_ERR[i])
			pflux_epeak=np.append(pflux_epeak,gbm.PFLX_COMP_EPEAK[i])
			pflux_epeak_neg_err=np.append(pflux_epeak_neg_err,gbm.PFLX_COMP_EPEAK_NEG_ERR[i])
			pflux_epeak_pos_err=np.append(pflux_epeak_pos_err,gbm.PFLX_COMP_EPEAK_POS_ERR[i])
			pflux_beta=np.append(pflux_beta,0.)
			pflux_beta_neg_err=np.append(pflux_beta_neg_err,0.)
			pflux_beta_pos_err=np.append(pflux_beta_pos_err,0.)

		elif 'SBPL' in gbm.PFLX_BEST_FITTING_MODEL[i]:
			pflux=np.append(pflux,float(gbm.PFLX_SBPL_ERGFLUX[i]))
			pflux_err=np.append(pflux_err,float(gbm.PFLX_SBPL_ERGFLUX_ERROR[i]))
			pflux_alpha=np.append(pflux_alpha,gbm.PFLX_SBPL_INDX1[i])
			pflux_alpha_neg_err=np.append(pflux_alpha_neg_err,gbm.PFLX_SBPL_INDX1_NEG_ERR[i])
			pflux_alpha_pos_err=np.append(pflux_alpha_pos_err,gbm.PFLX_SBPL_INDX1_POS_ERR[i])
			pflux_epeak=np.append(pflux_epeak,gbm.PFLX_SBPL_BRKEN[i])
			pflux_epeak_neg_err=np.append(pflux_epeak_neg_err,gbm.PFLX_SBPL_BRKEN_NEG_ERR[i])
			pflux_epeak_pos_err=np.append(pflux_epeak_pos_err,gbm.PFLX_SBPL_BRKEN_POS_ERR[i])
			pflux_beta=np.append(pflux_beta,gbm.PFLX_SBPL_INDX2[i])
			pflux_beta_neg_err=np.append(pflux_beta_neg_err,gbm.PFLX_SBPL_INDX2_NEG_ERR[i])
			pflux_beta_pos_err=np.append(pflux_beta_pos_err,gbm.PFLX_SBPL_INDX2_POS_ERR[i])

		elif 'PLAW' in gbm.PFLX_BEST_FITTING_MODEL[i]:
			pflux=np.append(pflux,float(gbm.PFLX_PLAW_ERGFLUX[i]))
			pflux_err=np.append(pflux_err,float(gbm.PFLX_PLAW_ERGFLUX_ERROR[i]))
			pflux_alpha=np.append(pflux_alpha,gbm.PFLX_PLAW_INDEX[i])
			pflux_alpha_neg_err=np.append(pflux_alpha_neg_err,gbm.PFLX_PLAW_INDEX_NEG_ERR[i])
			pflux_alpha_pos_err=np.append(pflux_alpha_pos_err,gbm.PFLX_PLAW_INDEX_POS_ERR[i])
			pflux_epeak=np.append(pflux_epeak,0.)
			pflux_epeak_neg_err=np.append(pflux_epeak_neg_err,0.)
			pflux_epeak_pos_err=np.append(pflux_epeak_pos_err,0.)
			pflux_beta=np.append(pflux_beta,0.)
			pflux_beta_neg_err=np.append(pflux_beta_neg_err,0.)
			pflux_beta_pos_err=np.append(pflux_beta_pos_err,0.)

		else:
			pflux=np.append(pflux,0.)
			pflux_err=np.append(pflux_err,0.)
			pflux_alpha=np.append(pflux_alpha,0.)
			pflux_alpha_neg_err=np.append(pflux_alpha_neg_err,0.)
			pflux_alpha_pos_err=np.append(pflux_alpha_pos_err,0.)
			pflux_epeak=np.append(pflux_epeak,0.)
			pflux_epeak_neg_err=np.append(pflux_epeak_neg_err,0.)
			pflux_epeak_pos_err=np.append(pflux_epeak_pos_err,0.)
			pflux_beta=np.append(pflux_beta,0.)
			pflux_beta_neg_err=np.append(pflux_beta_neg_err,0.)
			pflux_beta_pos_err=np.append(pflux_beta_pos_err,0.)

		utc=Time(float(gbm.TRIGGER_TIME[i]),format='mjd',scale='utc')
		tdiff=utc-met0
		tdiff.format='sec'
		met=np.append(met,tdiff.value)
		utc.format='iso'
		trigtime=np.append(trigtime,utc.value)

		i=i+1

	rtable=Table([gbmname,trigtime,met,ra,dec,t90,t90_err,t90_start,fluence,fluence_err,\
		#flux,flux_err,\
		mo,alpha,alpha_neg_err,alpha_pos_err,\
		epeak,epeak_neg_err,epeak_pos_err,beta,beta_neg_err,beta_pos_err,\
		pflux_mo,pflux,pflux_err,pflux_alpha,pflux_alpha_neg_err,pflux_alpha_pos_err,\
		pflux_epeak,pflux_epeak_neg_err,pflux_epeak_pos_err,\
		pflux_beta,pflux_beta_neg_err,pflux_beta_pos_err,\
		pfluxph64,pfluxph256,pfluxph1024],\
		names=['GBMNAME','TRIG_TIME','TRIG_MET','RA','Dec','T90',\
		'T90_err','T90_start','FLUENCE','FLUENCE_err',\
#		'PEAK_FLUX','PEAK_FLUX_err',
		'FLNC_BEST_FITTING_MODEL',\
		'FLNC_ALPHA','FLNC_ALPHA_neg_err','FLNC_ALPHA_pos_err',\
		'FLNC_EPEAK','FLNC_EPEAK_neg_err','FLNC_EPEAK_pos_err',\
		'FLNC_BETA','FLNC_BETA_neg_err','FLNC_BETA_pos_err',\
		'PFLX_BEST_FITTING_MODEL','PFLX','PFLX_err',\
		'PFLX_ALPHA','PFLX_ALPHA_neg_err','PFLX_ALPHA_pos_err',\
		'PFLX_EPEAK','PFLX_EPEAK_neg_err','PFLX_EPEAK_pos_err',\
		'PFLX_BETA','PFLX_BETA_neg_err','PFLX_BETA_pos_err',\
		'PFLX_PH_64','PFLX_PH_256','PFLX_PH_1024'],\
		dtype=('S12','S23','f4','f4','f4','f4','f4','f4','f4','f4','S20','f4','f4','f4','f4','f4','f4','f4','f4','f4','S20','f4','f4','f4','f4','f4','f4','f4','f4','f4','f4','f4','f4','f4','f4'))

#	cat.add_columns(rtable.columns.values())

#	cat.write('LAT_Cat_Table.html',format='ascii.html',overwrite=True)
#	cat.write('LAT_Cat_Table.dat',format='ascii',overwrite=True)

	return rtable

def load_GRBOX(nointernet=False):

	url='http://www.astro.caltech.edu/grbox/grboxtxt.php?form=submitted&starttime=041222&endtime=220101&sort=time&reverse=y&showindex=y&showt90=y&showra=y&showdec=y&showz=y&showut=y&xor=y&ref=y&observatory=t&obsdate=2017-08-18&posfmt=dec&xrtpos=gcn&format=txt'
	filename='grboxtxt.txt'
	if nointernet==False: urllib.urlretrieve(url,filename)
	grbox=ascii.read(filename,format='fixed_width',\
		names=['GRB','UT','T90','RA','DEC','z','det'],data_start=1,\
		col_starts=(0,8,17,23,36,49,55))

	ngrbox=len(grbox['GRB'])

	met0='2001-01-01 00:00:00'
#	gname=[]
#	gfrac=[]
	met=[]
	for i in range(ngrbox):

		if grbox['UT'][i] != '':
			year=grbox['GRB'][i][0:2]
			month=grbox['GRB'][i][2:4]
			day=grbox['GRB'][i][4:6]
			if ':' in grbox['UT'][i][0:2]: q=-1
			else: q=0
			hr=float(grbox['UT'][i][0:2+q])
			mn=float(grbox['UT'][i][3+q:5+q])
			sec=grbox['UT'][i][6+q:8+q]
			if sec != '': sec=float(sec)
			else: sec=0.
			utc='20'+year+'-'+month+'-'+day+' '+grbox[i]['UT']

			times=[met0,utc]
			t=Time(times,format='iso',scale='utc')
			tdiff=t[1]-t[0]
			tdiff.format='sec'
			met=np.append(met,tdiff.value)
		else: met=np.append(met,0.)

	# add redshifts not in GRBOX (from Fong et al. 2015)
	grbs=np.array(['111117A','100625A','100206A','100117A','080905A','070729','070809','071227','090515','101219A'])
	z=np.array([1.3,0.452,0.407,0.915,0.122,0.8,0.473,0.381,0.403,0.718])
	m1,m2=match_catalogs_name(grbox['GRB'],grbs)
	grbox['z'][m1]=z[m2]

#	w=np.where(grbox['z'].mask == True)
#	grbox['z']=0.

	w=np.where(grbox['z'].mask == False)
	grbox['z'][w]=np.array([x.replace('?','') for x in grbox['z'][w]])

	rtable=Table([grbox['GRB'],grbox['UT'],met,grbox['RA'],grbox['DEC'],grbox['z'],grbox['det']],\
		names=['GRB','TRIG_UT','TRIG_MET','RA','Dec','z','afterglow'])#,\
		#		dtype=('S7','S23','f4','f4','f4','f4','S4'))

	return rtable

def load_BAT():

	filename='/Users/jracusin/Swift/BATCAT/batcat_summary_general.txt'
	bat=ascii.read(filename)#,\
#		names=['GRB','trigid','trigmet','trigtime','g_ra','g_dec','f_im','i_snr','t90','t90err',\
#		't50','t50err','evstart','evstop','pcode','comment'],delimiter='|')

	return bat






