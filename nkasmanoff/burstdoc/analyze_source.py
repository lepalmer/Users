#!/usr/bin/env python
"""
------------------------------------------------------------------------
Script to parse values from .sim files and make effective area plots:
Author: Regina Caputo (regina.caputo@nasa.gov)
Date: February 17, 2017
Usage Examples:
import analyze_source
#make get effective area for a list of sim files
analyze_source.getAeff('sim/FarFieldPointSource_*Cos1.0*.sim',triggers=10000,r_sphere=300)
#This returns 3 arrays: energy, aeff and costh

#plot effective area vs. energy with or without GBM data!
analyze_source.plotAeff('sim/FarFieldPointSource_*Cos1.0*.sim',WithGBM=True,save=True)

#plot effective area vs. angle
analyze_source.plotAeffVsAngle('sim/FarFieldPointSource_100.000keV_Cos*.sim',save=True)

------------------------------------------------------------------------
"""

import os
import time
import sys
import fileinput
import numpy
import math
import glob

try:
	import matplotlib.pyplot as plot
	from mpl_toolkits.mplot3d import Axes3D
	import matplotlib.mlab as mlab
	import matplotlib.gridspec as gridspec
	from matplotlib.ticker import MultipleLocator, FormatStrFormatter
	from matplotlib.ticker import AutoMinorLocator
	from matplotlib.colors import LogNorm
	
	# Set the default title font dict
	titleFormat = {'fontsize': 12, 'fontweight' : plot.rcParams['axes.titleweight'], 'verticalalignment': 'baseline', 'horizontalalignment': 'center'}

except:
	print("\n**** Warning: matplotlib not found. Do not try to make plots or bad things will happen! ****")

try:
	import iminuit 
	import probfit
except:
	print("\n**** Warning: iminuit not found. Do not try to fit plots or bad things will happen! ****")


def getDetailsFromFilename(filename):

    '''Function to get the energy and angle from a filename.
    Really should be meta data.'''
    
    details = {}
    info = filename.split('_')
    details['keV'] = info[1][:-3]
    
    angle = info[2].split('.')
    details['Cos'] = "{}.{}".format(angle[0][3:], angle[1])
    
    return details

def getAngle(filename):
    file_text_name = os.path.splitext(os.path.basename(filename))  
    file_last_num = os.path.basename(file_text_name[0]).split('_') 
    deet=file_last_num[2].split('.')
    num = "{}.{}".format(deet[0][3:], deet[1])
    result = float(num)	    
    return result

def getSphere(filename):

    '''Function to get the radius of the surrounding sphere from a geometery file'''
    
    f=open(filename,'r')
    lines=f.readlines()

    for line in lines:
        if 'SurroundingSphere' in line:
            lineContents = line.split()
            print("Radius of Sphere: ", lineContents[1])
            return float(lineContents[1])

def sourceToDict(filename):

	'''Turn source file into a dictionary'''

	geo = {}
	with open(filename) as f:
		for line in f:
			lineContents = line.split()
			if len(lineContents) > 1:
				geo[lineContents[0]] = lineContents[1:]
	return geo

def getKey(filename):
    file_text_name = os.path.splitext(os.path.basename(filename))  
    file_last_num = os.path.basename(file_text_name[0]).split('_') 
    num=file_last_num[1].split('.')
    result = int(num[0])
    return result

def cosFun(x, a, b):
	#takes input angle in degrees
	return a*math.cos(math.radians(x))**b

def cosNorm(x, a):
	#takes input angle in degrees
	return a*math.cos(math.radians(x))

def parse(filename, sourceTheta=None):

    print('\nParsing: %s' % filename)

    #read all the lines in the file
    command = 'wc %s' % filename
    output = os.popen(command).read()
    totalNumberOfLines = int(output.split()[0])

    lineNumber=0
    
    # initialize everything
    energy_deposited = []
    energy_escaped = []
    energy_nonsensitive = []
    

    with open(filename) as f:
	    for line in f:

		    try:
			    sys.stdout.write("Progress: %d%%   \r" % (lineNumber/totalNumberOfLines * 100) )
			    sys.stdout.flush()
		    except:
			    pass
	    
		    if 'TS' in line:
			    
			    # Split the line
			    lineContents = line.split()
			    print("Total generated particles: ", lineContents[1])
			    
		    if 'ED' in line:
			    lineContents = line.split()
			    energy_deposited.append(lineContents[1])

		    if 'EC' in line:
			    lineContents = line.split()
			    energy_escaped.append(lineContents[1])

		    if 'NS' in line:
			    lineContents = line.split()
			    energy_nonsensitive.append(lineContents[1])
		    

		    lineNumber = lineNumber + 1

    #print len(energy_deposited)

    return energy_deposited, energy_escaped, energy_nonsensitive


def passEres(filename, alpha=2.57): #2.57 from 10% at 662 keV scaling at 1/sqrt(E)
	
	ed,es,ns=parse(filename)
	good=0.0
	mod = 0.0
	frac=1.0
	mod_frac=1.0
	escape = 30.0 #escape photon energy in keV for CsI

	for i in range(len(ed)):
		tot = float(ed[i])+float(es[i])+float(ns[i])
		sigma = 0.0
		ediff = tot-float(ed[i])
		ediff2=ediff-escape
		if ed[i] != 0:
			#Note: sigma = FWHM(or eres)/2.35
			sigma = float(ed[i])*alpha/math.sqrt(float(ed[i]))/2.35
		
		if math.fabs(ediff) < sigma:
			good=good+1.0
		if math.fabs(ediff2) < sigma:
			mod=mod+1.0
		
	frac=float(good)/float(len(ed))
	mod_frac=(float(mod)+float(good))/float(len(ed))

	print(filename, "has this fraction of events that pass Eres cut: ", frac, mod_frac)

	return frac, mod_frac


def getGBMData():

    filename='/data/slag2/BurstCube/BurstCube/Simulations/GEANT3/gbm_effective_area.dat'

    f=open(filename,'r')
    lines=f.readlines()
    aeff=[]
    energy=[]

    for i in range(len(lines)):
        data = lines[i].split('\t')
        if i>1:
            energy.append(float(data[0]))
            aeff.append(float(data[1]))
    
    return energy, aeff
    

def getGenPart(filename):
    
    f=open(filename,'r')
    lines=f.readlines()

    #print "reading:", filename, "... has: ",len(lines), " lines"
    last=lines[len(lines)-1]
    hold=last.split(' ')
    gen=hold[1].split('\n')
    #print "The last line is: ",last
    print(int(gen[0]))
    return int(gen[0])
    
        
def CalculateAeff(filename, triggers, r_sphere=None):

    generated_particles = getGenPart(filename)
    sourceDict = sourceToDict(filename)

    if not r_sphere:
        geo_file = sourceDict['Geometry'][0]
        r_sphere = getSphere(geo_file)


    Aeff = r_sphere**2*math.pi*float(triggers)/generated_particles

    #print filename, ", Aeff: ", Aeff
    return Aeff
    

def getAeff(directory, triggers, r_sphere, sortAngle=False):

	if sortAngle:
		filenames = sorted(glob.glob(directory),key=getAngle)
	else:
		filenames = sorted(glob.glob(directory),key=getKey)
		
	full_details={}
        #lists
	energy = []
	aeff = []
	aeff_eres = []
	aeff_eres_modfrac = []
	ang = []
    
	for fn in filenames:
          
                #if '10.000keV' in fn:
		#    continue
                #elif '15.000keV' in fn:
		#    continue

		details = getDetailsFromFilename(fn)
		details['Aeff'] = CalculateAeff(fn,triggers,r_sphere)
		frac, mod_frac=passEres(fn,alpha=2.57) #assumes energy resolution of 10% at 662 keV
		
		
		for key, value in details.iteritems():
                        #print key, value
			if key is 'keV':
				energy.append(value)
			elif key is 'Aeff':
				aeff.append(value)
				aeff_eres.append(float(value)*float(frac))
				aeff_eres_modfrac.append(float(value)*float(mod_frac))
			elif key is 'Cos':
				ang.append(float(value))
			else:
				print("key not found")
				
        #print aeff[4], aeff_eres[4]

	return energy, aeff, ang, aeff_eres, aeff_eres_modfrac

def plotAeff(files, comparison=False, WithGBM=False, save=False, rcirc=None):

    GBM_e=[]
    GBM_aeff=[]

    energy, aeff, ang, aeff_eres, aeff_eres_modfrac=getAeff(files, 10000.,rcirc)
    plot.figure(figsize=(8,6))
    plot.scatter(energy, aeff, color='black')
    plot.plot(energy, aeff, color='black', alpha=0.5, linestyle='--', lw=2, label='BurstCube')

    plot.scatter(energy, aeff_eres, color='blue')
    plot.plot(energy, aeff_eres, color='blue', alpha=0.5, linestyle='--', lw=2, label='BurstCube with E$_{\mathrm{res}}$')

    plot.scatter(energy, aeff_eres_modfrac, color='red')
    plot.plot(energy, aeff_eres_modfrac, color='red', alpha=0.5, linestyle='--', lw=2, label='BurstCube with E$_{\mathrm{res}}$ + escape')

    if comparison: 
	    energy2, aeff2, ang2=getAeff('sim/9.4x9.4cmCube/FarFieldPointSource_*Cos1.0*.sim', 10000.,rcirc)
	    #plot.figure(figsize=(8,6))
	    plot.scatter(energy2, aeff2, color='blue')
	    plot.plot(energy2, aeff2, color='blue', alpha=0.5, linestyle='--', lw=2, label='1 of 4')

	    energy3, aeff3, ang3=getAeff('sim/5.8x5.8cmCube/FarFieldPointSource_*Cos1.0*.sim', 10000.,rcirc)
	    #plot.figure(figsize=(8,6))
	    plot.scatter(energy3, aeff3, color='purple')
	    plot.plot(energy3, aeff3, color='purple', alpha=0.5, linestyle='--', lw=2, label='1 of 9 thick')

    plot.xscale('log')
    plot.xlabel('Energy (keV)', fontsize=16)
    plot.gca().set_xlim([5.,10000.])

    plot.yscale('log')
    plot.gca().set_ylim([1.,200.])
    plot.ylabel('Effective Area (cm$^2$)', fontsize=16)

    legend = plot.legend(loc='lower center',prop={'size':16},numpoints=1,frameon=False)

    if WithGBM:
        print("with GBM!")
        GBM_e, GBM_aeff=getGBMData()
        plot.plot(GBM_e, GBM_aeff, color='green', alpha=0.75, linestyle='-', lw=2, label='GBM NaI')
        
    if save:
        plot.savefig('EffectiveArea_vs_E.png')
        plot.savefig('EffectiveArea_vs_E.pdf')

    plot.show()


def plotAeffVsAngle(files, comparison=False, save=False, doFit= False, rcir=None):
    
    energy, aeff, ang, aeff_eres, aeff_eres_modfrac=getAeff(files, 10000.,rcir,sortAngle=True)

    #print ang
    angle=[]
    for i in range(len(ang)):
	    angle.append(round(numpy.degrees(numpy.arccos(ang[i]))))
    aeff_err=[x * 0.01 for x in aeff]

    if doFit:
	    #print iminuit.describe(cosFun)
	    #print angle
	    chi2 = probfit.Chi2Regression(cosFun, numpy.asarray(angle), numpy.asarray(aeff), numpy.asarray(aeff_err))
	    #print iminuit.describe(chi2)
	    minuit = iminuit.Minuit(chi2, a=80., b=1., error_a=1, error_b=0.01, limit_a=(60.,100.), limit_b=(0.2,1.0))
	    minuit.migrad()
	    print(minuit.values)
	    print(minuit.errors)
	    ((data_edges, datay), err, (total_pdf_x, total_pdf_y), parts) = chi2.draw(minuit);

    if not doFit:
	    chi2Cos = probfit.Chi2Regression(cosNorm, numpy.asarray(angle), numpy.asarray(aeff), numpy.asarray(aeff_err))
	    minuitCos = iminuit.Minuit(chi2Cos, a=80., error_a=1, limit_a=(60.,100.))
	    minuitCos.migrad()
	    print(minuitCos.values)
	    print(minuitCos.errors)
	    ((data_edges_cos, datay_cos), err_cos, (total_pdf_x_cos, total_pdf_y_cos), parts_cos) = chi2Cos.draw(minuitCos);
	    

    plot.figure(figsize=(8,6)) 
    #plot.errorbar(angle, aeff, yerr=aeff_err, color='black',fmt='o',label='BurstCube')
    plot.errorbar(angle, aeff, yerr=aeff_err, color='black',fmt='o')

    if comparison:
	    energy2, aeff2, ang2, aeff_eres2, aeff_eres_modfrac2=getAeff('sim/thin/FarFieldPointSource_100.000keV_Cos*.sim', 10000.,rcir)
	    angle2=[]
	    for i in range(len(ang2)):
		    angle2.append(round(numpy.degrees(numpy.arccos(ang2[i]))))
	    aeff_err2=[x * 0.01 for x in aeff2]

	    #plot.scatter(angle2, aeff2, color='blue', label='thin')
	    #plot.errorbar(angle2, aeff2, yerr=aeff_err2, color='blue',fmt='o',label='thin')
	    plot.errorbar(angle2, aeff2, yerr=aeff_err2, color='blue',fmt='o')

	    if doFit:
		    chi2_thin = probfit.Chi2Regression(cosFun, numpy.asarray(angle2), numpy.asarray(aeff2), numpy.asarray(aeff_err2))
		    minuit_thin = iminuit.Minuit(chi2_thin, a=80., b=1., error_a=1, error_b=0.01, limit_a=(60.,100.), limit_b=(0.2,1.0))
		    minuit_thin.migrad()
		    print(minuit_thin.values)
		    print(minuit_thin.errors)
		    ((data_edges_thin, datay_thin), err_thin, (total_pdf_x_thin, total_pdf_y_thin), parts_thin) = chi2.draw(minuit_thin);

	    
	    energy3, aeff3, ang3, aeff_eres3, aeff_eres_modfrac3=getAeff('sim/thick/FarFieldPointSource_100.000keV_Cos*.sim', 10000.,rcir)
	    angle3=[]
	    for i in range(len(ang3)):
		    angle3.append(round(numpy.degrees(numpy.arccos(ang2[i]))))
	    aeff_err3=[x * 0.01 for x in aeff3]

	    #plot.scatter(angle3, aeff3, color='red', label='thick')
	    #plot.errorbar(angle3, aeff3, yerr=aeff_err3, color='red',fmt='o',label='thick')
	    plot.errorbar(angle3, aeff3, yerr=aeff_err3, color='red',fmt='o')

	    if doFit:
		    chi2_thick = probfit.Chi2Regression(cosFun, numpy.asarray(angle3), numpy.asarray(aeff3), numpy.asarray(aeff_err3))
		    minuit_thick = iminuit.Minuit(chi2_thick, a=80., b=1., error_a=1, error_b=0.01, limit_a=(60.,100.), limit_b=(0.2,1.0))
		    minuit_thick.migrad()
		    print(minuit_thick.values)
		    print(minuit_thick.errors)
		    ((data_edges_thick, datay_thick), err_thick, (total_pdf_x_thick, total_pdf_y_thick), parts_thick) = chi2.draw(minuit_thick);


    if doFit:
	    
	    plot.figure(figsize=(8,6)) 
	    
	    plot.errorbar(angle, aeff, yerr=aeff_err, color='black',fmt='o')

	    function = r"Nominal fit: %.1f * cos($\theta$)$^{%.2f}$" % (minuit.values['a'],minuit.values['b'])
	    plot.plot(total_pdf_x, total_pdf_y, color='black', lw=2, label= function)
	    
	    if comparison: 
		    plot.errorbar(angle2, aeff2, yerr=aeff_err2, color='blue',fmt='o')
		    function2 = r"Thin fit: %.1f * cos($\theta$)$^{%.2f}$" % (minuit_thin.values['a'],minuit_thin.values['b'])
		    plot.plot(total_pdf_x_thin, total_pdf_y_thin, color='blue', lw=2, label= function2)

		    plot.errorbar(angle3, aeff3, yerr=aeff_err3, color='red',fmt='o')
		    function3 = r"Thick fit: %.1f * cos($\theta$)$^{%.2f}$" % (minuit_thick.values['a'],minuit_thick.values['b'])
		    plot.plot(total_pdf_x_thick, total_pdf_y_thick, color='red', lw=2, label= function3)


    plot.gca().set_xlim([0.,80.])
    plot.xlabel('Incident Angle (deg)', fontsize=16)

    plot.gca().set_ylim([20.,100.])
    plot.ylabel('Effective Area (cm$^2$)', fontsize=16)

    legend = plot.legend(loc='lower center',scatterpoints=1,prop={'size':16},frameon=False)

    plot.grid(True)

    if save:
	    plot.savefig('EffectiveArea_vs_Ang.png')
	    plot.savefig('EffectiveArea_vs_Ang.pdf')

    plot.show()
    
