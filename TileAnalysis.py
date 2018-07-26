#!/usr/bin/python
# Aran Garcia-Bellido aran@fnal.gov (July 2018)

import pandas
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize
from scipy import stats
import uncertainties as unc
import uncertainties.unumpy as unp

Data = pandas.read_excel('RochesterAnalysis/Raddam_Data.xlsx',names=['TileName','SampleRun','Size','Dose (Mrad)',' Duration (day)','Comment','DoseRate (rad/h)','TestBeamPE','TestBeamPEError','Bi207Trial1PE','Bi207Trial2PE','Bi207Trial3PE','Bi207AveragePE','Trial1PMTCurr','Trial2PMTCurr','Trial3PMTCurr','PMTCurrAverage (nA)','TBA','TBA2'])

#df=Data['Bi207Trial2PE'].dropna() # Bi207Trial1PE is the name of a column
#print(df)
#s3=Data[Data['Dose (Mrad)']==1.7] # the full table only for sample3 tiles (with 1.7 Mrad)
#print s3['Bi207Trial2PE'].dropna() 

'''
	Plotting Testbeam PE vs Bi207 PE, fit to linear curve
'''
print("Plotting Testbeam PE vs Bi207 PE")
beamtiles=Data[Data['TestBeamPE']>0.]
#beamtiles=beamtiles.sort_values(by=['Bi207AveragePE'])
y=beamtiles['Bi207AveragePE']
yerr=0.05*y
x=beamtiles['TestBeamPE']
xerr=beamtiles['TestBeamPEError']
names=beamtiles['TileName']
print(beamtiles[['TileName','SampleRun','TestBeamPE','TestBeamPEError','Bi207AveragePE']])

plt.errorbar(x, y, yerr, xerr=xerr, linestyle='None', ecolor='k', marker='o', markerfacecolor='black', markersize=5.2, label='Tile')
for i in range(0,len(names)):
	plt.text(x.values[i],y.values[i],names.values[i], fontsize=12, color='g')

#slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
#plt.plot(x,slope*x+intercept,'r',label=r'Fit: $y={0:.2f}x {1:.2f}$'.format(slope,intercept))

def func_pol1(x, a, b):
    return a*x+b

fitres = optimize.curve_fit(func_pol1, x, y, sigma=yerr, p0=[4.,-16.], full_output=True, absolute_sigma=True)
popt=fitres[0]; pcov=fitres[1]
# Need to use correlated errors to draw sigma bands
a, b = unc.correlated_values(popt, pcov)
print('Fit results: \na = {0:.3f} \nb = {1:.3f}'.format(a,b))
redchisq = (fitres[2]['fvec']**2).sum()/(len(fitres[2]['fvec'])-len(popt))
print("chi2/Ndof = %6.3f" % redchisq)
#print len(fitres[2]['fvec'])
finex=np.arange(min(x),max(x),0.2)
plt.ylim([0.,18.])
plt.xlim([0.,55.])
plt.plot(finex, func_pol1(finex, *popt), 'b-',linewidth=2,label=r'Fit: $y={0:.2f}x+{1:.2f}$'.format(*popt)) # fit line
py = func_pol1(finex,a,b)
nom = unp.nominal_values(py)
std = unp.std_devs(py)
#plt.fill_between(finex, nom-1*std, nom+1*std, facecolor='b',alpha=0.3,label='1$\sigma$ band of fit') # bands
#plt.fill_between(finex, nom-2*std, nom+2*std, facecolor='b',alpha=0.2,label='2$\sigma$ band of fit') # bands

plt.legend(loc='upper left', frameon=False, framealpha=1, numpoints=1)

#plt.xlim(0.,10.)
plt.title("")
plt.ylabel(r"Bi207 Source [# PE]")
plt.xlabel(r"Test Beam [# PE]")
plt.grid()
plt.savefig("BeamPE_vs_SourcePE.png")
plt.figure()

''' Light yield: L(d) = L_0 exp(-d/D) = L_0 exp(-d*mu)
d=dose (Mrad); L_0 = initial light yield; D=constant (Mrad); mu=1/D
D can be obtained from weighted mean of mu in each R (doserate) bin.
We can fit: D=a*R^b
But some plots show saturation at high R, so we can use: D=(a*R^b)/(1+c*R^b)
In Fig. 13 of note DN-18-001, you can see: D=sqrt(R)/(a+b*sqrt(R))
https://indico.cern.ch/event/726619/contributions/2990245/attachments/1643942/2626684/DN-18-001_temp.pdf
'''

def func_nexp(x,a,b):
	return a*np.exp(-x/b)
	
# plot light yield L vs dose d to get constant D
def plot_lightyield_vs_dose(d,title=''):
	print(" Calculating D for " + title)
	plt.figure()
	y=d['Bi207AveragePE']
	yerr=0.7 # for some reason if I use 0.05*y, the plt.errorbar crashes with the contains N2 condition above...?
	x=d['Dose (Mrad)']
	names=d['TileName']
	print(d[['TileName','SampleRun','Dose (Mrad)','Comment','DoseRate (rad/h)','Bi207AveragePE']])

	#plt.plot(x, y, linewidth=2.0, linestyle='-', color='black', marker='o', markerfacecolor='r', markersize=5.2)
	plt.errorbar(x, y, yerr=yerr, linestyle='None', ecolor='k', marker='o', markerfacecolor='black', markersize=5.2, label='Tile')
	for i in range(0,len(names)):
		plt.text(x.values[i],y.values[i],names.values[i], fontsize=12, color='g')

	fitres = optimize.curve_fit(func_nexp, x, y, sigma=yerr, full_output=True, absolute_sigma=True)
	popt=fitres[0]; pcov=fitres[1]
	# Need to use correlated errors to draw sigma bands
	a, b = unc.correlated_values(popt, pcov)
	print('Fit results: \na = {0:.3f} \nb = {1:.3f}'.format(a,b))
	redchisq = (fitres[2]['fvec']**2).sum()/(len(fitres[2]['fvec'])-len(popt))
	print("chi2/Ndof = %6.3f" % redchisq)
	#print len(fitres[2]['fvec'])
	finex=np.arange(min(x),max(x),0.2)
	plt.plot(finex, func_nexp(finex, *popt), 'b-',linewidth=2,label=r'Fit: $y={0:.2f}e^{{-x/{1:.2f}}}$'.format(*popt)) # fit line
	plt.legend(loc='upper right', frameon=False, framealpha=1, numpoints=1)
	plt.xscale('log') # make it log scale!
	plt.xlim([5.e-2,20.])
	plt.title(title)
	plt.ylabel(r"Bi207 Source [# PE]")
	plt.xlabel(r"Dose d [Mrad]")
	plt.grid()
	plt.savefig("SourcePELightYield_vs_Dose"+title+".png")
	return a,b

Data['Comment'].to_string()
nocomment=Data[(Data['Bi207AveragePE']>0.)&(Data['Comment']==' ')]
nitrogen =Data[(Data['Bi207AveragePE']>0.)&(Data['Comment'].str.contains('N2')==True)]
oxygen   =Data[(Data['Bi207AveragePE']>0.)&(Data['Comment'].str.contains('O2')==True)]

a, normalD = plot_lightyield_vs_dose(nocomment,'Normal tiles')
a, nitroD  = plot_lightyield_vs_dose(nitrogen,'Nitrogen tiles')
a, oxyD    = plot_lightyield_vs_dose(oxygen,'Oxygen tiles')

# plot Dose Constant D (Mrad) vs Dose Rate (krad/h) 



