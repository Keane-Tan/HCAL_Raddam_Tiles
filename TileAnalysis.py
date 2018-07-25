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

s3=Data[Data['Dose (Mrad)']==1.7] # the full table only for sample3 tiles (with 1.7 Mrad)
#print s3['Bi207Trial2PE'].dropna() 

beamtiles=Data[Data['TestBeamPE']>0.]
#beamtiles=beamtiles.sort_values(by=['Bi207AveragePE'])
x=beamtiles['Bi207AveragePE']
print x
xerr=0.05*x
y=beamtiles['TestBeamPE']
print y
yerr=beamtiles['TestBeamPEError']
names=beamtiles['TileName']
print(names)

plt.errorbar(x, y, yerr, xerr=xerr, linestyle='None', ecolor='k', marker='o', markerfacecolor='black', markersize=5.2, label='Tile')
for i in range(0,len(names)):
	plt.text(x.values[i],y.values[i],names.values[i], fontsize=12, color='g')

#slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
#plt.plot(x,slope*x+intercept,'r',label=r'Fit: $y={0:.2f}x {1:.2f}$'.format(slope,intercept))

def func_pol1(x, a, b):
    return a*x+b

fitres = optimize.curve_fit(func_pol1, x, y, sigma=yerr*10, p0=[4.,-16.], full_output=True, absolute_sigma=True)
popt=fitres[0]; pcov=fitres[1]
# Need to use correlated errors to draw sigma bands
a, b = unc.correlated_values(popt, pcov)
print('Fit results: \na = {0:.3f} \nb = {1:.3f}'.format(a,b))
redchisq = (fitres[2]['fvec']**2).sum()/(len(fitres[2]['fvec'])-len(popt))
print("chi2/Ndof = %6.3f" % redchisq)
#print len(fitres[2]['fvec'])
finex=np.arange(4.,16.,0.2)
plt.ylim([0.,55.])
plt.plot(finex, func_pol1(finex, *popt), 'b-',linewidth=2,label=r'Fit: $y={0:.2f}x {1:.2f}$'.format(*popt)) # fit line
py = func_pol1(finex,a,b)
nom = unp.nominal_values(py)
std = unp.std_devs(py)
plt.fill_between(finex, nom-1*std, nom+1*std, facecolor='b',alpha=0.3,label='1$\sigma$ band of fit') # bands
plt.fill_between(finex, nom-2*std, nom+2*std, facecolor='b',alpha=0.2,label='2$\sigma$ band of fit') # bands

plt.legend(loc='upper left', frameon=False, framealpha=1, numpoints=1)

#plt.xlim(0.,10.)
plt.title("")
plt.xlabel(r"Bi207 Source [# PE]")
plt.ylabel(r"Test Beam [# PE]")
plt.grid()
plt.savefig("BeamPE_vs_SourcePE.png")
