#!/usr/bin/python
# Aran Garcia-Bellido aran@fnal.gov (May 2018)

# To run on Windows, you can open the cmd.exe app and type:
# cd radtest
# python TileCharacterization.py 53_2_May23.Spe -o
# Or you can select those commands here and copy (ctrl+C) and then right click on the cmd.exe terminal

# Things we could improve: 
# 1) Make a big text file with each tile name and the most optimal values for -d -n -m -z in a row.
#    That way once we have collected all measurements, we can read that text file on the fly and we 
#    don't have to find and enter the best values every time... (work in progress: need to make the program read the best values)
# 2) Make two separate command line arguments for the gain calculation: the peak to start from and 
#    how many peaks to count. Something like -gn 1 and -gm 1, for example, would only calculate the
#    gain from the peak#2-peak#1. And then the avgPE would still be calculated after -m 8
# 3) The whole passing of global variables is pretty awful. We could create a class called TileParameters
#    that stores all the info that is needed in the functions and contains the current values of the 
#    parameters, specially mpd when we are optimizing...
# 4) Write down whenever the -o option doesn't work on a file, and try to improve the optimize_peak_finding
#    function so that it would work next time.
# 5) Add an uncertainty to the avgPE calculation! Error propagation on the formula? 
# 6) For the extremely low yield tiles: try subtracting the pedestal??? Normalization could be a big issue! 

import sys, os, argparse, linecache
import numpy as np
import matplotlib.pyplot as plt
from detect_peaks import detect_peaks, smooth_spectrum

_MinPeakNum=8 # The number of PE peaks to cut (and the first peak to begin counting the gain)
NPeaks=8      # The range of PE peaks to count the gain. 
# The next two parameters set the range for the x-axis ZOOM plot of the peaks (where we want 20 peaks or so to show up)
_MinZoomADC=0      # Sometimes this may be higher, depending on the full spectrum
_MaxZoomADC=900    # This is our estimate of where the first 20 peaks or so will appear in the x-range
_MinPeakADCDist=20 # Used to find peaks and valleys: minimum distance between peaks
basename = '' # Input file's name without extension; used for writing the data in OptimalTilePara.txt
avgPE = 0

def DecodeArguments():
	global _MinPeakADCDist, _MinZoomADC, _MaxZoomADC, _MinPeakNum, NPeaks
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=''' 
python TileCharacterization.py "C:\User\test 13 may 8.Spe" 
python TileCharacterization.py test_14_may_8.Spe -d 30 -z 800 -s
To plot the pedestal run: python TileCharacterization.py May_22_Sample_84_1.Spe May_22_Sample_84_pedestal.Spe -d 40 -z 1000
''')
	parser.add_argument('f', metavar='N', type=str, nargs='+',help='Input Spe file(s)')
	# We could input here many files, like: test*Spe and loop over them...
	parser.add_argument('-d', nargs=1, help='Set the minimum distance between peaks, in ADC counts')
	parser.add_argument('-z', nargs=1, help='Set the maximum zoom range, in ADC counts. Where we expect the first 20 peaks to be contained.')
	parser.add_argument('-m', nargs=1, help='Peak number to start counting after the noise peaks')
	parser.add_argument('-n', nargs=1, help='Number of peaks to calculate the range from')
	parser.add_argument('-s', action='store_true', help='Do smoothing of the spectrum', default=False)
	parser.add_argument('-o', action='store_true', help='Optimize parameters automatically. Try to start from low mpd (-d 20). It will not find correct values if mpd is larger than gain.', default=False)
	
	if len(sys.argv)==1:
	    parser.print_help()
	    sys.exit(1)
	args = parser.parse_args()
	if args.d:
		_MinPeakADCDist=int(args.d[0])
	if args.z:
		_MaxZoomADC=int(args.z[0])
	if args.m:
		_MinPeakNum=int(args.m[0])
	if args.n:
		NPeaks=int(args.n[0])
	return args
	
def optimize_peak_finding(y):
	''' Make sure we get a proper order of peaks and valleys. Determine the zoom size where the first 20 peaks or will show up, and the proper minimum ADC distance between peaks (mpd for detect_peaks). 
	Could pass the gain and make sure the mpd is not bigger than the gain, or calculate the gain on the first two peaks...'''
	# To do: if there is no peak in ~2*mpd, we should reduce the mpd
	# Look at smoothing https://pythonhosted.org/scikits.datasmooth/regularsmooth.html
	global _MinPeakADCDist, _MinZoomADC, _MaxZoomADC, _MinPeakNum, NPeaks
	NonzeroIndeces=y.nonzero(); xmin=np.amin(NonzeroIndeces); xmax=np.amax(NonzeroIndeces); width=xmax-xmin
	xminzoom=xmin-0.02*width; xmaxzoom=xmax
	p=detect_peaks(y,edge='falling',mpd=_MinPeakADCDist)
	v=detect_peaks(y,edge='falling',mpd=_MinPeakADCDist,valley=True)
	Np=p.size; Nv=v.size
	if (p.size>10 and v.size>10): print p[:10], v[:10];
	# Check whether there are two consecutive peaks before any valley, then increase mpd by 2:
        NTrial = 0 
	while True:
		NBad=0
		for i in np.arange(0,min(p.size,v.size,20)-1): # loop up to 20 peaks or smallest number
			if p[0]>v[0] and p[i] < v[i]:
				print("Wrong! Peak is before valley! i=%i p[i]=%i v[i]=%i" % (i, p[i], v[i]))
				NBad+=1
			if p[i]>v[i+1]:
				NBad+=1
				print("Wrong! Two consecutive peaks with no valleys in between!")
			if y[p[i]]<y[p[i+1]]:
                                print("Wrong! This peak height is smaller than next peak height! i=%i y[p[i]]=%i y[p[i+1]]=%i"%(i,y[p[i]],y[p[i+1]]))
                                NBad+=1
		if NBad>1 and NTrial<10: # One error is ok, specially in the hump region. If two errors, then rerun:
			_MinPeakADCDist+=2
			print("Optimize: Going again NBad=%i, mpd=%i"%(NBad, _MinPeakADCDist))
                        NTrial+=1
			return optimize_peak_finding(y)
		else:
			break
	
	mpd=_MinPeakADCDist
	if (p.size>20):
		xmaxzoom=p[20]+1 # the location of the 20th peak
	else:
		xmaxzoom=p[p.size-1]+1 # or as high as we can go
	print('Optimized parameters: _MinZoomADC=%i _MaxZoomADC=%i _MinPeakADCDist=%i NPeaks=%i' % (xminzoom, xmaxzoom, mpd,NPeaks))
	return xminzoom, xmaxzoom, mpd

def calculate_gain(counts,bins,MinGainPNum=8,NGainPeaks=8):
	''' The gain is defined as the distance in bins (=ADC counts) between two peaks. 
	Here we use the average distance between the 8th peak and 18th (8+NGainPeaks) and divide by NGainPeaks. 
	Remember the first peak, peaks[0], is always the pedestal.
	The decision of where to start counting, depends on running a dark-count spectrum and seeing 
	what is the maximum number of phe peaks you can get. In our case, the dark-count spectrum produces up to 6 or 7 discernible phe peaks.'''
	g=0.
	peaks=detect_peaks(counts,edge='falling',mpd=_MinPeakADCDist) # make sure this is called with the same parameters as outside this function!
	g=float(bins[peaks[MinGainPNum+NGainPeaks]]-bins[peaks[MinGainPNum]])/float(NGainPeaks)
	print("Calculated photo-electron gain using peaks %2ith to %2ith = %4.2f" % (MinGainPNum,MinGainPNum+NGainPeaks,g))
	return g

def plot_full_spectrum(y,basename, yped=np.array([0])):
	''' Find max and min ranges and plot the complete spectrum. 
	Will need to zoom in for interesting regions later. If you pass the pedestal counts, it will overlay that spectrum. It will find all the peaks in this range too.'''
	global _MinPeakADCDist
	plt.plot(y, 'b', lw=1, label=basename+'.Spe') # Draws the spectrum
	plt.yscale('log')
	NonzeroIndeces=y.nonzero(); xmin=np.amin(NonzeroIndeces); xmax=np.amax(NonzeroIndeces)
	width=xmax-xmin
	ymax=np.amax(y)
	plt.axis([xmin-0.02*width, xmax, 1, ymax+0.5*(ymax-1.)])
	plt.xlabel('Charge: ADC channel #', fontsize=14)
	plt.ylabel('Counts', fontsize=14)
	p=detect_peaks(y[:xmax],edge='falling',mpd=_MinPeakADCDist)
	if p.size:
            peaknum=np.arange(p.size)
            #print(peaknum)
            plt.plot(p, y[p], '+', mfc=None, mec='r', mew=2, ms=8,
                    label='%d Peaks (mpd=%i)' % (p.size,_MinPeakADCDist))
            for i in peaknum: # write the peak number besides the cross
                plt.text(p[i],y[p[i]],peaknum[i])
	#plt.title(basename+'.Spe')
	if len(yped) > 2:
		plt.plot(yped, 'gray', lw=1, linestyle='-', label='Pedestal Run') # Draws the spectrum
	plt.legend(loc='best', framealpha=.5, numpoints=1)
	#plt.show()
	plt.savefig("figs/"+basename+"_full.png")
	plt.clf(); plt.cla()

def plot_cuts(y,p,v,basename,description='',mpd=_MinPeakADCDist,xmin=_MinZoomADC,xmax=_MaxZoomADC,ymin=1,ymax=2.e6,cut=1,avgpe=0.,pegain=0.):
	''' Adapted from _plot() in detect_peaks module 
	y are the counts
	p are the indeces for the peaks
	v are the indeces for the valleys'''
	plt.plot(y, 'b', lw=1) # Draws the spectrum
        if p.size:
            peaknum=np.arange(p.size)
            plt.plot(p, y[p], '+', mfc=None, mec='r', mew=2, ms=8,
                    label='%d %s mpd=%i' % (p.size, 'Peaks',mpd))
            for i in peaknum: # write the peak number besides the cross
                plt.text(p[i],y[p[i]],peaknum[i])
            
        if v.size:
            peaknum=np.arange(v.size)
            plt.plot(v, y[v], '+', mfc=None, mec='g', mew=2, ms=8,
                    label='%d %s mpd=%i' % (v.size, 'Valleys',mpd))
        plt.legend(loc='best', framealpha=.5, numpoints=1,title=basename+'.Spe')
        plt.xlabel('Charge: ADC channel #', fontsize=14)
        plt.ylabel('Counts', fontsize=14)
        plt.title(description)
        plt.yscale('log')  # plt.grid()
	plt.axis([xmin, xmax, ymin, ymax])
        x=np.arange(xmin,xmax)
        plt.fill_between(x[cut:],y[x[cut:]],facecolor='b',alpha=0.3)
	plt.text(x[cut],ymin+5,' Avg = %4.2f pe/MIP' % avgpe,fontsize=18,color='b')
	# draw lines for gain: 
	xleft=x[p[_MinPeakNum]]; xright=p[_MinPeakNum+NPeaks]; xcenter=xleft+(xright-xleft)/2
	ylo=y[p[_MinPeakNum]]; yhi=10000 # ymax-0.5*ymax
	plt.plot([xleft,xleft], [ylo,yhi], 'k--', lw=1)
	plt.plot([xright,xright], [ylo,yhi], 'k--', lw=1)
	plt.annotate('', xy=[xleft,yhi], xytext=[xright,yhi], arrowprops={'arrowstyle': '<->'}, ha='center')
	plt.annotate('Gain=%4.2f'%pegain, xy=[xleft,yhi], xytext=[xcenter,yhi], ha='center',va='top')
	plt.savefig("figs/"+basename+"_cuts.png")
	plt.show()
	plt.clf(); plt.cla()

def calculate_avgPE_and_PE(y,PEgain,bins,label,xindex=0,qpercentile=0.99):
	# Calculate last 1% final counts:
	totPE=(bins*y)/PEgain
	cstotPE=np.cumsum(totPE)  # cumulative sum
	cstotPE_LastQPercentile=np.where(cstotPE>qpercentile*np.sum(totPE)) # a subarray with the indeces of cstotPE where it is bigger than the last 1% percentile (Qperc=0.99) of totPE
	Qindex=cstotPE_LastQPercentile[0][0]       # The first index that is bigger than Q percentile
	NonzeroIndeces=totPE.nonzero(); Xmax=np.amax(NonzeroIndeces)
	print(label+'Top %i percentile found after bin %i of non-zero max of: %i'%((1.-qpercentile)*100,Qindex,Xmax))
	print(label+'PE counts of the last %i percentile: %i. [ %i/%i=%5.4f ]' % ( (1.-qpercentile)*100, np.sum(totPE[Qindex:]), np.sum(totPE[Qindex:]),np.sum(totPE[0:]),np.sum(totPE[Qindex:])/np.sum(totPE[0:]) ) )
        #print('Avg Photoelectrons per event from %i bin= %4.2f' % (Qindex,np.sum(PE_fromCut[Qindex:])))
        print(label+'%i Percentile channel in PE =%4.2f'%((1.-qpercentile)*100,Qindex/PEgain))

def saveTileParam(): # this function saves the optimal tile parameters
	f = open('OptimalTilePara.txt','r+')
	f1 = f.readlines()

	if any(basename in s for s in f1):
    		print('Optimal tile parameters already existed in the file.')
	else:
	    f1.append(basename+'\t%i\t%i\t%i\t%i\t%4.2f\n' % (_MinPeakADCDist,_MaxZoomADC,_MinPeakNum,NPeaks,avgPE))
	    firstRow = f1[0] 
	    f1.remove(firstRow) # get rid of the column titles of the table
	    f1.sort() 		# sort the entries in the file
	    f1.insert(0,firstRow) # put the column titles of the table back
	    print("Optimal tile parameters are saved.")
	    f.seek(0) 		# this together with .truncate() allows us to overwrite the text file
	    f.writelines(f1)
	    f.truncate()

	f.close()

def main():
	global _MinPeakADCDist, _MinZoomADC, _MaxZoomADC, _MinPeakNum, NPeaks, basename, avgPE
	args=DecodeArguments()
	inputfilename=args.f[0] # inputfilename="C:/User/test 13 may 8.Spe"
	# Read txt file, skipping the first 12 rows (header) and the last 15 rows (footer)
	# Could also make it read only 8192 entries starting after the row: 0 8191
	data = np.genfromtxt(inputfilename,delimiter=" ",names=['counts'],skip_header=12,skip_footer=15)
	counts=data['counts']
	if args.s:
		counts=smooth_spectrum(counts)
	descriptionline=linecache.getline(inputfilename, 2).strip() # just read one line in textfile
	# Try to find pedestal spectrum in the second argument of the files
	yped=np.array([0])
	if len(args.f)>1:
		pedfile=args.f[1]
		ped_data = np.genfromtxt(pedfile,delimiter=" ",names=['counts'],skip_header=12,skip_footer=15)
		yped=ped_data['counts']
		print(' Will use pedestal file %s with %i rows of data' % (pedfile,len(yped)))
	bins=np.arange(1,len(counts)+1) # start at 1 and go to 8192
	cb=counts*bins
	print(" Read %i rows of data counts from file %s" % (len(counts),inputfilename))
	print(" Tile: %s" % descriptionline)
	print(" Initial values: _MinPeakADCDist=%i _MinZoomADC=%i _MaxZoomADC=%i _MinPeakNum=%i NPeaks=%i" % (_MinPeakADCDist, _MinZoomADC, _MaxZoomADC,_MinPeakNum, NPeaks))
	dirname, fname = os.path.split(inputfilename) # extract path (directories) and filename
	basename, ext = os.path.splitext(fname) # extract extension from filename
	
	# Try to optimize the main parameters:
	if args.o:
		_MinZoomADC, _MaxZoomADC, _MinPeakADCDist = optimize_peak_finding(counts)
	
	plot_full_spectrum(counts,basename,yped)
	#print data['counts'][10]
	#print(data['row'], data['x'],data['xerr']) 
	peaks=detect_peaks(counts[:_MaxZoomADC],edge='falling',mpd=_MinPeakADCDist) # peaks is the index i of the peak in counts[i]
	# Sanity check: 
	if peaks.size < _MinPeakNum:
		print('ERROR: # of peaks found = %i up to _MaxZoomADC is smaller than parameter _MinPeakNum=%i'%(peaks.size,_MinPeakNum))
		sys.exit(2)
	if peaks.size <= _MinPeakNum+NPeaks:
		print('EROR: # of peaks found = %i up to _MaxZoomADC is <= than _MinPeakNum+NPeaks+1=%i needed to calculate the gain.'%(peaks.size,_MinPeakNum+NPeaks))
		print('Either increase _MaxZoomADC, reduce _MinPeakNum+NPeaks, or reduce _MinPeakADCDist')
		sys.exit(3)
	#print peaks.size, _MaxZoomADC, peaks[_MinPeakNum], peaks[_MinPeakNum+NPeaks]
	#if _MaxZoomADC < peaks[_MinPeakNum+NPeaks]:
	#	_MaxZoomADC *= 2*_MaxZoomADC
	#	print('Changing _MaxZoomADC to %i' % _MaxZoomADC)
	
	#print("Peak#  Idx  ADC")
	#for i in np.arange(len(peaks)):
	#	print("%5i %4i %4i" %(i,peaks[i],bins[peaks[i]]))
	
	PEgain=calculate_gain(counts,bins,_MinPeakNum,NPeaks)
	
	valleys=detect_peaks(counts[:_MaxZoomADC],edge='falling',mpd=_MinPeakADCDist,valley=True)
	#print("Min#  Idx  ADC")
	#for i in np.arange(len(valleys)):
	#	print("%5i %4i %4i" %(i,valleys[i],bins[valleys[i]]))
	print 'Peak positions:', peaks
	print 'Valley positions:',valleys
	# Calculate average photoelectrons per MIP: We start counting on the valley after the _MinPeakNum.
	ValleyCut=valleys[_MinPeakNum+1] # this gives errors when valleys are not calculated correctly. 
	# We can just start counting at half the distance after the _MinPeakNum:
	Cut=peaks[_MinPeakNum]+int((peaks[_MinPeakNum+1]-peaks[_MinPeakNum])/2)
	#if np.abs(Cut-ValleyCut)>int(_MinPeakADCDist/2):
	#	print('Double check the position of the Cut: ')
	print('Using cut at index=%3i counts=%4i bin=%4i'%(Cut,counts[Cut],bins[Cut]))
	# Calculate photo electron average:
	PE=(bins*counts)/PEgain/(np.sum(counts[Cut:])+1)
	avgPE=np.sum(PE[Cut:])
	print('Avg Photoelectrons per event = %4.2f' % avgPE)
	# Calculate error on the mean:
	#sigma=
	plot_cuts(counts[:_MaxZoomADC],peaks,valleys,basename,descriptionline,mpd=_MinPeakADCDist,xmax=_MaxZoomADC,cut=Cut,avgpe=avgPE,pegain=PEgain)
	
	# Calculate total count from 0:
	Cut=0
	totPE=(bins*counts)/PEgain
	totPE_fromCut=np.sum(totPE[Cut:])
	print('Total PE counts from %i bin = %i'%(Cut,totPE_fromCut))
        PE_fromCut=(bins*counts)/PEgain/(np.sum(counts[Cut:])+1)
        avgPE_fromCut=np.sum(PE_fromCut[Cut:])
        print('Avg Photoelectrons per event from %i bin= %4.2f' % (Cut,avgPE_fromCut))
	print('The location of the first yped peak is: %i' % (np.argmax(yped)))
	print('The location of the first counts peak is: %i' %  (np.argmax(counts)))

	# Calculate last 1% final counts:
	calculate_avgPE_and_PE(counts,PEgain,bins, '')

        # Subtract pedestal from counts:
        PedSubCounts=counts-yped
        #print PedSubCounts[0:40]
	plot_full_spectrum(PedSubCounts,basename+'_pedsub_')
	calculate_avgPE_and_PE(PedSubCounts,PEgain,bins, 'PedSub')
	print('The location of the first pedSub trough is: %i' %  (np.argmin(PedSubCounts)))        
	## Smoothing really changes the peak location ability
	#logfile = open('smooth.Spe', 'w')
	#sy=smooth_spectrum(counts)
	#for i in np.arange(0,counts.size):
		#logfile.write('%i %i \n'%(sy[i],counts[i]))
	#logfile.close()
	#plot_full_spectrum(sy[:_MaxZoomADC],'caca')
	
	saveTileParam()
main()
