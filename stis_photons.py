from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import leastsq,newton
from numpy.lib.recfunctions import append_fields




def localise(fln_tag,fln_x1d,fln_dsp,fln_moc,del_slit = 10,del_line = 40,cen_wav = 1400,del_wav = 15,
            transform = False, verbose=True, plot_diagnostics=True,output='stis_timetag'):
    """Find Perform the best localisation of the slit
    """

    if transform: photons = fits.getdata(fln_tag)
    if not transform: photons = np.load(output)

    opt_elem = fits.getval(fln_tag,'OPT_ELEM')
    spec=fits.getdata(fln_x1d)
    #spec2d=fits.getdata('../spectroscopy/stis/oci802020_x2d.fits')
    wav = fits.getdata(fln_dsp)    ## DSP file
    moc = fits.getdata(fln_moc)	## MOC file
    minwav, maxwav = fits.getval(fln_tag,'MINWAVE'),fits.getval(fln_tag,'MAXWAVE')
    moff1,moff2 = fits.getval(fln_tag,'MOFFSET1'),fits.getval(fln_tag,'MOFFSET2')
    a2center = fits.getval(fln_tag,'CENTERA2')
    cen_wav = fits.getval(fln_tag,'CENWAVE')
    ss = (wav['OPT_ELEM'] == opt_elem) * (wav['A2CENTER'] == a2center)
    tt = moc['OPT_ELEM'] == opt_elem
    a = wav['COEFF'][ss] + moc['COEFF1'][tt]*moff1 + moc['COEFF2'][tt]*moff2
    a = a[0]
    #print(a[0])
    if verbose: print( 'HST-STIS Wavelength assignment routine')

    num_phot = len(photons)
    if verbose: print( 'Number of events in file: '+str(num_phot))

    if plot_diagnostics:
        mask = np.random.randint(num_phot,size=int(.1*num_phot))
        fig=plt.figure(num='Slit Location')
        plt.clf()
        fig.add_subplot(311)
        plt.title('HST-STIS\nFile: '+fln_tag)
        #plt.scatter(photons['AXIS1'][mask],photons['AXIS2'][mask],marker='.',s=5,alpha=0.3)

        xedges = np.arange(0,2048,1)
        yedges = np.arange(0,2048,4.0)

        H, xedges, yedges = np.histogram2d(photons['AXIS1'],photons['AXIS2'], bins=(xedges, yedges))
        #H = np.rot90(H)
        #H = np.flipud(H)
        #Hmasked = np.ma.masked_where(H==0,H)
        #Hmasked = np.ma.masked_where(H==0,H)
        im = plt.imshow(np.arcsinh(H.T), interpolation='nearest', origin='low', aspect='auto',
                extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],cmap='viridis')

        xslit = np.arange(0,1024*2,1.)
        slit = np.interp(xslit,np.arange(0,1024,1)*2.0,spec['EXTRLOCY'][0]*2.0)

        plt.plot(xslit,slit,'r')
        plt.plot(xslit,slit+del_slit,'r--')
        plt.plot(xslit,slit-del_slit,'r--')

        plt.ylim(np.mean(slit)-50,np.mean(slit)+50)
        plt.xlim(0,2024)
        plt.xlabel('Pixel AXIS1')
        plt.ylabel('Pixel AXIS2')
        #stop
        plt.subplots_adjust(top=0.9)


        fig.add_subplot(313)
        plt.plot(spec['WAVELENGTH'][0],spec['FLUX'][0],'k')
        plt.xlim(minwav,maxwav)
        plt.ylim(0,3.5e-15)
        plt.xlabel('Wavelength, $\AA$')
        plt.ylabel('Flux, erg s$^{-2}$ cm$^{-2}$ $\AA^{-1}$')
        plt.axvline(x=cen_wav,color='b',linestyle='--')
        plt.axvline(x=cen_wav+del_wav,color='b',linestyle='--')
        plt.axvline(x=cen_wav-del_wav,color='b',linestyle='--')
        plt.tight_layout()
        plt.subplots_adjust(top=0.95)
        plt.draw()
        plt.show()
    if transform:
    	wave2=[]
    	cc=0
    	for s in photons['AXIS1']/2.0:
    		def disptab(x):
    			return a[0]+a[1]*x + a[2]*x**2+a[3] + a[4]*x +a[5]*x +a[6]*x**2-s
    		wave2.append(newton(disptab,cen_wav))
    		cc+=1
    		Printer('Calculating wavelengths: %5.2f%%'%(float(cc)/float(num_phot)*100))
    	wave2=np.array(wave2)

    	print( '...Done...')
    	def extract(axis1,axis2,xslit,slit):
    		temp = axis2
    		cc=0
    		for i,j in zip(xslit,slit):
    			mm = axis1 == i
    			temp[mm] = axis2[mm] - int(j) + 1000
    			Printer('Correcting Pixel values: %5.2f%%'%(float(cc)/len(slit)*100))
    			cc+=1
    		return temp
    	new_x = extract(photons['AXIS1'],photons['AXIS2'],xslit[:-1],slit[:-1])
    	print ('...Done...')
    fig.add_subplot(312)
    if transform:
    	#plt.scatter(wave2[mask],new_x[mask],marker='.',s=5,alpha=0.3)
        xedges = np.arange(minwav,maxwav,1.0)
        yedges = np.arange(900,1100,2.0)

        H, xedges, yedges = np.histogram2d(wave2,new_x, bins=(xedges, yedges))
        #H = np.rot90(H)
        #H = np.flipud(H)
        #Hmasked = np.ma.masked_where(H==0,H)
        #Hmasked = np.ma.masked_where(H==0,H)
        im = plt.imshow(np.arcsinh(H.T), interpolation='nearest', origin='low', aspect='auto',
                extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],cmap='viridis')
    else:
    	plt.scatter(photons['WAVELENGTH'][mask],photons['AXIS2_X'][mask],marker='.',s=5,alpha=0.3)
    plt.xlim(minwav,maxwav)
    plt.ylim(1000-50,1000+50)
    plt.axvline(x=cen_wav,color='b',linestyle='-')
    plt.axvline(x=cen_wav+del_wav,color='b',linestyle='--')
    plt.axvline(x=cen_wav-del_wav,color='b',linestyle='--')
    plt.axhline(y=1000,color='r',linestyle='-')
    plt.axhline(y=1000+del_slit,color='r',linestyle='--')
    plt.axhline(y=1000-del_slit,color='r',linestyle='--')

    plt.xlabel('Wavelength, $\AA$')
    plt.ylabel('Pixel AXIS2')
    plt.draw()
    fig.add_subplot(311)

    fig.add_subplot(312)
    if transform:
    	mask2=( ( new_x < 1000+del_slit  ) * ( new_x > 1000-del_slit ) ) * ( ( wave2 > cen_wav-del_wav  ) * ( wave2 < cen_wav+del_wav ) )
    	mask_ran=np.random.randint(sum(mask2),size=int(.1*sum(mask2)  ) )
    	plt.scatter(wave2[mask2][mask_ran],new_x[mask2][mask_ran],
                    marker='.',s=5,alpha=0.3,color='y')
    else:
    	mask2=( ( photons['AXIS2_X'] < 1000+del_slit  ) * ( photons['AXIS2_X'] > 1000-del_slit ) ) * ( ( photons['WAVELENGTH'] > cen_wav-del_wav  ) * ( photons['WAVELENGTH'] < cen_wav+del_wav ) )
    	mask_ran=np.random.randint(sum(mask2),size=int(.1*sum(mask2)  ) )
    	plt.scatter(photons['WAVELENGTH'][mask2][mask_ran],
                    photons['AXIS2_X'][mask2][mask_ran], marker='.', s=5,
                    alpha=0.3,color='y')

    fig.add_subplot(311)
    plt.scatter(photons['AXIS1'][mask2][mask_ran],photons['AXIS2'][mask2][mask_ran],
                marker='.',s=5,alpha=0.3,color='g')
    plt.draw()

    plt.tight_layout()
    plt.subplots_adjust(top=0.93)
    if transform:
    	photons=append_fields(photons,'TIME2',photons['TIME'].astype('f8'),usemask=False)
    	photons=append_fields(photons,'WAVELENGTH',wave2,usemask=False)
    	photons=append_fields(photons,'AXIS2_X',new_x,usemask=False)
    	np.save(output+'.npy',photons)
    	np.savetxt(output+'.txt',np.transpose([photons['TIME2'].astype('f8'),
                    photons['WAVELENGTH'].astype('f8'),photons['AXIS2_X']]),
    				delimiter=',',fmt=['%.6f','%.4f','%.0f'])


class Printer():
    """
    Print things to stdout on one line dynamically
    """

    def __init__(self,data):

        sys.stdout.write("\r\x1b[K"+data.__str__())
        sys.stdout.flush()
