from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import leastsq,newton
from os import system
import sys

def localise(fln_tag, fln_x1d, fln_dsp, fln_moc, del_slit = 10, del_line = 40,
             cen_wav=1400, del_wav = 15, verbose=True, plot_diagnostics=True,
             output='stis_timetag',clobber=True):
    """Find Perform the best localisation of the slit
    """

    data = fits.open(fln_tag)[1].data
    head = fits.open(fln_tag)[0].header

    opt_elem = fits.getval(fln_tag,'OPT_ELEM')
    spec=fits.getdata(fln_x1d)
    wav = fits.getdata(fln_dsp)    ## DSP file
    moc = fits.getdata(fln_moc)	## MOC file
    minwav = fits.getval(fln_tag,'MINWAVE')
    maxwav = fits.getval(fln_tag,'MAXWAVE')
    moff1 = fits.getval(fln_tag,'MOFFSET1')
    moff2 = fits.getval(fln_tag,'MOFFSET2')
    a2center = fits.getval(fln_tag,'CENTERA2')
    cen_wav = fits.getval(fln_tag,'CENWAVE')
    ss = (wav['OPT_ELEM'] == opt_elem) * (wav['A2CENTER'] == a2center)
    tt = moc['OPT_ELEM'] == opt_elem
    a = wav['COEFF'][ss] + moc['COEFF1'][tt]*moff1 + moc['COEFF2'][tt]*moff2
    a = a[0]

    xslit = np.arange(0,1024*2,1.)
    slit = np.interp(xslit,np.arange(0,1024,1)*2.0,spec['EXTRLOCY'][0]*2.0)

    if verbose: print( 'HST-STIS Wavelength assignment routine')

    num_phot = len(data)
    if verbose: print( 'Number of events in file: '+str(num_phot))
    wave2 = np.zeros(num_phot)
    cc=0
    for i,s in enumerate(data['AXIS1']/2.0):
        def disptab(x):
            return a[0]+a[1]*x + a[2]*x**2+a[3] + a[4]*x +a[5]*x +a[6]*x**2-s
        wave2[i] = newton(disptab,cen_wav)
        cc+=1
        if verbose:
            Printer('Calculating wavelengths: %5.2f%%'%(float(cc)/float(num_phot)*100))

    if verbose: print( '...Done...')
    def extract(axis1,axis2,xslit,slit):
        temp = np.zeros(axis2.size)
        cc=0
        for i,j in zip(xslit,slit):
            mm = axis1 == i
            temp[mm] = axis2[mm] - int(j) + 1000
            if verbose:
                Printer('Correcting Pixel values: %5.2f%%'%(float(cc)/len(slit)*100))
            cc+=1
        return temp
    new_x = extract(data['AXIS1'],data['AXIS2'],xslit[:-1],slit[:-1])
    if verbose: print('...Done...')

    ## SAVE NEW FITS FILE
    #if clobber: system('rm '+output+'.fits')

    cols = []
    cols.append(fits.Column(name='WAVELENGTH', format='1E',unit = 'angstrom',
                bscale = 1.0, bzero=0.0, array=wave2 ))
    cols.append(fits.Column(name='AXIS2_X', format='1E',unit = 'pixels',
                bscale = 1.0, bzero=0.0, array=new_x ))
    orig_cols = data.columns
    new_cols = fits.ColDefs(cols)
    c = orig_cols + new_cols

    head['HISTORY'] = 'Assigned wavelengths with stis_photons'
    #head['history'] = 'http://alymantara.github.io'

    #hdu = fits.BinTableHDU.from_columns(c, header=head)
    #hdu.writeto()

    new_hdul = fits.HDUList()
    new_hdul.append(fits.PrimaryHDU(header=head))
    new_hdul.append(fits.BinTableHDU.from_columns(c, header=fits.open(fln_tag)[1].header, name='DATA'))
    new_hdul.append(fits.BinTableHDU(fits.open(fln_tag)[2].data, header=fits.open(fln_tag)[2].header, name='GTI'))
    new_hdul.writeto(output+'.fits', clobber=clobber)

    if plot_diagnostics: plotter(output+'.fits',fln_x1d,fln_dsp,fln_moc,
                del_slit = 10,del_line = 40,cen_wav = 1400,del_wav = 15)

def plotter(fln_tag,fln_x1d,fln_dsp,fln_moc,del_slit = 10,del_line = 40,
            cen_wav = 1400,del_wav = 15):
    """Documentation coming soon
    """
    try:
        vlo = fits.getval(fln_tag,'TTYPE5',1)
    except:
        print("Please run localise first to plot diagnostics")
        return
    data = fits.getdata(fln_tag)
    data_ori = fits.getdata(fln_x1d[:-8]+'tag.fits')
    opt_elem = fits.getval(fln_tag,'OPT_ELEM',1)
    spec=fits.getdata(fln_x1d)

    wav = fits.getdata(fln_dsp)    ## DSP file
    moc = fits.getdata(fln_moc)	## MOC file
    minwav = fits.getval(fln_tag,'MINWAVE')
    maxwav = fits.getval(fln_tag,'MAXWAVE')
    moff1 = fits.getval(fln_tag,'MOFFSET1')
    moff2 = fits.getval(fln_tag,'MOFFSET2')
    a2center = fits.getval(fln_tag,'CENTERA2',1)
    cen_wav = fits.getval(fln_tag,'CENWAVE',1)
    ss = (wav['OPT_ELEM'] == opt_elem) * (wav['A2CENTER'] == a2center)
    tt = moc['OPT_ELEM'] == opt_elem
    a = wav['COEFF'][ss] + moc['COEFF1'][tt] * moff1 + moc['COEFF2'][tt] * moff2
    a = a[0]

    num_phot = len(data)
    mask = np.random.randint(num_phot,size=int(.1*num_phot))

    ##### DO THE PLOTTING
    fig=plt.figure(num='Slit Location',figsize=(6,8))
    plt.clf()
    fig.add_subplot(311)
    plt.title('HST-STIS\nFile: '+ fln_tag)

    xedges = np.arange(0,2048,1)
    yedges = np.arange(0,2048,4.0)

    H, xedges, yedges = np.histogram2d(data_ori['AXIS1'],data_ori['AXIS2'],
                        bins=(xedges, yedges))

    im = plt.imshow(np.arcsinh(H.T), interpolation='nearest', origin='low',
                    aspect='auto', extent=[xedges[0], xedges[-1], yedges[0],
                    yedges[-1]],cmap='binary')

    xslit = np.arange(0,1024*2,1.)
    slit = np.interp(xslit,np.arange(0,1024,1)*2.0,spec['EXTRLOCY'][0]*2.0)

    plt.plot(xslit,slit,'r')
    plt.plot(xslit,slit+del_slit,'r--')
    plt.plot(xslit,slit-del_slit,'r--')

    plt.ylim(np.mean(slit)-50,np.mean(slit)+50)
    plt.xlim(0,2024)
    plt.xlabel('Pixel AXIS1')
    plt.ylabel('Pixel AXIS2')

    fig.add_subplot(313)
    plt.plot(spec['WAVELENGTH'][0],spec['FLUX'][0],'k')
    plt.xlim(minwav,maxwav)
    plt.ylim(0,3.5e-15)
    plt.xlabel('Wavelength / $\AA$')
    plt.ylabel('Flux, erg s$^{-2}$ cm$^{-2}$ $\AA^{-1}$')
    plt.axvline(x=cen_wav,color='b',linestyle='--')
    plt.axvline(x=cen_wav+del_wav,color='b',linestyle='--')
    plt.axvline(x=cen_wav-del_wav,color='b',linestyle='--')

    fig.add_subplot(312)
    xedges = np.arange(minwav,maxwav,1.0)
    yedges = np.arange(900,1100,2.0)

    H, xedges, yedges = np.histogram2d(data['WAVELENGTH'],data['AXIS2_X'],
                        bins=(xedges, yedges))
    im = plt.imshow(np.arcsinh(H.T), interpolation='nearest', origin='low',
                    aspect='auto', extent=[xedges[0], xedges[-1], yedges[0],
                    yedges[-1]], cmap='binary')
    plt.scatter(data['WAVELENGTH'][mask],data['AXIS2_X'][mask],marker='.',
                s=5,alpha=0.3)
    plt.xlim(minwav,maxwav)
    plt.ylim(1000-50,1000+50)
    plt.axvline(x=cen_wav,color='b',linestyle='-')
    plt.axvline(x=cen_wav+del_wav,color='b',linestyle='--')
    plt.axvline(x=cen_wav-del_wav,color='b',linestyle='--')
    plt.axhline(y=1000,color='r',linestyle='-')
    plt.axhline(y=1000+del_slit,color='r',linestyle='--')
    plt.axhline(y=1000-del_slit,color='r',linestyle='--')

    mask2 = ( ( data['AXIS2_X'] < 1000+del_slit  ) * \
            ( data['AXIS2_X'] > 1000-del_slit ) ) * \
            ( ( data['WAVELENGTH'] > cen_wav-del_wav  ) * \
            ( data['WAVELENGTH'] < cen_wav+del_wav ) )
    mask_ran=np.random.randint(sum(mask2),size=int(.1*sum(mask2)  ) )
    plt.scatter(data['WAVELENGTH'][mask2][mask_ran],
                data['AXIS2_X'][mask2][mask_ran], marker='.', s=5,
                alpha=0.3,color='y')
    plt.xlabel('Wavelength / $\AA$')
    plt.ylabel('Pixel AXIS2')

    fig.add_subplot(311)
    plt.scatter(data_ori['AXIS1'][mask2][mask_ran],
                data_ori['AXIS2'][mask2][mask_ran],marker='.', s=5,
                alpha=0.3, color='g')
    plt.draw()
    plt.show()
    plt.tight_layout(h_pad=-3)
    plt.subplots_adjust(top=0.93)

class Printer():
    """
    Print things to stdout on one line dynamically
    """
    def __init__(self,data):

        sys.stdout.write("\r\x1b[K"+data.__str__())
        sys.stdout.flush()
