# stis_photons
Assing wavelengths to HST/STIS TIME_TAG event files



<<<<<<< HEAD
##Example:

 ```
 import stis_photons as stis

fln_tag = 'o4oj01c30_tag.fits'	## TAG-TIME file fln_x1d = 'o4oj01c30_x1d.fits'	## Spectra associated with TIME-TAG File fln_moc = 'h4s1350io_moc.fits'	## MOC file used in the data reduction fln_dsp = 'm7p16110o_dsp.fits' ## DSP file used in the data reduction output = 'stis_timetag'	## Output numpy file

stis.localise(fln_tag,fln_x1d,fln_dsp,fln_moc, output=output, transform = True, verbose=True, plot_diagnostics=True, del_slit = 10,del_line = 40, del_wav = 15)

```

##Acknowledgments:
=======
Example

::code::
import stis_photons as stis

fln_tag  = 'o4oj01c30_tag.fits'		## TAG-TIME file
fln_x1d  = 'o4oj01c30_x1d.fits'		## Spectra associated with TIME-TAG File
fln_moc  = 'h4s1350io_moc.fits'		## MOC file used in the data reduction
fln_dsp  = 'm7p16110o_dsp.fits'   ## DSP file used in the data reduction
output 	 = 'stis_timetag'	## Output numpy file

stis.localise(fln_tag,fln_x1d,fln_dsp,fln_moc, output=output, transform = True,
            verbose=True, plot_diagnostics=True,
            del_slit = 10,del_line = 40, del_wav = 15)
::code::
>>>>>>> origin/master
