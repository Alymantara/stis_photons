# stis_photons
Assign wavelengths to HST/STIS TIME_TAG event files

## Barycentre Photons
It is necessary to retrieve the ephemeris of HST to correct the time stamp
of each photon to the solar-system barycenter. This will correct for (from
ODELAYTIME help)
* general relativistic effects (up to 2 milliseconds),
* displacement of the telescope from the center of the Earth (up
    to 20 milliseconds),
* and displacement of the Earth from the solar-system barycenter
    (up to 500 seconds).

### Retrieve ephemeris.
Go directly to [HST Starview](https://starview.stsci.edu/web/). On the left hand
side, choose >ALL>ENGINEERING. In the search form, choose in the ARchive class:
ORB, to search for orbital ephemeris. Click search and look for the files that
contain your observations (every _ORB_ file usually contains three days worth of
ephemeris). From the dataset than you retrieve within the HST archive, the one
we need is the _FIT_ file.

### Correction
Once in iraf/pyraf, load >stdas/hst_cal/stis. In there, there will be a routine
called ODELAYTIME. To correct each _tag.fits_, simply run the routine as an
example:

```
odelay odfa01030_tag.fits obs_ephem=p6d0000r.fit
```
This will overwrite the original file (I usually keep a raw copy of these files
just in case.)


## stis_photons example:

 ```
import stis_photons

fln_tag = 'o4oj01020_tag.fits'	## TAG-TIME file
fln_x1d = 'o4oj01020_x1d.fits'	## Spectra associated with TIME-TAG File
fln_moc = 'h4s1350io_moc.fits'	## MOC file used in the data reduction
fln_dsp = 'm7p16110o_dsp.fits' ## DSP file used in the data reduction
output = 'stis_timetag20'	## Output numpy and text file

stis_photons.localise(fln_tag,fln_x1d,fln_dsp,fln_moc, output=output,
              verbose=True, plot_diagnostics=True, del_slit = 10,del_line = 40,
              del_wav = 15)
```
If you have already the fits file with wavelengths, you can also plot the 
diagnostics:
```
stis_photons.plotter('stis_timetag.fits',fln_x1d,fln_dsp,fln_moc)
```
