This program downloads fits files with fields from SDSS DR12 for specified coordinates.
It needs Python >= 2.6 to run. Also numpy and astropy modules have to be installed.


INSTALLATION

The program does not require any installation process. Just download and unpack
it anywhere you want.


USAGE

The program requires a table with coordinates of desired objects. The table must
contain at least three columns: the first one -- names of objects (they need not to belong to
any catalogue, just string for filenames, it has to be unique), the second one -- RA
in degrees, the third one -- DEC in degrees. If the package is launched with '-a' option, then
the fourth parameter (radius, see below) must be in the input file. By default the name of this file is
'coordinates.dat', but you can specify an arbitrary name (see KEYS section).
    The only one command line argument that is needed is a set of filters. The set of
filters must by specified as a solid string of lowercase letters, for example urz or ugriz.


KEYS

    -i, --input FILENAME
specified a file with table of coordinates (in case of its name is not 'coordinates.dat')

    -a, --adjacent            
If the object is located close to the edge of the field (or if object size is big) it may
be split onto two or more adjacent fields. Setting the -a parameter force the program to
search and download all fields within specified radius around the object. The radius is
stored in fourth column of the input file (units -- arcminutes). 
Note that about one-third of SDSS is overlapped by two or more fields and with -a parameter
the program will download them all. The maximum value of radius is 60 arcminutes.
Be careful with big values of R: it may take a large amount of disk space.

    -s, --swarp
If you turn -a parameter on, you will receive several fields for some objects (they
will be stored in different files with _0, _1, _2, etc endings). You may want to concatenate
all this pieces into one image by using Emmanuel Bertins SWarp package
(http://www.astromatic.net/software/swarp). With -s parameter turned on, the program 
will automatically concatenate fields using SWarp package when downloading of all fields is finished.
SWarp package have to be installed on your system.

    -c, --convert
Different fields in SDSS have somewhat different photometric parameters (GAIN and m0) which makes
incorrect direct concatenating. Before concatenating one may want to convert images to identical
photometric calibration. With -c parameter turned on package will convert all images from SDSS NMgy
units to DN values. (Read http://data.sdss3.org/datamodel/files/BOSS_PHOTOOBJ/frames/RERUN/RUN/CAMCOL/frame.html
for details)

    -f, --free
One can free some disc space by removing original separated fields after concatenation. With this
parameter turned on package will automatically delete these fields. One can keep them by
omitting this parameter.

    -t, --trim
Crop (trim) image of galaxy to size given in the fourth column of the input file. Result is square FITS
file with size 2*r and with galaxy located at the centre of the image.
Note: to run this option one must have the astropy module installed. Try 'pip install astropy' to install it.
See astropy.org for details.

    -p, --ps
Download psFields files (if you want to get a PSF data). Will be stored into downloads/ps directory.

    --scatter
This option forces the package to place all files related to every object into separate
directory with name according to the object name. In other words, when downloading is finished
there will be as many sub-directories in 'downloads' directory as objects in the input file.
The default behaviour (without --scatter option) is to place in the same folder all images
taken in the same filter.

    --add_urls  %file_with_urls%
If for some reason SDSS server does not return all fields in the desired region, one can
add this fields by specifying their full http addresses in the given file. The file must be
organised as follows: object name (the same as in the file with coordinates) followed by
list of urls separated by white spaces (thus one line per object). The filter character
must be replaced with * symbol, script will automatically download these fields for every
specified filter. For example, the next line is a correct line for this file:

m81  http://data.sdss3.org/sas/dr12/boss/photoObj/frames/301/4294/4/frame-*-004294-4-0234.fits.bz2 http://data.sdss3.org/sas/dr12/boss/photoObj/frames/301/4294/5/frame-*-004294-5-0234.fits.bz2

    --add_fields %file_with_fields%
The same as for add_urls option, but run, rerun, camcol, field data are specified
for objects. Format is object name with colon, followed by integer values for run,
rerun, camcol and field gouped for every field and separated by semicolons:

obj1: 	4381 301 5 105; 4381 301 5 105

It contains urls for two additional fields for object 'm81'. Note the asterisk
symbols at filter positions in given urls. 

Important notes about --add_urls option.
1) This option requires wget package to be installed on your system.
2) Use this option with caution. Specified additional fields must be near to the main ones.
SWarp package does not check proximity of the fields to be coadded and creates output
file as large as it needed to fit all the fields. If some of the fields are far from other,
the output file can be HUGE. The main list of fields is guarantied to be inside of the
region specified by user, but additional fields are not being checked.

OUTPUT

The fields will be downloaded in the 'downloads' directory, each filter in separate subdirectory
(downloads/g, downloads/r, etc.). SWarp-scripts (if -s option is activated) will appear
in the same directory with downloaded files.
   The program generates also two ASCII files: 'fiedls.dat' with lists of the downloaded fields
for each object and 'errors_404.dat' with a list of objects which have not been downloaded
due to "not found" or "connection" errors.


EXAMPLES

    python ./sdss_downloader.py gri -i mygalaxies.dat -a -s -c -f

downloades all fields within specified radius for objects from file 'mygalaxies.dat'
in filters g, r and i and concatenates adjacent fields. Original fields will
be removed, when concatenation is done.