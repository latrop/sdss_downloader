This program downloads fits files with fields from SDSS DR9 for specified coordinates.
It needs only Python >= 2.6 to run, nothing else.


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

    -s, --script	     
If you turn the -a parameter on, you will receive several fields for some objects (they
will be stored in different files with _0, _1, _2... endings). You may want to concatenate
all this pieces into one image by using Emmanuel Bertins SWarp program 
(http://www.astromatic.net/software/swarp). With -s parameter turned on, the program 
generates a simpe Bash script which runs the SWarp program for all this fields. Just run it
when downloading is completed (the script takes the downloaded files as it is, you must
not rename or unpack them). SWarp package have to be installed on your system.

    -p, --ps
Download psFields files (if you want to get a PSF data). Will be stored into downloads/ps directory.

OUTPUT

The fields will be downloaded in the 'downloads' directory, each filter in separate subdirectory
(downloads/g, downloads/r, etc.). SWarp-scripts (if -s option is activated) will appear
in the same directory with downloaded files.
   The program generates also two ASCII files: 'fiedls.dat' with lists of the downloaded fields
for each object and 'errors_404.dat' with a list of objects which have not been downloaded
due to "not found" or "connection" errors.


EXAMPLES

    python ./sdss_downloader.py gri -i mygalaxies.dat -a -s

downloades all fields within specified radius for objects from file 'mygalaxies.dat'
in filters g, r and i and creates scripts to concatenate adjacent fields.