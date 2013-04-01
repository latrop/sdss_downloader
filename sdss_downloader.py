#! /usr/bin/env python


#                                
#  Author: Sergey Savchenko (savchenko.s.s@gmail.com)
#


from threading import Thread

from math import hypot
from time import time, sleep
import urllib2
import sys
import os

import argparse

def findField2(objRA, objDEC, radius):
    if radius < 0.001:
        request = "http://skyserver.sdss3.org/dr9/en/tools/search/x_radial.asp?ra=%1.5f&dec=%1.5f&radius=2&min_u=0&max_u=20&min_g=0&max_g=20&min_r=0&max_r=20&min_i=0&max_i=20&min_z=0&max_z=20&entries=top&topnum=10&format=csv" % (objRA, objDEC)
    else:
        request = "http://skyserver.sdss3.org/dr9/en/tools/search/x_radial.asp?ra=%1.5f&dec=%1.5f&radius=%1.2f&min_u=0&max_u=20&min_g=0&max_g=20&min_r=0&max_r=20&min_i=0&max_i=20&min_z=0&max_z=20&entries=top&topnum=500&format=csv" % (objRA, objDEC, radius)
    u = urllib2.urlopen(request)
    table = u.read().split("\n")
    if len(table) < 2:
        return [(-1, -1, -1, 1)], -1
    # Find the nearest object and the corresponding field
    minDist = 1e10
    for line in table[1:-1]:
        objParams = line.split(',')
        objID  = int(objParams[0])
        run    = int(objParams[1])
        rerun  = int(objParams[2])
        camcol = int(objParams[3])
        field  = int(objParams[4])
        ra     = float(objParams[7])
        dec    = float(objParams[8])
        dist   = hypot(objRA -ra, objDEC - dec)
        if dist < minDist:
            minDist = dist
            optObjID = objID
            optRun = run
            optRerun = rerun
            optCamcol = camcol
            optField = field
    fList = [(optRun, optRerun, optCamcol, optField)]
    if radius < 0.001:
        return fList, optObjID
    else:
        for line in table[1:-1]:
            objParams = line.split(',')
            run    = int(objParams[1])
            rerun  = int(objParams[2])
            camcol = int(objParams[3])
            field  = int(objParams[4])
            if (run, rerun, camcol, field) not in fList:
                fList.append((run, rerun, camcol, field))
    return fList, optObjID

def getUrl(run, rerun, camcol, field, band):
    return "http://data.sdss3.org/sas/dr9/boss/photoObj/frames/%s/%i/%i/frame-%s-%.6i-%i-%.4i.fits.bz2" % (
        rerun, run, camcol, band, run, camcol, field)


def testUrlExists(url):
    try:
        u = urllib2.urlopen(url)
        code = u.code
        meta = u.info()
        file_size = int(meta.getheaders("Content-Length")[0])
        if code == 200:
            return file_size
        return -1
    except urllib2.HTTPError:
        return -1
        

def download(url, passband, file_name):
    try:
        u = urllib2.urlopen(url)
    except urllib2.HTTPError:
        print ""
        return -1
    meta = u.info()
    file_size = int(meta.getheaders("Content-Length")[0])
    f = open("./downloads/%s/%s_%s.fits.bz2" % (passband, file_name, passband), 'wb')
    buff = u.read(file_size)
    f.write(buff)
    f.close()

def threadsAlive(listOfThreads):
    for th in listOfThreads:
        if th.isAlive():
            return True
    return False

parser = argparse.ArgumentParser(description="Download fits-files of fields for specified coordinates")
parser.add_argument("filters", help="List of filters (for example gri or uiz)")
parser.add_argument("-i", "--input", default="coordinates.dat", help="File with coordinates to download")
parser.add_argument("-r", "--radius", default=0.0, help="Download ajacent fields if any exist in r arcmin (max 60 arcmin)", type=float)
parser.add_argument("-s", "--script", action="store_true", default=False, help="Create a script for concatenation by SWarp.")
args = parser.parse_args()

if args.radius > 60:
    print "Radius must by less than 60 arcmin"
    sys.exit(1)

bandlist = args.filters
for band in bandlist:
    if not os.path.exists("./downloads/%s" % (band)):
        os.makedirs("./downloads/%s" % (band))

if args.script:
    for band in bandlist:
        f = open("./downloads/%s/concatenate.sh" % (band), "w")
        f.truncate(0)
        f.close()


listOfCoords = open(args.input).readlines()
counter = 0
errFile = open("errors_404.dat", "w", buffering=0)
errFile.truncate(0)
with open("fields.dat", "w", buffering=0) as outFile:
    outFile.truncate(0)
    for line in listOfCoords:
        if line.startswith("#") or line.startswith("!"):
            continue
        startTime = time()
        counter += 1
        params = line.split()

        if len(params) == 3:
            galName = params[0]
            ra = float(params[1])
            dec = float(params[2])
            print "\033[33m Downloading field for %1.5f %1.5f: '%s'  (%i/%i) \033[0m" % (ra, dec, galName, counter, len(listOfCoords))
        else:
            print "%s must contain three columns: name, RA and DEC" % (args.input)
            sys.exit(1)
        objFieldList, objID = findField2(ra, dec, args.radius)

        if len(objFieldList) > 1:
            print "There are %i files for this object" % (len(objFieldList))
        outFile.write("%s  %1.6f  %1.6f  " % (galName, ra, dec))
        for ifield in xrange(len(objFieldList)):
            if len(objFieldList) > 1:
                print "Downloading (%i/%i)" % (ifield + 1, len(objFieldList))
                curGalName = galName + "_" + str(ifield)
            else:
                curGalName = galName

            run, rerun, camcol, field = objFieldList[ifield]
            # Check if fits files exist for all filters:
            print " Checking file existense:"
            allExist = True
            urls = {}
            for band in bandlist:
                print "                    %s" % (band),
                url = getUrl(run, rerun, camcol, field, band)
                answer = testUrlExists(url)
                if answer == -1:
                    allExist = False
                    print "\033[31m [Fail!] \033[0m\n"
                    break
                print "\033[32m [OK] \033[0m   (%i bytes)" % (answer)
                urls[band] = url
            if not allExist:
                errFile.write("%s  %1.6f  %1.6f \n" % (galName, ra, dec))
                continue
            downloadThreads = []
            for band in bandlist:
                dth = Thread(target=download, args=(urls[band], band, curGalName))
                dth.daemon = True
                dth.start()
                downloadThreads.append(dth)
            print " Downloading",
            while threadsAlive(downloadThreads):
                sys.stdout.write(".")
                sys.stdout.flush()
                sleep(0.333)
                sys.stdout.write(".")
                sys.stdout.flush()
                sleep(0.333)
                sys.stdout.write(".")
                sys.stdout.flush()
                sleep(0.333)
                sys.stdout.write("\b\b\b   \b\b\b")
                sys.stdout.flush()
                sleep(0.333)
            print "         \033[34m [DONE] \033[0m    (%1.2f sec)\n\n" % (time()-startTime)
            outFile.write(" %s.fits " % (curGalName))
        outFile.write("\n")
        if args.script and (len(objFieldList) > 1):
            for band in bandlist:
                f = open("./downloads/%s/concatenate.sh" % (band), "a")
                for i in xrange(len(objFieldList)):
                    f.write("bunzip2 %s_%i_%s.fits.bz2 \n" % (galName, i, band))
                f.write("swarp ")
                for i in xrange(len(objFieldList)):
                    f.write(" %s_%i_%s.fits[0] " % (galName, i, band))
                f.write("\n")
                f.write("mv ./coadd.fits ./%s_%s.fits \n" % (galName, band))
                for i in xrange(len(objFieldList)):
                    f.write("rm %s_%i_%s.fits \n" % (galName, i, band))
                f.write("rm coadd.weight.fits \n \n")
                f.close()
