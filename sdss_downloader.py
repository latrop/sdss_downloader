#! /usr/bin/env python

################################################################################################
#  Author: Sergey Savchenko (savchenko.s.s@gmail.com)
#                                                                                              #
#  Input: file "coordinates.dat" with two columns RA and DEC in degrees. One row per galaxy.   #
#                                                                                              #
#  Using: first command line argument is a string with list of filters.                        #
#                                                                                              #
#  Example: $> downloader_dr9_th.py gri                                                        #
#                                                                                              #
################################################################################################

from threading import Thread

from math import hypot
from time import time, sleep
import urllib2
import sys
import os

def findField2(objRA, objDEC):
    radius = 30
    request = "http://skyserver.sdss3.org/dr9/en/tools/search/x_radial.asp?ra=%1.5f&dec=%1.5f&radius=10&min_u=0&max_u=20&min_g=0&max_g=20&min_r=0&max_r=20&min_i=0&max_i=20&min_z=0&max_z=20&entries=top&topnum=10&format=csv" % (objRA, objDEC)
    u = urllib2.urlopen(request)
    table = u.read().split("\n")
    minDist = 1e10
    for line in table[1:-1]:
        objid, run, rerun, camcol, field, obj, t, ra, dec, u, g, r, i, z, E_u, E_g, E_r, E_i, Err_z = [float(i) for i in line.split(',')]
        dist = hypot(objRA -ra, objDEC - dec)
        if dist < minDist:
            minDist = int(dist)
            optRun = int(run)
            optRerun = int(rerun)
            optCamcol = int(camcol)
            optField = int(field)
    return optRun, optRerun, optCamcol, optField

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
    if file_name is None:
        file_name = url.split('/')[-1].split('.')[0]
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


bandlist = sys.argv[1]

for band in bandlist:
    if not os.path.exists("./downloads/%s" % (band)):
        os.makedirs("./downloads/%s" % (band))

listOfCoords = open("coordinates.dat").readlines()
counter = 0
errFile = open("errors_404.dat", "w", buffering=0)
errFile.truncate(0)
with open("coords2fields.dat", "w", buffering=0) as outFile:
    outFile.truncate(0)
    for line in listOfCoords:
        startTime = time()
        counter += 1
        params = line.split()
        if len(params) == 3:
            galName = params[0]
            ra = float(params[1])
            dec = float(params[2])
            print "\033[33m Downloading field for %1.5f %1.5f: '%s'  (%i/%i) \033[0m" % (ra, dec, galName, counter, len(listOfCoords))
        else:
            ra = float(params[0])
            dec = float(params[1])
            galName = None
            print "\033[33m Downloading field for %1.5f %1.5f  (%i/%i) \033[0m" % (ra, dec, counter, len(listOfCoords))
        run, rerun, camcol, field = findField2(ra, dec)
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
            errFile.write("%1.6f  %1.6f   %s \n" % (ra, dec, url.split("/")[-1][8:-9]))
            continue
        downloadThreads = []
        for band in bandlist:
            dth = Thread(target=download, args=(urls[band], band, galName))
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
        if galName is None:
            galName = url.split("/")[-1][8:-9]
        print "         \033[34m [DONE] \033[0m    (%1.2f sec)\n\n" % (time()-startTime)
        outFile.write("%1.6f  %1.6f   %s \n" % (ra, dec, galName))

