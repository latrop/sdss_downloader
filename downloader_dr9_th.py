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

from time import time, sleep
import urllib2
import sys
import os

from numpy import load, hypot, where
from scipy.ndimage import minimum_position

def findField(objRA, objDEC):
    dists = hypot(objRA - gridRA, objDEC - gridDEC)
    inds = dists.argsort()
    for i in inds:
        run = runList[i]
        rerun = rerunList[i]
        camcol = camcolList[i]
        field = fiedlList[i]
        u = getUrl(run, rerun, camcol, field, 'g')
        if dists[i] > 0.1392:
            run, rerun, camcol, field = -1, -1, -1, -1
            break
        if (testUrlExists(u) > 0):
            break
    return run, rerun, camcol, field


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


gridRA = load("./data/ra.npy")
gridDEC = load("./data/dec.npy")
runList = load("./data/run.npy")
rerunList = load("./data/rerun.npy")
camcolList = load("./data/camcol.npy")
fiedlList = load("./data/field.npy")
mu_startList = load("./data/mu_start.npy")
mu_endList = load("./data/mu_end.npy")
nu_startList = load("./data/nu_start.npy")
nu_endList = load("./data/nu_end.npy")




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
        run, rerun, camcol, field = findField(ra, dec)
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

