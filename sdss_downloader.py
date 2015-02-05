#! /usr/bin/env python

#                                
#  Author: Sergey Savchenko (savchenko.s.s@gmail.com)
#

from threading import Thread
import subprocess
from math import hypot, log10, sqrt
from time import time, sleep
import urllib2
import sys
import os
import shutil
import argparse

import numpy as np
import pyfits

def findField2(objRA, objDEC, radius):
    if radius < 0.001:
        request = "http://skyserver.sdss.org/dr12/en/tools/search/x_radial.aspx?ra=%1.5f&dec=%1.5f&radius=2&min_u=0&max_u=20&min_g=0&max_g=20&min_r=0&max_r=20&min_i=0&max_i=20&min_z=0&max_z=20&entries=top&topnum=100&format=csv" % (objRA, objDEC)
    else:
        request = "http://skyserver.sdss.org/dr12/en/tools/search/x_radial.aspx?ra=%1.5f&dec=%1.5f&radius=%1.2f&min_u=0&max_u=20&min_g=0&max_g=20&min_r=0&max_r=20&min_i=0&max_i=20&min_z=0&max_z=20&entries=top&topnum=10000&format=csv" % (objRA, objDEC, radius)
    u = urllib2.urlopen(request)
    table = u.read().split("\n")
    if len(table) < 2:
        return [(-1, -1, -1, 1)], -1
    # Find the nearest object and the corresponding field
    minDist = 1e10
    for line in table:
        # ignore comments, header and footer of table
        if (len(line)==0) or (not line[0].isdigit()):
            continue
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
        for line in table:
            # ignore comments, header and footer of table
            if (len(line)==0) or (not line[0].isdigit()):
                continue
            objParams = line.split(',')
            run    = int(objParams[1])
            rerun  = int(objParams[2])
            camcol = int(objParams[3])
            field  = int(objParams[4])
            if (run, rerun, camcol, field) not in fList:
                fList.append((run, rerun, camcol, field))
    return fList, optObjID


def getUrl(run, rerun, camcol, field, band):
    return "http://dr12.sdss3.org/sas/dr12/boss/photoObj/frames/%s/%i/%i/frame-%s-%.6i-%i-%.4i.fits.bz2" % (
        rerun, run, camcol, band, run, camcol, field)


def getUrlPs(run, rerun, camcol, field):
    return "http://dr12.sdss3.org/sas/dr12/env/PHOTO_REDUX/%s/%i/objcs/%i/psField-%.6i-%i-%.4i.fit" % (
        rerun, run, camcol, run, camcol, field)


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
    if passband == "ps":
        f = open("./downloads/%s/%s_%s.fits" % (passband, file_name, passband), 'wb')
    else:
        f = open("./downloads/%s/%s_%s.fits.bz2" % (passband, file_name, passband), 'wb')
    buff = u.read(file_size)
    f.write(buff)
    f.close()


def threadsAlive(listOfThreads):
    for th in listOfThreads:
        if th.isAlive():
            return True
    return False


### some sdss functions are below
def prep_ima(gal):
    new_gal = "./prep_ima_tmp.fits"
    # This function re-writes pixel values given in NMgy to ADU
    hdulist = pyfits.open(gal, do_not_scale_image_data=True, mode='update')
    img = hdulist[0].data
    (dimy,dimx) = img.shape
    cimg = np.zeros(shape=(dimy,dimx))
    nrowc = dimy
    calib = hdulist[1].data
    for i in range(dimy):
        cimg[i] = calib 
    dn = img/cimg
    shutil.copy(gal, new_gal) 
    hdulist1 = pyfits.open(new_gal, do_not_scale_image_data=True, mode='update')
    img_new = hdulist1[0].data
    for i in range(dimy):
        for k in range(dimx):
            img_new[i,k] = dn[i,k]  #new fits-file img_new with ADU    
    hdulist1.flush()
    os.remove(gal)
    shutil.move(new_gal, gal)

def gain_dark_SDSS(camcol,band,run):
    if band=='u':
        if camcol==1:
            gaidark = [1.62,9.61]
        if camcol==2:
            if run<1100:
                gaidark = [1.595,12.6025]        
            else:
                gaidark = [1.825,12.6025]
        if camcol==3:
            gaidark = [1.59,8.7025]
        if camcol==4:
            gaidark = [1.6,12.6025]
        if camcol==5:
            gaidark = [1.47,9.3025]
        if camcol==6:
            gaidark = [2.17,7.0225]
    if band=='g':
        if camcol==1:
            gaidark = [3.32,15.6025]
        if camcol==2:
            gaidark = [3.855,1.44]        
        if camcol==3:
            gaidark = [3.845,1.3225]
        if camcol==4:
            gaidark = [3.995,1.96]
        if camcol==5:
            gaidark = [4.05,1.1025]
        if camcol==6:
            gaidark = [4.035,1.8225]
    if band=='r':
        if camcol==1:
            gaidark = [4.71,1.8225]
        if camcol==2:
            gaidark = [4.6,1.00]        
        if camcol==3:
            gaidark = [4.72,1.3225]
        if camcol==4:
            gaidark = [4.76,1.3225]
        if camcol==5:
            gaidark = [4.725,0.81]
        if camcol==6:
            gaidark = [4.895,0.9025]
    if band=='i':
        if camcol==1:
            gaidark = [5.165,7.84]
        if camcol==2:
            if run<1500:
                gaidark = [6.565,5.76]
            if run>1500:
                gaidark = [6.565,6.25]            
        if camcol==3:
            gaidark = [4.86,4.6225]
        if camcol==4:
            if run<1500:
                gaidark = [4.885,6.25]
            if run>1500:
                gaidark = [4.885,7.5625]    
        if camcol==5:
            gaidark = [4.64,7.84]
        if camcol==6:
            gaidark = [4.76,5.0625]

    if band=='z':
        if camcol==1:
            gaidark = [4.745,0.81]
        if camcol==2:
            gaidark = [5.155,1.0]        
        if camcol==3:
            gaidark = [4.885,1.0]
        if camcol==4:
            if run<1500:
                gaidark = [4.775,9.61]
            if run>1500:
                gaidark = [4.775,12.6025]    
        if camcol==5:
            if run<1500:
                gaidark = [3.48,1.8225]
            if run>1500:
                gaidark = [3.48,2.1025]    
        if camcol==6:
            gaidark = [4.69,1.21]
    return gaidark[0],gaidark[1]

def SDSS_dr8(gal_image):
    #http://data.sdss3.org/datamodel/files/BOSS_PHOTOOBJ/frames/RERUN/RUN/CAMCOL/frame.html
    hdulist = pyfits.open(gal_image)
    prihdr = hdulist[0].header
    run = int(prihdr['RUN'])
    band = str(prihdr['FILTER'])
    camcol = int(prihdr['CAMCOL'])
    kkk = prihdr['NMGY']
    m0 = 22.5 - 2.5*log10(kkk)        
    GAIN,var = gain_dark_SDSS(camcol,band,run)
    read_out_noise = sqrt(var)*GAIN    # should be <~5 e
    return GAIN, read_out_noise, m0

# parsing the argument line here
parser = argparse.ArgumentParser(description="Download fits-files of fields for specified coordinates")
parser.add_argument("filters", help="List of filters (for example gri or uiz)")
parser.add_argument("-i", "--input", default="coordinates.dat", help="File with coordinates to download")
parser.add_argument("-a", "--adjacent", action="store_true", 
                    default=False, help="Download adjacent fields if any exist")
parser.add_argument("-s", "--swarp", action="store_true", default=False, help="Concatenate fields using SWarp package")
parser.add_argument("-c", "--convert", action="store_true", default=False, help="Convert fields to the same photometry zeropoint")
parser.add_argument("-f", "--free", action="store_true", default=False, help="Remove individual fields after concatenation")
parser.add_argument("-t", "--trim", action="store_true", default=False, help="Crop image to galaxy size")
parser.add_argument("-p", "--ps", action="store_true", default=False, help="Download psField files")
args = parser.parse_args()


# Make dirs for all bands
bandlist = args.filters
for band in bandlist:
    if not os.path.exists("./downloads/%s" % (band)):
        os.makedirs("./downloads/%s" % (band))

# Make dir for ps
if args.ps:
    if not os.path.exists("./downloads/ps/"):
        os.makedirs("./downloads/ps/")

listOfCoords = open(args.input).readlines()
counter = 0
errFile = open("errors_404.dat", "w", buffering=0)
errFile.truncate(0)
with open("fields.dat", "w", buffering=0) as outFile:
    outFile.truncate(0)
    for line in listOfCoords:
        if line.startswith("#") or line.startswith("!"):
            continue
        counter += 1
        params = line.split()

        if (len(params) in (3, 4)) or ((len(params) == 4) and (args.adjacent is True)):
            galName = params[0]
            ra = float(params[1])
            dec = float(params[2])
            if args.adjacent is True:
                r_adj = float(params[3])
            else:
                r_adj = 0.0
            print "\033[33m Downloading field for %1.5f %1.5f: '%s'  (%i/%i) \033[0m" % (ra, dec, galName, counter, len(listOfCoords))
        else:
            print "Invalid number of columns in input file %s" % args.input
            sys.exit(1)
        objFieldList, objID = findField2(ra, dec, r_adj)

        if len(objFieldList) > 1:
            print "There are %i files for this object" % (len(objFieldList))
        outFile.write("%s  %1.6f  %1.6f  " % (galName, ra, dec))
        for ifield in xrange(len(objFieldList)):
            startTime = time()
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
            if args.ps:
                # Check if ps file exists
                print "                   ps",
                urlPs = getUrlPs(run, rerun, camcol, field)
                answer = testUrlExists(urlPs)
                if answer == -1:
                    allExist = False
                    print "\033[31m [Fail!] \033[0m\n"
                else:
                    print "\033[32m [OK] \033[0m   (%i bytes)" % (answer)
            if not allExist:
                errFile.write("%s  %1.6f  %1.6f \n" % (galName, ra, dec))
                continue
            downloadThreads = []
            # Downloading files in all passbands
            for band in bandlist:
                dth = Thread(target=download, args=(urls[band], band, curGalName))
                dth.daemon = True
                dth.start()
                downloadThreads.append(dth)
            # Downloading ps file
            if args.ps:
                dth = Thread(target=download, args=(urlPs, "ps", curGalName))
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

        # Concatenating fields
        if args.swarp and (len(objFieldList) > 1):
            for band in bandlist:
                GAINList = []
                READOUTList = []
                m0List = []
                for i in xrange(len(objFieldList)):
                    subprocess.call("bunzip2 ./downloads/%s/%s_%i_%s.fits.bz2" % (band, galName, i, band), shell=True)
                    if args.convert:
                        prep_ima("./downloads/%s/%s_%i_%s.fits" % (band, galName, i, band))
                        GAIN, READOUT, m0 = SDSS_dr8("./downloads/%s/%s_%i_%s.fits" % (band, galName, i, band))
                        GAINList.append(GAIN)
                        READOUTList.append(READOUT)
                        m0List.append(m0)
                print "Running SWarp for %s band..." % (band)
                callSwarpString = "swarp -verbose_type quiet "+" ".join([("./downloads/%s/%s_%i_%s.fits[0]"%(band, galName, i, band))
                                                     for i in xrange(len(objFieldList))])
                subprocess.call(callSwarpString, shell="True")
                shutil.move("./coadd.fits", "./downloads/%s/%s_%s.fits"%(band, galName, band))
                os.remove("coadd.weight.fits")
                os.remove("swarp.xml")
                if args.free:
                    for i in xrange(len(objFieldList)):
                        os.remove("./downloads/%s/%s_%i_%s.fits"%(band, galName, i, band))
                # store mean keywords to coadded file
                hdu = pyfits.open("./downloads/%s/%s_%s.fits"%(band, galName, band), do_not_scale_image_data=True, mode="update")
                header = hdu[0].header
                header["GAIN"] = np.mean(GAINList)
                header["READOUT"] = np.mean(READOUTList)
                header["M0"] = np.mean(m0List)
                hdu.flush()
                # Crop images
                if args.trim:
                    print "Cropping..."
                    # we need astropy to run this option
                    try:
                        from astropy import wcs
                    except ImportError:
                        print "Error: astropy module was not found."
                        print "Cropping is impossible."
                        print "Install astropy or remove -t option."
                        exit(1)
                    hdu = pyfits.open("./downloads/%s/%s_%s.fits"%(band, galName, band))
                    data = hdu[0].data
                    header = hdu[0].header
                    wcs = wcs.WCS("./downloads/%s/%s_%s.fits"%(band, galName, band))
                    xGalPix, yGalPix = wcs.wcs_world2pix([[ra, dec]], 1)[0]
                    size = min(r_adj*60.0/0.396127, xGalPix, yGalPix, data.shape[1]-xGalPix, data.shape[0]-yGalPix)
                    xCropMin = xGalPix - size
                    xCropMax = xGalPix + size
                    yCropMin = yGalPix - size
                    yCropMax = yGalPix + size
                    data = data[yCropMin:yCropMax, xCropMin:xCropMax]
                    #print name, int(xGalPix), int(yGalPix), "|", int(xCropMin), int(xCropMax), int(yCropMin), int(yCropMax), sizes[name]
                    outHDU = pyfits.PrimaryHDU(data=data, header=header)
                    filename = "./downloads/%s/%s_%s_trim.fits"%(band, galName, band)
                    if os.path.exists(filename):
                        os.remove(filename)
                    outHDU.writeto(filename)
