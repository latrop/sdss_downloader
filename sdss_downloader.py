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
import glob
import shutil
import argparse
import bz2

import numpy as np
from astropy.io import fits as pyfits

try:
    from astropy.wcs import WCS
    wcsOK = True
except ImportError:
    wcsOK = False


def move(src, dst):
    if not os.path.isfile(src):
        print "File %s not found and cannot be moved"
        return
    shutil.copy(src, dst)
    os.remove(src)



    
def findField2(objRA, objDEC, radius):
    request = "http://skyserver.sdss.org/dr13/en/tools/search/x_results.aspx?"
    request += "searchtool=Radial&uband=&gband=&rband=&iband=&zband=&jband=&hband=&kband="
    request += "&TaskName=Skyserver.Search.Radial&ReturnHtml=true&whichphotometry=optical&"
    request += "coordtype=equatorial&ra=%1.5f&dec=%1.5f" % (objRA, objDEC)
    if radius < 0.01:
        request += "&radius=2"
    else:
        request += "&radius=%1.2f" % (radius)
    request += "&min_u=0&max_u=20&min_g=0&max_g=20&min_r=0&max_r=20&min_i=0&max_i=20"
    request += "&min_z=0&max_z=20&min_j=0&max_j=20&min_h=0&max_h=20&min_k=0&max_k=20"
    request += "&format=csv&TableName=&limit=99999"

    u = urllib2.urlopen(request)
    table = u.read().split("\n")
    optRun = None
    optRerun = None
    optCamcol = None
    optField = None
    optObjID = None
    # Find the nearest object and the corresponding field
    minDist = 1e10
    for line in table:
        # ignore comments, header and footer of table
        if (len(line) == 0) or (not line[0].isdigit()):
            continue
        objParams = line.split(',')
        objID = int(objParams[0])
        run = int(objParams[1])
        rerun = int(objParams[2])
        camcol = int(objParams[3])
        field = int(objParams[4])
        ra = float(objParams[7])
        dec = float(objParams[8])
        dist = hypot(objRA - ra, objDEC - dec)
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
            if (len(line) == 0) or (not line[0].isdigit()):
                continue
            objParams = line.split(',')
            run = int(objParams[1])
            rerun = int(objParams[2])
            camcol = int(objParams[3])
            field = int(objParams[4])
            if (run, rerun, camcol, field) not in fList:
                fList.append((run, rerun, camcol, field))
    return fList, optObjID


def getUrl(run, rerun, camcol, field, band):
    u = "http://dr13.sdss.org/sas/dr13/eboss/photoObj/frames/"
    u += "%s/%i/%i/frame-%s-%.6i-%i-%.4i.fits.bz2" % (rerun, run, camcol, band,
                                                      run, camcol, field)
    return u


def getUrlPs(run, rerun, camcol, field):
    u = "http://dr12.sdss3.org/sas/dr12/env/PHOTO_REDUX/"
    u += "%s/%i/objcs/%i/psField-%.6i-%i-%.4i.fit" % (rerun, run, camcol,
                                                      run, camcol, field)
    return u


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
    fName = "./downloads/%s_%s.fits" % (file_name, passband)
    if passband != "ps":  # PSF-files are not compressed on SDSS servers
        fName += ".bz2"
    f = open(fName, 'wb')
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
    (dimy, dimx) = img.shape
    cimg = np.zeros(shape=(dimy, dimx))
    calib = hdulist[1].data
    for i in range(dimy):
        cimg[i] = calib
    dn = img/cimg
    shutil.copy(gal, new_gal)
    hdulist1 = pyfits.open(new_gal, do_not_scale_image_data=True,
                           mode='update')
    img_new = hdulist1[0].data
    for i in range(dimy):
        for k in range(dimx):
            img_new[i, k] = dn[i, k]  # new fits-file img_new with ADU
    hdulist1.flush()
    os.remove(gal)
    move(new_gal, gal)


def change_m0(fitsName, oldM0Value, refM0):
    # function changes value of magnitude zeropoint to a given one
    hdulist = pyfits.open(fitsName)
    data = hdulist[0].data
    header = hdulist[0].header
    deltaM0 = refM0 - oldM0Value
    data = data * (10.0**(0.4*deltaM0))
    outHDU = pyfits.PrimaryHDU(data=data, header=header)
    outHDU.writeto("tmp.fits")
    os.remove(fitsName)
    move("tmp.fits", fitsName)


def gain_dark_SDSS(camcol, band, run):
    if band == 'u':
        if camcol == 1:
            gaidark = [1.62, 9.61]
        if camcol == 2:
            if run < 1100:
                gaidark = [1.595, 12.6025]
            else:
                gaidark = [1.825, 12.6025]
        if camcol == 3:
            gaidark = [1.59, 8.7025]
        if camcol == 4:
            gaidark = [1.6, 12.6025]
        if camcol == 5:
            gaidark = [1.47, 9.3025]
        if camcol == 6:
            gaidark = [2.17, 7.0225]
    if band == 'g':
        if camcol == 1:
            gaidark = [3.32, 15.6025]
        if camcol == 2:
            gaidark = [3.855, 1.44]
        if camcol == 3:
            gaidark = [3.845, 1.3225]
        if camcol == 4:
            gaidark = [3.995, 1.96]
        if camcol == 5:
            gaidark = [4.05, 1.1025]
        if camcol == 6:
            gaidark = [4.035, 1.8225]
    if band == 'r':
        if camcol == 1:
            gaidark = [4.71, 1.8225]
        if camcol == 2:
            gaidark = [4.6, 1.00]
        if camcol == 3:
            gaidark = [4.72, 1.3225]
        if camcol == 4:
            gaidark = [4.76, 1.3225]
        if camcol == 5:
            gaidark = [4.725, 0.81]
        if camcol == 6:
            gaidark = [4.895, 0.9025]
    if band == 'i':
        if camcol == 1:
            gaidark = [5.165, 7.84]
        if camcol == 2:
            if run < 1500:
                gaidark = [6.565, 5.76]
            if run > 1500:
                gaidark = [6.565, 6.25]
        if camcol == 3:
            gaidark = [4.86, 4.6225]
        if camcol == 4:
            if run < 1500:
                gaidark = [4.885, 6.25]
            if run > 1500:
                gaidark = [4.885, 7.5625]
        if camcol == 5:
            gaidark = [4.64, 7.84]
        if camcol == 6:
            gaidark = [4.76, 5.0625]

    if band == 'z':
        if camcol == 1:
            gaidark = [4.745, 0.81]
        if camcol == 2:
            gaidark = [5.155, 1.0]
        if camcol == 3:
            gaidark = [4.885, 1.0]
        if camcol == 4:
            if run < 1500:
                gaidark = [4.775, 9.61]
            if run > 1500:
                gaidark = [4.775, 12.6025]
        if camcol == 5:
            if run < 1500:
                gaidark = [3.48, 1.8225]
            if run > 1500:
                gaidark = [3.48, 2.1025]
        if camcol == 6:
            gaidark = [4.69, 1.21]
    return gaidark[0], gaidark[1]


def SDSS_dr8(gal_image):
    # http://data.sdss3.org/datamodel/files/BOSS_PHOTOOBJ/
    # frames/RERUN/RUN/CAMCOL/frame.html
    hdulist = pyfits.open(gal_image)
    prihdr = hdulist[0].header
    run = int(prihdr['RUN'])
    band = str(prihdr['FILTER'])
    camcol = int(prihdr['CAMCOL'])
    kkk = prihdr['NMGY']
    m0 = 22.5 - 2.5*log10(kkk)
    GAIN, var = gain_dark_SDSS(camcol, band, run)
    read_out_noise = sqrt(var)*GAIN  # should be <~5 e
    return GAIN, read_out_noise, m0


def bunzip(zipName):
    bzipFile = bz2.BZ2File(zipName)
    data = bzipFile.read()
    outName = zipName[:-4]
    if os.path.isfile(outName):
        os.remove(outName)
    outFile = open(outName, 'wb')
    outFile.write(data)
    outFile.close()
    os.remove(zipName)


def reduce_to_same_m0(listOfImages):
    GAINList = []
    READOUTList = []
    m0List = []
    for i, fName in enumerate(listOfImages):
        bunzip("%s.bz2" % fName)
        prep_ima(fName)
        GAIN, READOUT, m0 = SDSS_dr8(fName)
        if i == 0:
            refM0 = m0
        else:
            change_m0(fName, m0, refM0)
        GAINList.append(GAIN)
        READOUTList.append(READOUT)
        hdu = pyfits.open(fName, do_not_scale_image_data=True,
                          mode="update")
        header = hdu[0].header
        header["M0"] = refM0
        header["EXPTIME"] = 1.0
        hdu.flush()
    return GAINList, READOUTList, m0List, refM0



# parsing the argument line here
parser = argparse.ArgumentParser(
    description="Download fits-files of fields for specified coordinates")
parser.add_argument("filters", help="List of filters (for example gri or uiz)")
parser.add_argument("-i", "--input", default="coordinates.dat",
                    help="File with coordinates to download")
parser.add_argument("-a", "--adjacent", action="store_true", default=False,
                    help="Download adjacent fields if any exist")
parser.add_argument("-s", "--swarp", action="store_true", default=False,
                    help="Concatenate fields using SWarp package")
parser.add_argument("-c", "--convert", action="store_true", default=False,
                    help="Convert fields to the same photometry zeropoint")
parser.add_argument("-f", "--free", action="store_true", default=False,
                    help="Remove individual fields after concatenation")
parser.add_argument("-t", "--trim", action="store_true", default=False,
                    help="Crop image to galaxy size")
parser.add_argument("-p", "--ps", action="store_true", default=False,
                    help="Download psField files")
parser.add_argument("--scatter", action="store_true", default=False,
                    help="Put every object in a separate directory")
parser.add_argument("--add_urls", default=None,
                    help="File with additional urls of fields for objects")
parser.add_argument("--add_fields", default=None,
                    help="File wuth additional run,rerun,camcol,fields data for objects")
args = parser.parse_args()


bandlist = args.filters
# Make dirs for all bands and psf in case all files for the same 
# colour are in the same directories (scatter option is turned off)
if not args.scatter:
    for band in bandlist:
        if not os.path.exists("./downloads/%s" % (band)):
            os.makedirs("./downloads/%s" % (band))
    if args.ps:
        if not os.path.exists("./downloads/ps/"):
            os.makedirs("./downloads/ps/")
else:
    # if every object will be placed in the separate directory
    # (scatter option is on) create just the main downloads directory for now
    if not os.path.exists("./downloads/"):
        os.makedirs("./downloads/")


if args.swarp:
    # Check what name has SWarp package on this system
    rCode = subprocess.call("which swarp >/dev/null", shell=True)
    if rCode == 0:
        swarpName = "swarp"
    else:
        rCode = subprocess.call("which SWarp >/dev/null", shell=True)
        if rCode == 0:
            swarpName = "SWarp"
        else:
            print "\033[31m Error: SWarp was not found on your system.\033[0m"
            print "\033[31m The command has to be either 'swarp' or 'SWarp'\033[0m"
            print "\033[31m Intall SWarp package or try to run this script without -s option.\033[0m"
            exit(1)


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

        if (len(params) in (3, 4)) or ((len(params) == 4)
                                       and (args.adjacent is True)):
            galName = params[0]
            ra = float(params[1])
            dec = float(params[2])
            if args.adjacent is True:
                r_adj = float(params[3])
            else:
                r_adj = 0.0
            msg = "\033[33m Downloading field for "
            msg += "%1.5f %1.5f: '%s'  (%i/%i) \033[0m" % (ra, dec,
                                                           galName, counter,
                                                           len(listOfCoords))
            print msg
        else:
            print "Invalid number of columns in input file %s" % args.input
            sys.exit(1)
        objFieldList, objID = findField2(ra, dec, r_adj)
        if objID is None:
            print "\033[31m Error! No object was found at given coordinates.\033[0m"
            print "\033[31m This area is probably outside of the SDSS footprint.\033[0m"
            errFile.write("%s  %1.6f  %1.6f \n" % (galName, ra, dec))
            continue

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
                dth = Thread(target=download,
                             args=(urls[band], band, curGalName))
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
            msg = "         \033[34m [DONE] \033[0m    "
            msg += "(%1.2f sec)\n\n" % (time()-startTime)
            print msg
            outFile.write(" %s.fits " % (curGalName))
        outFile.write("\n")
        
        # If there are additional urls
        if (args.add_fields is not None) or (args.add_urls is not None):
            addNames = {band:[] for band in bandlist}
        if args.add_urls is not None:
            for line in open(args.add_urls):
                if line.split()[0] == galName:
                    listOfAddUrls = line.split()[1:]
                    thereAreAddFields = True
                    for band in bandlist:
                        for url in listOfAddUrls:
                            urlBand = url.replace("*", band)
                            outNameAdd = "./downloads/%s_%s" % (galName, urlBand.split("/")[-1])
                            addNames[band].append(outNameAdd[:-4])
                            print "Downloading additional field %s" % (outNameAdd)
                            subprocess.call("wget -nv -O %s %s" % (outNameAdd, urlBand), shell=True)
                    break
        if args.add_fields is not None:
            for line in open(args.add_fields):
                if line.split(":")[0] == galName:
                    listOfAddFields = line.split(":")[1]
                    for runData in listOfAddFields.split(";"):
                        run = int(runData.split()[0])
                        rerun = int(runData.split()[1])
                        camcol = int(runData.split()[2])
                        field = int(runData.split()[3])
                        for band in bandlist:
                            url = getUrl(run, rerun, camcol, field, band)
                            outNameAdd = "./downloads/%s_%s" % (galName, url.split("/")[-1])
                            addNames[band].append(outNameAdd[:-4])
                            print "Downloading additional field %s" % (outNameAdd)
                            subprocess.call("wget -nv -O %s %s" % (outNameAdd, url), shell=True)
            

        # Concatenating fields
        if args.swarp and (len(objFieldList) > 1):
            for band in bandlist:
                listOfImages = ["./downloads/%s_%i_%s.fits" % (galName, i, band) for i in xrange(len(objFieldList))]
                if (args.add_urls is not None) or (args.add_fields is not None):
                    listOfImages.extend(addNames[band])
                if args.convert:
                    GAINList, READOUTList, m0List, refM0 = reduce_to_same_m0(listOfImages)
                print "Running SWarp for %s band..." % (band)
                callSt = "%s -verbose_type quiet -BACK_TYPE MANUAL " % (swarpName)
                callSt += " ".join(["%s[0]" % (s) for s in listOfImages])
                subprocess.call(callSt, shell="True")
                move("./coadd.fits", "./downloads/%s_%s.fits" % (galName, band))
                os.remove("coadd.weight.fits")
                os.remove("swarp.xml")
                if args.free:
                    for fN in listOfImages:
                        os.remove(fN)
                # store mean keywords to coadded file
                if args.convert:
                    hdu = pyfits.open("./downloads/%s_%s.fits" % (galName, band),
                                      do_not_scale_image_data=True,
                                      mode="update")
                    header = hdu[0].header
                    header["GAIN"] = np.mean(GAINList)
                    header["READOUT"] = np.mean(READOUTList)
                    header["M0"] = refM0
                    header["EXPTIME"] = 1.0
                    hdu.flush()

        # Convert singlefield images
        if args.convert and (len(objFieldList) == 1):
            for band in bandlist:
                fName = "./downloads/%s_%s.fits" % (galName, band)
                bunzip("%s.bz2" % (fName))
                prep_ima(fName)
                GAIN, READOUT, m0 = SDSS_dr8(fName)
                hdu = pyfits.open(fName, do_not_scale_image_data=True,
                                  mode="update")
                header = hdu[0].header
                header["M0"] = m0
                header["EXPTIME"] = 1.0
                hdu.flush()

        # Crop images
        if args.trim and (not wcsOK):
            print "Astropy module was not found. Images cannot be trimmed"
        elif args.trim and wcsOK:
            print "Cropping..."
            # At first determine the common size of cropping area
            # (so all images will be of the same size)
            pixelCoords = {}
            cropSizes = []
            for band in bandlist:
                fName = "./downloads/%s_%s.fits" % (galName, band)
                hdu = pyfits.open(fName)
                data = hdu[0].data
                wcs = WCS(fName)
                xGalPix, yGalPix = wcs.wcs_world2pix([[ra, dec]], 1)[0]
                size = min(r_adj*60.0/0.396127, xGalPix, yGalPix,
                           data.shape[1]-xGalPix, data.shape[0]-yGalPix)
                hdu.close()
                pixelCoords[band] = (int(xGalPix), int(yGalPix))
                cropSizes.append(size)
            commonCropSize = int(min(cropSizes))
            for band in bandlist:            
                fName = "./downloads/%s_%s.fits" % (galName, band)
                hdu = pyfits.open(fName)
                data = hdu[0].data
                header = hdu[0].header
                xCropMin = pixelCoords[band][0] - commonCropSize
                xCropMax = pixelCoords[band][0] + commonCropSize
                yCropMin = pixelCoords[band][1] - commonCropSize
                yCropMax = pixelCoords[band][1] + commonCropSize
                data = data[yCropMin:yCropMax, xCropMin:xCropMax]
                outHDU = pyfits.PrimaryHDU(data=data, header=header)
                fOutName = "./downloads/%s_%s_trim.fits" % (galName, band)
                if os.path.exists(fOutName):
                    os.remove(fOutName)
                outHDU.writeto(fOutName)
                hdu.close()
                if args.free:
                    os.remove(fName)

        # Downloading and processing are finished. Now we have to place files to folders
        # according to scatter option
        if args.scatter:
            # every object has separate filder, where all files related to this object
            # are located
            fileList = glob.glob("./downloads/%s*.fits" % (galName))
            dst = "./downloads/%s/" % (galName)
            if not os.path.exists(dst):
                os.mkdir(dst)
            for src in fileList:
                move(src, dst)
        else:
            # scatter option is off, so all images taken in the same passband
            # will be in the same folder.
            for band in bandlist:
                fileList = glob.glob("./downloads/%s_*%s.fits" % (galName, band))
                if args.trim:
                    fileList.append("./downloads/%s_%s_trim.fits" % (galName, band))
                dst = "./downloads/%s/" % (band)
                for src in fileList:
                    if os.path.exists(src):
                        move(src, dst)
            if args.ps:
                # move psFields
                fileList = glob.glob("./downloads/%s_*ps.fits" % (galName))
                dst = "./downloads/ps/"
                for src in fileList:
                    if os.path.exists(src):
                        move(src, dst)
