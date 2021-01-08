import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import optparse
import glob
import subprocess
import os
from binbla import check_mask

def origlen(CANDFILE):
    cand = np.genfromtxt(CANDFILE)
    return len(cand[:,0])

def siftSN(CANDFILE):
    cand = np.genfromtxt(CANDFILE)
    filt = cand[(cand[:,0] >= 8.0)]
    return len(filt)

def siftDM(CANDFILE):
    cand = np.genfromtxt(CANDFILE)
    filt = cand[(cand[:,5] >= 20.0)]
    return len(filt)

def siftW(CANDFILE):
    cand = np.genfromtxt(CANDFILE)
    filt = cand[(cand[:,3] < 9)] # width of 55 ms or less
    return len(filt)

def siftMASK(CANDFILE):
    l = 0
    cand = np.genfromtxt(CANDFILE)
    lencands = len(cand[:,0])
    beams = np.loadtxt('beams.txt')
    x = beams[:,0]; y = beams[:,1]
    for i in range(lencands):
        B = int(cand[i,-1]) - 1
        CM = int(cand[i,10]) 
        the_truth = check_mask(x,y,B,CM)
        if the_truth:
            l += 1
    return l

def siftcands(CANDFILE):
    cand = np.genfromtxt(CANDFILE)
    filtcand = cand[(cand[:,0] >= 8.0) & (cand[:,5] >= 20) & (cand[:,3] < 9)]
    return filtcand

def find_phil(PATH):
    phils = []
    for dirpath,dirnames,filenames in os.walk(PATH):
        for f in filenames:
            ff = os.path.join(dirpath,f)
            if os.path.splitext(ff)[-1] == '.fil':
                phils.append(ff)
    return phils

def find_cands(PATH):
    c = []
    for dirpath,dirnames,filenames in os.walk(PATH):
        for f in filenames:
            ff = os.path.join(dirpath,f)
            if ff.endswith('all.cand'):
                c.append(ff)
    return c

if __name__=="__main__":
    desc = "bla, I describe things"
    parser = optparse.OptionParser(description=desc)
    parser.add_option('--source',dest='source',type='string')
    parser.add_option('--search',dest='search',type='string',
        default='/beegfsEDD/PAF/PAF/SEARCH')
    parser.add_option('--results',dest='results',type='string',
        default='/beegfsEDD/PAF/PAF/RESULTS')
    (opts,args) = parser.parse_args()
    bigdir = opts.results
    otherbigdir = opts.search
    Source = glob.glob(os.path.join(bigdir,opts.source+"*")) # find results dirs for the source
    for d in Source: # iterate through results directories
        D = find_cands(d)[0] # candidate file
        outdir = os.path.join(os.path.dirname(D),"CANDS")
        if not os.path.exists(outdir):
            try:
                subprocess.check_call(["mkdir",outdir])
            except OSError as error:
                print error
        thename = D.split(os.sep)[5]
        t = thename.split('_')[-1]
        tt = "{}-{}-{}T{}:{}:{}".format(t[0:4],t[4:6],t[6:8],t[9:11],t[11:13],t[13:15])
        s = thename.split('_')[0]
        tmpdir = os.path.join(otherbigdir,t+"*"+s) # filterbank dir
        fildir = glob.glob(tmpdir)[0]
        fils = find_phil(fildir)
        numcands = origlen(D)
        sncands = siftSN(D)
        dmcands = siftDM(D)
        wcands = siftW(D) 
        maskcands = siftMASK(D)
        candidates = siftcands(D)
        totplotted = 0
        gcands = os.path.join(outdir,'zcands.txt')
        # create good cand file
        with open(gcands,'w') as f:
            f.write("#S/N\tsampidx\ttoff\t\tfiltidx\tDMidx\tDM\tassoc\tearly\tlate\tnbeams\tmask\tmaxbeam\tmaxSN\tbeam\n") 
        for cand in candidates:
            if cand[0] != cand[12]:
                print "this event is stronger in another beam, skip"
                continue
            beams = np.loadtxt('beams.txt')
            x = beams[:,0]; y = beams[:,1]
            B = int(cand[-1]) 
            CM = int(cand[10])
            if check_mask(x,y,B,CM):
                totplotted += 1
                for fil in fils:
                    findme = "BEAM_{:0>2d}.fil".format(B) 
                    name = os.path.basename(fil)
                    if name.find(findme) > 0:
                        thefil = fil
                        break
                with open(gcands,'a') as f:
                    f.write("{:.2f}\t{:d}\t{:.5f}\t{:d}\t{:d}\t{:.2f}\t{:d}\t{:d}\t{:d}\t{:d}\t{:d}\t{:d}\t{:.2f}\t{:d}".format(cand[0],int(cand[1]),cand[2],int(cand[3]),int(cand[4]),cand[5],int(cand[6]),int(cand[7]),int(cand[8]),int(cand[9]),int(cand[10]),int(cand[11]),cand[12],int(cand[13])))
                    f.write("\n")
                DM = str(cand[5])
                SAMP = str(int(cand[1]))
                W = str(int(2**cand[3]))
                BEAM = str(B)
                try:
                    subprocess.check_call(["python","candpanel.py",
                        "--PAF",
                        "--data",thefil,
                        "--downsamp",str(4), # feels like a good downsamp factor for the PAF
                        "--out",outdir,
                        "--dm",DM,
                        "--samp",SAMP,
                        "--beamno",BEAM])
                except subprocess>CalledProcessError as err:
                    print err
        statout = os.path.join(outdir,'stats.txt')
        with open(statout,'w') as f:
            f.write("#{}".format(os.path.basename(D)))
            f.write("\n")
            f.write("total cands: {}\n".format(numcands))
            f.write("SN cands: {}\n".format(sncands))
            f.write("DM cands: {}\n".format(dmcands))
            f.write("Width cands: {}\n".format(wcands))
            f.write("Mask cands: {}\n".format(maskcands))
            f.write("Total remaining: {}\n".format(totplotted))
