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
    filt = cand[(cand[:,3] < 8)] # PAFMAR19 equates to 22 ms or less
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
    filtcand = cand[(cand[:,0] >= 8.0) & (cand[:,5] >= 20) & (cand[:,3] < 8)]
    return filtcand

if __name__=="__main__":
    desc = "bla, I describe things"
    parser = optparse.OptionParser(description=desc)
    parser.add_option('--candfile',dest='candfile',type='string')
    parser.add_option('--datadir',dest='datadir',type='string')
    (opts,args) = parser.parse_args()
    dir1 = os.path.dirname(opts.candfile)
    plotpath = os.path.join(dir1,"plots")
    numcands = origlen(opts.candfile)
    sncands = siftSN(opts.candfile)
    dmcands = siftDM(opts.candfile)
    wcands = siftW(opts.candfile) 
    maskcands = siftMASK(opts.candfile)
    #print "bla {}".format(numcands)
    #print "bla {}".format(sncands)
    #print "bla {}".format(dmcands)
    #print "bla {}".format(wcands)
    #print "bla {}".format(maskcands)
    candidates = siftcands(opts.candfile)
    totplotted = 0
    fils = glob.glob(os.path.join(opts.datadir,"*.fil"))
    gcands = os.path.join(os.path.dirname(opts.candfile),'zcands.txt')
    with open(gcands,'w') as f:
        f.write("#S/N\tsampidx\ttoff\t\tfiltidx\tDMidx\tDM\tassoc\tearly\tlate\tnbeams\tmask\tmaxbeam\tmaxSN\tbeam\n") 
    for cand in candidates:
        if cand[0] != cand[12]:
            print "this event is stronger in another beam, skip"
            continue
        beams = np.loadtxt('beams.txt')
        x = beams[:,0]; y = beams[:,1]
        B = int(cand[-1]) - 1
        CM = int(cand[10])
        if check_mask(x,y,B,CM):
            totplotted += 1
            for fil in fils:
                findme = 'beam{}.fil'.format("%02d"%B)
                #findme = 'beam_{}_'.format(B+1)
                name = os.path.basename(fil)
                if name.find(findme) > 0:
                    thefil = fil
                    break
            #print "beam {}".format(B)
            #print "fil {}".format(thefil)
            #thefil = '/beegfs/some/test/path/meow_beam{}.fil'.format("%02d"%B)
            with open(gcands,'a') as f:
                #f.write(str(cand))
                f.write("{:.2f}\t{:d}\t{:.5f}\t{:d}\t{:d}\t{:.2f}\t{:d}\t{:d}\t{:d}\t{:d}\t{:d}\t{:d}\t{:.2f}\t{:d}".format(cand[0],int(cand[1]),cand[2],int(cand[3]),int(cand[4]),cand[5],int(cand[6]),int(cand[7]),int(cand[8]),int(cand[9]),int(cand[10]),int(cand[11]),cand[12],int(cand[13])))
                f.write("\n")
            DM = str(cand[5])
            SAMP = str(int(cand[1]))
            W = str(int(2**cand[3]))
            BEAM = str(B)
            subprocess.check_call(["python","generalplotter.py","--ftop","1459","--fchan","0.44907","--nchan","512","--samptime","216e-03","--data",thefil,"--out",plotpath,"--dm",DM,"--samp",SAMP,"-w",W,"--beamno",BEAM,"--tscrunch","2"])
    #print "bla {}".format(totplotted)
    statout = os.path.join(os.path.dirname(opts.candfile),'stats.txt')
    with open(statout,'w') as f:
        f.write("#{}".format(os.path.basename(opts.candfile)))
        f.write("\n")
        f.write("total cands: {}\n".format(numcands))
        f.write("SN cands: {}\n".format(sncands))
        f.write("DM cands: {}\n".format(dmcands))
        f.write("Width cands: {}\n".format(wcands))
        f.write("Mask cands: {}\n".format(maskcands))
        f.write("Total remaining: {}\n".format(totplotted))


    



