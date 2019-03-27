import numpy as np
import optparse

def siftSN(CANDFILE):
    cand = np.genfromtxt(CANDFILE)
    filtSN = cand[(cand[:,0] >= 10.0)]
    return filtSN

def siftDM(CANDFILE):
    cand = np.genfromtxt(CANDFILE)
    filtDM = cand[(cand[:,5] >= 25.0)]
    return filtDM

def siftMASK(CANDFILE):
    cand = np.genfromtxt(CANDFILE)
    filtMASK = cand[(cand[:,10] <= 65536)]
    return filtMASK

def siftWIDTH(CANDFILE):
    cand = np.genfromtxt(CANDFILE)
    filtWIDTH = cand[(cand[:,3] < 10)]
    return filtWIDTH

def siftBEAMS(CANDFILE):
    cand = np.genfromtxt(CANDFILE)
    filtBEAMS = cand[(cand[:,9] < 4)]
    return filtBEAMS

def siftall(CANDFILE):
    cand = np.genfromtxt(CANDFILE)
    filtcand = cand[(cand[:,0] >= 10.0) & (cand[:,5] >= 25.0) & (cand[:,10] <= 65536) & (cand[:,3] < 10) & (cand[:,9] < 4)]
    return filtcand

if __name__=="__main__":
    desc = "bla"
    parser = optparse.OptionParser(description=desc)
    parser.add_option('--candfile',dest='candfile',type='string')
    (opts,args) = parser.parse_args()
    cand = np.genfromtxt(opts.candfile)
    print "candsize {}".format(np.shape(cand))
    print "{}".format(np.shape(cand)[0])
    cands = siftSN(opts.candfile)
    print "sift SN size {}".format(np.shape(cands))
    print "{}".format(np.shape(cands)[0])
    cands = siftDM(opts.candfile)
    print "sift DM size {}".format(np.shape(cands))
    print "{}".format(np.shape(cands)[0])
    cands = siftMASK(opts.candfile)
    print "sift MASK size {}".format(np.shape(cands))
    print "{}".format(np.shape(cands)[0])
    cands = siftWIDTH(opts.candfile)
    print "sift WIDTH size {}".format(np.shape(cands))
    print "{}".format(np.shape(cands)[0])
    cands = siftBEAMS(opts.candfile)
    print "sift BEAMS size {}".format(np.shape(cands))
    print "{}".format(np.shape(cands)[0])
    cands = siftall(opts.candfile)
    print "sift all size {}".format(np.shape(cands))
    print "{}".format(np.shape(cands)[0])
