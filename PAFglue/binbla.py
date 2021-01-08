import numpy as np

beams = np.loadtxt('beams.txt')
x = beams[:,0]; y = beams[:,1]



def check_mask(x,y,BEAM,CANDMASK):
    origx = x[BEAM]; origy = y[BEAM]
    #dist = np.sqrt((origx-x[0:22])**2 + (origy-y[0:22])**2)
    dist = np.sqrt((origx-x[0:32])**2 + (origy-y[0:32])**2)
    good_beams = np.asarray(np.where(dist < .13)).flatten()
    bad_beams = np.asarray(np.where(dist > .13)).flatten()
    #good_mask = np.binary_repr(np.sum(2**good_beams),width=22)
    #bad_mask = np.binary_repr(np.sum(2**bad_beams),width=22)
    #cmask = np.binary_repr(CANDMASK,width=22)
    good_mask = np.binary_repr(np.sum(2**good_beams),width=32)
    bad_mask = np.binary_repr(np.sum(2**bad_beams),width=32)
    cmask = np.binary_repr(CANDMASK,width=32)
    #print "beam {}".format(BEAM)
    #print dist
    #print "cmask {}".format(cmask)
    #print "gmask {}".format(good_mask)
    #print "bmask {}".format(bad_mask)
    bm = int(bad_mask,2)
    cm = int(cmask,2)
    andy = bin(bm & cm) 
    #print "candmask and bad mask {}".format(andy)
    if andy == bin(0):
        return True
        #print "and is zero, this is a good boy"
    else:
        return False
        #print "candmask ands with badmask, booo"
        #print "cmask {}".format(cmask)
        #print "gmask {}".format(good_mask)
        #print "bmask {}".format(bad_mask)
        #print "candmask and bad mask {}".format(andy)
    #print "--------------------"

if __name__=="__main__":
    cands = np.loadtxt('2018-08-29-05:02:29_all.cand')
    numcands = len(cands[:,0])
    #tstcand = cands[2512,:]
    #beam = int(tstcand[-1]) - 1
    #candmask = int(tstcand[10])
    #check_mask(x,y,beam,candmask)

    for i in range(numcands):
        print "cand no {}".format(i)
        tstcand = cands[i,:]
        beam = int(tstcand[-1]) - 1
        candmask = int(tstcand[10])
        check_mask(x,y,beam,candmask)
'''
for xy in zip(x[0:22],y[0:22]):
    origx = x[beam]; origy = y[beam]
    dist = np.sqrt((origx-x[0:22])**2 + (origy-y[0:22])**2)
    good_beams = np.asarray(np.where(dist < .13)).flatten()
    print dist
    bad_beams = np.asarray(np.where(dist > .13)).flatten()
    good_mask = np.binary_repr(np.sum(2**good_beams),width=22)
    bad_mask = np.binary_repr(np.sum(2**bad_beams),width=22)
    gm = int(good_mask,2)
    bm = int(bad_mask,2)
    xor = bin(gm ^ bm)
    andy = bin(gm & bm)
    print "beam {}".format(beam)
    print "xor {}".format(xor)
    print "len xor {}".format(len(xor))
    print "and {}".format(andy)
    if andy == bin(0):
        print "and is zero"
    print "good beams {}".format(good_beams)
    print "good beam mask {}".format(good_mask)
    print "len {}".format(len(good_mask))
    print "bad beams {}".format(bad_beams)
    print "bad beam mask {}".format(bad_mask)
    print "len {}".format(len(bad_mask))
    beam +=1
'''
