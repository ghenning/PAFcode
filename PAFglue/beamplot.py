import numpy as np
import matplotlib.pyplot as plt

a = np.loadtxt('beams.txt')

x = a[:,0]
y = a[:,1]

fig = plt.figure()
ax = fig.add_subplot(1,1,1)
plt.plot(x[0:22],y[0:22],'.')
#plt.plot(x,y,'.')
b = 0
for xy in zip(x[0:22],y[0:22]):
#for xy in zip(x,y):
    origx = x[b]; origy = y[b]
    dist = np.sqrt((origx-x[0:22])**2 + (origy-y[0:22])**2)
    print b
    print dist
    good_beams = np.where(dist < .13)
    bad_beams = np.where(dist > .13)
    good_beams = np.asarray(good_beams).flatten()
    bad_beams = np.asarray(bad_beams).flatten()
    good_mask = np.binary_repr(np.sum(2**(good_beams)),width=22) 
    bad_mask = np.binary_repr(np.sum(2**np.asarray(bad_beams)),width=22)
    gm = int(good_mask,2)
    bm = int(bad_mask,2)
    xor = bin(gm ^ bm)
    andy = bin(gm & bm)
    print "xor {}".format(xor)
    print "len xor {}".format(len(xor))
    print "and {}".format(andy)
    if andy == bin(0):
        print "andy is zero"
    print "np shape {}".format(np.shape(good_beams))
    print "good beams {}".format(good_beams)
    print "good beam mask {}".format(good_mask)
    print "len {}".format(len(good_mask))
    print "bad beams {}".format(bad_beams)
    print "bad beam mask {}".format(bad_mask)
    print "len {}".format(len(bad_mask))
    #print "xor {}".format(xor)
    ax.annotate('%s' % b, xy=xy)
    b += 1
    #print b
    #ax.annotate('(%s,%s)' % xy, xy=xy,textcoords='data')
plt.show()
