import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import optparse
import glob
import subprocess
import os


def header(afile):
    inread = ""
    while True:
        tmp = afile.read(1)
        inread = inread + tmp
        flag = inread.find('HEADER_END')
        if flag != -1:
            break
    return inread

def grab_data(FILE,STARTSAMP,NUMSAMP,NCHAN,DTYPE,FIL):
    with open(FILE,'r') as F:
        if not FIL:
            headlen = 4096
        else:
            thehead = header(F)
            headlen = len(thehead)
        F.seek(headlen + NCHAN * STARTSAMP)
        data = np.fromfile(F,dtype=DTYPE,count=int(NCHAN*NUMSAMP))
        data = np.reshape(data,(-1,NCHAN)).T
    return data

def plotshit(DATA,NCHAN,SAVE,TIME):
    fig = plt.figure(1)
    ax = fig.add_subplot(1,1,1)
    ax.set_title("cmap")
    ax.set_ylabel("channel")
    ax.set_xlabel("sample")
    cax = ax.imshow(DATA,aspect='auto',interpolation='none')
    if SAVE:
        name = os.path.basename(os.path.splitext(opts.data)[0]) + "_time" + str(TIME) + "_cmap.png"
        outname = os.path.join(opts.out,name)
        print "name {}".format(name)
        plt.savefig(outname,format='png',dpi=200)

    STDs = np.std(DATA,axis=1)
    fig2 = plt.figure(2)
    ax = fig2.add_subplot(1,1,1)
    ax.set_title("median over channels")
    ax.set_ylabel("median")
    ax.set_xlabel("channel")
    MEDIANS = np.median(DATA,axis=1)
    print "len of medians {}".format(len(MEDIANS))
    # 90th percentile
    percs = np.zeros(NCHAN)
    for i in range(NCHAN):
        percs[i] = np.percentile(DATA[i,:],90)
    ax.plot(np.arange(NCHAN),MEDIANS)
    ax.fill_between(np.arange(NCHAN),MEDIANS+STDs,MEDIANS-STDs,facecolor='red',alpha=.5)
    ax.scatter(np.arange(NCHAN),percs,marker='.',color='magenta',s=1)
    if SAVE:
        name = os.path.basename(os.path.splitext(opts.data)[0]) + "_time" + str(TIME) +  "_std.png"
        outname = os.path.join(opts.out,name)
        print "name {}".format(name)
        plt.savefig(outname,format='png',dpi=200)
    #fig.show()
    plt.close('all')

if __name__=='__main__':
    desc = "I'm a descriptive description"
    parser = optparse.OptionParser(description=desc)
    parser.add_option('--nchan',dest='nchan',type='int',default=512)
    parser.add_option('--data',dest='data',type='str')
    parser.add_option('--samp',dest='samp',type='int',default=5000)
    parser.add_option('--dtype',dest='dtype',type='str',default=np.uint8)
    parser.add_option('--samptime',dest='samptime',type='float',default=216e-06)
    parser.add_option('--pad',dest='pad',type='int',default=900)
    parser.add_option('--out',dest='out',type='str',default='')
    (opts,args) = parser.parse_args()
    startsamp = opts.samp - opts.pad
    timeint = round(opts.samp * opts.samptime,2) # time for name of plot
    isfil = True
    ext = os.path.splitext(os.path.basename(opts.data))[-1]
    print "ext is {}".format(ext)
    if ext=='.dada':
        isfil = False
        print "im here"
        print "isfil {}".format(isfil)
    save = False
    if opts.out:
        save = True
        if not os.path.exists(opts.out):
            try:
                subprocess.check_call(["mkdir",opts.out])
            except OSError as error:
                print error
    data = grab_data(opts.data,startsamp,2*opts.pad,opts.nchan,opts.dtype,isfil) 
    plotshit(data,opts.nchan,save,timeint)
    
