import numpy as np
import sys
import glob
import optparse
import os
import subprocess
import struct

def header(afile):
    inread = ""
    while True:
        tmp = afile.read(1)
        inread = inread + tmp
        flag = inread.find('HEADER_END')
        if flag != -1:
            break
    return inread

def get_headparam(head, parlist):
    how2parse={
        'nchans': ('i',4),
        'tsamp': ('d',8),
        'foff': ('d',8),
        'fch1': ('d',8),
        'tstart': ('d',8),
        'ibeam': ('i',4)
    }
    n = 0
    for i in parlist:
        i1 = head.find(i)
        i2 = i1 + len(i)
        nbytes = how2parse[i][1]
        cstr = how2parse[i][0]

        val = struct.unpack(cstr, head[i2:i2+nbytes])[0]
        parlist[n] = val
        n += 1
        return parlist

def update_headparam(head,parlist,vallist):
    how2parse={
        'nchans': ('i',4),
        'tsamp': ('d',8),
        'foff': ('d',8),
        'fch1': ('d',8),
        'tstart': ('d',8),
        'ibeam': ('i',4)
        }
    n = 0
    for i in parlist:
        i1 = head.find(i)
        i2 = i1 + len(i)
        nbytes = how2parse[i][1]
        cstr = how2parse[i][0]

        val = vallist.pop(0)
        val = struct.pack(cstr,val)
        head = head[0:i2] + val + head[i2+nbytes:]
        n += 1
    return head

def change_tsamp(head,downsamp):
    old_ds = get_headparam(head,['tsamp'])[0]
    new_ds = old_ds * downsamp
    #print "old tsamp {}".format(old_ds)
    #print "new tsamp {}".format(new_ds)
    #print "len1 {}".format(len(head))
    head = update_headparam(head,['tsamp'],[new_ds])
    #print "len2 {}".format(len(head))
    return head

def read_n_ds(fil,downsamp,outfile,blocksize):
    with open(fil,'r') as X:
        print "working on file {}".format(fil)
        thehead = header(X)
        #print "header {}".format(thehead)
        #print "lenhead {}".format(len(thehead))
        nchan = get_headparam(thehead,['nchans'])[0]
        headlen = len(thehead)
        X.seek(headlen)
        head = change_tsamp(thehead,downsamp)
        with open(outfile,'w') as f:
            print "writing header to new file {}".format(outfile)
            f.write(head)
        COUNTER = 1
        datacounter = 0
        datacounter2 = 0
        #print "newheader {}".format(head)
        #print "newlenhead {}".format(len(head))
        print "downsampling"
        while 1:
            #print "downsampling chunk number {}".format(COUNTER)
            dataIn = np.fromfile(X,dtype=np.uint8,count=nchan*blocksize*downsamp)
            if len(dataIn) == 0: break
            datacounter += len(dataIn)
            #print "downsampling chunk"
            dataIn = downsamp_data(dataIn,downsamp,nchan)
            datacounter2 += len(dataIn)
            #print "downsampling chunk done"
            with open(outfile,'a') as f:
                #print "writing chunk to file {}".format(outfile)
                f.write(dataIn)
            COUNTER += 1
            #print datacounter/4096
            #print datacounter2/4096
        print "downsampling done"

def downsamp_data(data,downsamp,nchan):
    lendata = len(data)
    #print "lendata chunk {}".format(lendata)
    data = np.reshape(data, (-1,nchan)).T
    #print "shape of data {}".format(np.shape(data))
    #print "len of [:,0] {}".format(len(data[:,0]))
    #print "len of [0,:] {}".format(len(data[0,:]))
    #print "SHAPE 1 {}".format(np.shape(data))
    #print "NUMBER 1 {}".format(len(data[0,:]))
    newsamps = int(np.floor(len(data[0,:]) / downsamp))
    #print "new number of samps in chunk {}".format(newsamps)
    newdata = np.zeros((nchan,newsamps))
    #print "shape of newdata {}".format(np.shape(newdata))
    #print "downsampling actually happening"
    for tavg in np.arange(newsamps):
        subdata = data[:, tavg * downsamp : (tavg + 1) * downsamp]
        subavg = np.reshape(np.sum(subdata,axis=1), (nchan,1)) / downsamp
        newdata[:,tavg] = subavg[:,0]
    #print "newdata new shape {}".format(np.shape(newdata))
    #print "actual downsampling done"
    newdata = newdata.T
    #print "newdata final shape {}".format(np.shape(newdata))
    newdata = newdata.flatten()
    newdata = newdata.astype(np.uint8)
    return newdata
    # SOMEHOW RETURN THE DATA IN ORIGINAL FORM, I.E. FLATTEN SHIT

def test():
    thefile = 'FRB121102_180420_4sec.fil'
    outfile = 'dstest.fil'
    downsamp = 4
    blocksize = 100
    read_n_ds(thefile,downsamp,outfile,blocksize)

if __name__=="__main__":
    desc = "downsample filterbanks"
    parser = optparse.OptionParser(description=desc)
    parser.add_option('--dir',dest='d',type='string',
        help="directory of filterbanks")
    parser.add_option('--ds',dest='ds',type='int',default=4,
        help="downsample factor. Default 4")
    parser.add_option('--blocksize',dest='bs',type='int',default=100,
        help="blocksize. Default 100")
    (opts,args) = parser.parse_args()
    outdir = os.path.join(opts.d,'downsampled')
    try:
        os.mkdir(outdir)
    except OSError as error:
        print error
    thefiles = os.path.join(opts.d,"I*.fil")
    files = sorted(glob.glob(thefiles))
    #print "outdir {}".format(outdir)
    #print "downsample factor {}".format(opts.ds)
    #print "blocksize {}".format(opts.bs)
    for FIL in files:
        tmp = os.path.splitext(os.path.basename(FIL))
        outfile = os.path.join(outdir,tmp[0] + "_ds{}".format(opts.ds) + tmp[1])
        read_n_ds(FIL,opts.ds,outfile,opts.bs)
        #print "current fil {}".format(FIL)
        #print "outfile {}".format(outfile)
    #test()

###
# read all filterbanks in directory
# loop through each filterbank
# read header
# change tsamp in header to tsamp*DOWNSAMP
# read in data in chunks
# downsamp chunk
# write chunk to a new file
# done
###
