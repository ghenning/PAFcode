import numpy as np
import sys
import glob
import optparse
import os
from astropy.time import Time
from datetime import datetime, timedelta
import struct
import subprocess
from welford import Welford # needs welford.py in the same dir to calculate running statistics

# reads the header from the file and returns it
def header(afile):
    inread = ""
    while True:
        tmp = afile.read(1)
        inread = inread + tmp
        flag = inread.find('HEADER_END')
        if flag != -1:
            break
    return inread

# reads output from header, and returns header values you ask for
def get_headparam(head, parlist):
    how2parse={
        'nchans': ('i',4),
        'tsamp': ('d',8),
        'foff': ('d',8),
        'fch1': ('d',8),
        'tstart': ('d',8)
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

# reads output from header, and changes what you want
def update_headparam(head, parlist, vallist):
    how2parse = {
        'nchans': ('i',4),
        'tsamp': ('d',8),
        'foff': ('d',8),
        'fch1': ('d',8),
        'tstart': ('d',8)
        }
    n = 0
    for i in parlist:
        i1 = head.find(i)
        i2 = i1 + len(i)
        nbytes = how2parse[i][1]
        cstr = how2parse[i][0]

        val = vallist.pop(0)
        val = struct.pack(cstr, val)
        head = head[0:i2] + val + head[i2+nbytes:]
        n += 1
    return head

# finds the files you want
def find_files(directory,beamy):
    path = os.path.join(directory,"I*beam_" + str(beamy) + "_*.fil")
    files = sorted(glob.glob(path))
    return files

# just some test stuff, ignore
def printing_n_counting(DIR):
    files = find_files(DIR)
    for f in files:
        infile = open(f,'r')
        head = header(infile)
        headlen = len(head)
        print "HEADER LENGTH"
        print headlen
        print "Here's the HEADER"
        print head
        infile.close()

# creates the output file, which everything will be written to
def create_new_file(fil,outfile,BLA,starty):
    print " opening %s to make header of new file" % fil
    with open(fil,'r') as f:
        head = header(f)
    if BLA:
        head = update_headparam(head, ['tstart'], [starty])
    with open(outfile,'w') as f:
        f.write(head)
    print "header written to new file"

# garbage
#'''def add_to_new(fil,outfile):
#    blocksize = 5000 # make this optional?
#    with open(fil,'r') as X: 
#        headlen = len(header(X)) 
#        X.seek(headlen)
#        while 1:
#            dataIn = X.read(512*blocksize)
#            if len(dataIn) == 0: break
#            with open(outfile,'a') as f:
#                f.write(dataIn)'''

# loops through filterbanks and adds to the output files
def add_to_new(fils,counter,endMJD,outfile,logfile):
    RunStat = Welford()
    blocksize = 5000 # make this optional?
    # make nchan an option as well
    going = True
    beam, MJD = get_info(fils[counter])
    file_end_MJD = MJD
    biglen = 0
    #print "starting counter is %s" % counter
    nchans = 512 # make as an option, or read from header or something
    samptime = 54.e-6 # I think just setting this to zero is fine, because it grabs that info later
    # it doesn't pad on the first file, so it should be fine
    while going:
        with open(logfile,'a') as L:
            L.write("#Working on file %s\n" % fils[counter])
        # find mjd, compare to end MJD, if OK continue
        beam, MJD = get_info(fils[counter])
        print "counter is at %s " % counter
        print "the starting MJD is %s" % MJD
        print "checking if MJD of file %s fits" % fils[counter]
        if (endMJD - MJD) <= 0:
            print "out of MJD range, kill me"
            #going = False
            return biglen, samptime, RunStat.mean, RunStat.std
            #break
        print "MJD in range, let's do thiiiiiis!"
        MJDdiff = MJD - file_end_MJD
        #print MJDdiff
        print "difference between end of last file and starting of current: %s" % MJDdiff
        print "in seconds: %s" % (MJDdiff*86400)
        if MJDdiff >= samptime: # this needs to be less than some value, e.g. sample time
            print "there's a difference between start/end times, need to do some padding"
            numsamp_2_add = round(MJDdiff*86400/samptime)
            print "samples to add: %s " % numsamp_2_add
            print "padding with 512x%s of zeros" % round(numsamp_2_add)
            with open(logfile,'a') as L:
                timeintofile = (biglen/512) * samptime 
                L.write("#Padding %s samples, or %s seconds\n" % (numsamp_2_add, numsamp_2_add*samptime))
                L.write("#%s seconds into the file\n" % timeintofile)
                L.write("%s\t%s\t%s\n" % (numsamp_2_add,numsamp_2_add*samptime,timeintofile))
            with open(outfile,'a') as f:
                padnum = numsamp_2_add * 512
                current_mean = RunStat.mean
                print "current mean: %s " % current_mean
                current_std = RunStat.std
                print "current std: %s " % current_std
                noisecounter = 0
                diff = 1
                noiseblock = 500
                while diff != 0:
                    adding = noiseblock * nchans
                    noisecounter += adding
                    gauss_noise = np.random.normal(current_mean,current_std,int(adding))
                    #print "$%(*$&%$(^#(^%(#@*^%#)(&%^)#&%^$)&%"
                    #print "len of gauss noise less than 0"
                    #print len(gauss_noise[gauss_noise < 0])
                    #print "$%(*$&%$(^#(^%(#@*^%#)(&%^)#&%^$)&%"
                    gauss_noise[gauss_noise < 0] = np.random.choice([0,round(current_mean)])
                    gauss_noise = gauss_noise.astype(np.uint8)
                    f.write(gauss_noise)
                    diff = padnum - noisecounter
                    if diff < adding:
                        noiseblock = diff/nchans
                        print "diff in padsamps and counter"
                        print diff
                        print "I should only appear twice per padding!"
                #gauss_noise = np.random.normal(current_mean,current_std,padnum).astype(np.uint8)
                #gauss_noise = np.random.normal(current_mean,current_std,padnum)
                #print "$%(*$&%$(^#(^%(#@*^%#)(&%^)#&%^$)&%"
                #print "len of gauss noise less than 0"
                #print len(gauss_noise[gauss_noise < 0])
                #print "$%(*$&%$(^#(^%(#@*^%#)(&%^)#&%^$)&%"
                #gauss_noise[gauss_noise < 0] = np.random.choice([0,round(current_mean)])
                #gauss_noise = gauss_noise.astype(np.uint8)
                #padding = np.zeros(padnum,dtype=np.uint8)
                #f.write(padding)
                #f.write(gauss_noise)
                biglen += padnum
            print "padding done, crossing fingers it works"
            #print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
            #print "modulus of samples needed to add and 512 (nchan)"
            #damod = round(numsamp_2_add) % 512
            #print damod
            #print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
            # make a def for padding stuff, call it here
        with open(fils[counter],'r') as X:
            thehead = header(X)
            headlen = len(thehead)
            samptime = get_headparam(thehead,['tsamp'])[0]
            print "sample time is %s" % samptime
            X.seek(headlen)
            print "adding %s to the juicy new file" % fils[counter]
            datacounter = 0
            while 1:
                #dataIn = X.read(512*blocksize)
                dataIn = np.fromfile(X,dtype=np.uint8,count=512*blocksize)
                if len(dataIn) == 0: break
                # here I need to do something with the data to calculate a running median/std
                RunStat(dataIn[len(dataIn)/2:len(dataIn)/2+512])
                with open(outfile,'a') as f:
                    # IF I do anything with the data, I need to make sure it has the correct
                    # datatype. So I will need the following line before writing
                    # dataIn = dataIn.astype(np.uint8)
                    f.write(dataIn)
                    biglen += len(dataIn)
                datacounter += len(dataIn)
        if fils[-1] == fils[counter]:
            print "no more files left, kill kill kill"
            #going = False
            return biglen, samptime, RunStat.mean, RunStat.std
            #break
        file_duration = (datacounter/512) * samptime 
        print "the duration of that file was %s seconds" % file_duration
        file_duration_MJD = file_duration/86400
        file_end_MJD = MJD + file_duration_MJD
        print "end time of file in MJD is %s " % file_end_MJD
        counter += 1


def get_info(fil):
    beam = int(fil.split("_")[-2])
    MJD = float(fil.split("_")[-4])
    return beam, MJD

def find_start_file(MJDobs,MJDend,filenames):
    #print "I'm the last one"
    #print filenames[-1]
    counter = 0
    for NAME in filenames:
        #print NAME
        beam, MJD = get_info(NAME)
        diffy = MJD - MJDobs 
        #print diffy
        if diffy >= 0 and (MJDobs <= MJD and MJD <= MJDend):
            print "hi"
            startdiffy = diffy *86400 # / 60
            print "difference in start time of filterbank vs input start time: %s seconds" % startdiffy
            return NAME, counter
        counter += 1
    #return "NO FILE IN RANGE"
    return 0,0

def find_end_file(MJDobs,filenames):
    oldNAME = ""
    # WHAT IF IT REACHES THE END OF THE OBS LIST?
    for NAME in filenames:
        beam, MJD = get_info(NAME)
        diffy = MJD - MJDobs
        if diffy >= 0:
            print "hi"
            print diffy * 86400 / 60
            return oldNAME
        oldNAME = NAME
    return "something is off"

def thepad(origfil,outfile,presamps,postsamps,beam_mean,beam_std,logfile):
    blocksize = 5000
    nchans = 512
    prepadnum = presamps * nchans
    postpadnum = postsamps * nchans
    pre = 0 # for future looping?
    post = 0 # for future looping?
    # for future noise addition/calc?
    #with open(origfil,'r') as X:
        #thehead = header(X)
        #headlen = len(thehead)
        #samptime = get_headparam(thehead,['tsamp'])[0]
    print "doing pre-padding with %s samples" % presamps
    if presamps > 0:
        with open(outfile,'a') as f: # will I need to loop over this with large paddings? or when noise is added
            current_mean = beam_mean
            current_std = beam_std
            noisecounter = 0
            diff = 1
            noiseblock = 500
            while diff != 0:
                adding = noiseblock * nchans
                noisecounter += adding
                gauss_noise = np.random.normal(current_mean,current_std,int(adding))
                #print "$%(*$&%$(^#(^%(#@*^%#)(&%^)#&%^$)&%"
                #print "len of gauss noise less than 0"
                #print len(gauss_noise[gauss_noise < 0])
                #print "$%(*$&%$(^#(^%(#@*^%#)(&%^)#&%^$)&%"
                gauss_noise[gauss_noise < 0] = np.random.choice([0,round(current_mean)])
                gauss_noise = gauss_noise.astype(np.uint8)
                f.write(gauss_noise)
                diff = prepadnum - noisecounter
                if diff < adding:
                    noiseblock = diff/nchans
                    print "diff in padsamps and counter"
                    print diff
                    print "I should only appear twice per padding!"
            #gauss_noise = np.random.normal(current_mean,current_std,prepadnum).astype(np.uint8)
            #gauss_noise = np.random.normal(current_mean,current_std,prepadnum)
            #print "$%(*$&%$(^#(^%(#@*^%#)(&%^)#&%^$)&%"
            #print "len of gauss noise less than 0"
            #print len(gauss_noise[gauss_noise < 0])
            #print "$%(*$&%$(^#(^%(#@*^%#)(&%^)#&%^$)&%"
            #gauss_noise[gauss_noise < 0] = np.random.choice([0,round(current_mean)])
            #gauss_noise = gauss_noise.astype(np.uint8)
            #padding = np.zeros(prepadnum,dtype=np.uint8)
            #f.write(padding)     
            #f.write(gauss_noise)
    print "pre-padding done"
    print "writing data to padded file"
    with open(origfil,'r') as X:
        thehead = header(X)
        headlen = len(thehead)
        X.seek(headlen)
        while 1:
            #dataIn = X.read(nchans*blocksize)
            dataIn = np.fromfile(X,dtype=np.uint8,count=512*blocksize)
            if len(dataIn) == 0: break
            with open(outfile,'a') as f:
                f.write(dataIn)
    print "data writing done"
    print "doing post-padding with %s samples" % postsamps
    if postsamps > 0:
        with open(outfile,'a') as f:
            current_mean = beam_mean
            current_std = beam_std
            noisecounter = 0
            diff = 1
            noiseblock = 500
            while diff != 0:
                adding = noiseblock * nchans
                noisecounter += adding
                gauss_noise = np.random.normal(current_mean,current_std,int(adding))
                #print "$%(*$&%$(^#(^%(#@*^%#)(&%^)#&%^$)&%"
                #print "len of gauss noise less than 0"
                #print len(gauss_noise[gauss_noise < 0])
                #print "$%(*$&%$(^#(^%(#@*^%#)(&%^)#&%^$)&%"
                gauss_noise[gauss_noise < 0] = np.random.choice([0,round(current_mean)])
                gauss_noise = gauss_noise.astype(np.uint8)
                f.write(gauss_noise)
                diff = postpadnum - noisecounter
                if diff < adding:
                    noiseblock = diff/nchans
                    print "diff in padsamps and counter"
                    print diff
                    print "I should only appear twice per padding!"
            #gauss_noise = np.random.normal(current_mean,current_std,postpadnum).astype(np.uint8)
            #gauss_noise = np.random.normal(current_mean,current_std,postpadnum)
            #print "$%(*$&%$(^#(^%(#@*^%#)(&%^)#&%^$)&%"
            #print "len of gauss noise less than 0"
            #print len(gauss_noise[gauss_noise < 0])
            #print "$%(*$&%$(^#(^%(#@*^%#)(&%^)#&%^$)&%"
            #gauss_noise[gauss_noise < 0] = np.random.choice([0,round(current_mean)])
            #gauss_noise = gauss_noise.astype(np.uint8)
            #padding = np.zeros(postpadnum,dtype=np.uint8)
            #f.write(padding)
            #f.write(gauss_noise)
    with open(logfile,'a') as L:
        L.write("#pre-padded %s samples\n" % int(presamps))
        L.write("#which is %s seconds\n" % (presamps*samptime))
        L.write("%s\t%s\n" % (int(presamps),presamps*samptime))
        L.write("#post-padded %s samples\n" % int(postsamps))
        L.write("#which is %s seconds\n" % (postsamps*samptime))
        L.write("%s\t%s\n" % (int(postsamps),postsamps*samptime))
        L.write("#the prepadding needs to be taken into account when looking at the other padding\n")
    print "post-padding done"
    

if __name__=="__main__":
    desc = "example usage: PAFglue.py --dir /path/to/filterbanks/ --start 2018-06-14-12:30:00 \
        --duration 30 --start-beam 0 --nbeams 16 --outdir /path/to/filterbanks/B0355+54_180614\
        --outname B0355+54"
    parser = optparse.OptionParser(description=desc)
    parser.add_option('--dir',dest='d',type='string',
        help="Directory containing filterbank files")
    parser.add_option('--start',dest='start',type='string',
        help="Start time of observations in format YYYY-MM-DD-HH:MM:SS")
    parser.add_option('--duration',dest='dur',type='int',
        help="Duration of observations in minutes")
    parser.add_option('--nbeams',dest='nbeams',type='int',
        help="Number of beams to grab")
    parser.add_option('--outname',dest='outname',type='string',
        help="Name of the output file, e.g. B0355+54")
    parser.add_option('--start-beam',dest='stbeam',type='int',
        help="Number of starting beam")
    parser.add_option('--outdir',dest='outdir',type='string',
        help="Output directory of files (creates it if it doesn't exist)")
    (opts,args) = parser.parse_args()
    #printing_n_counting(opts.d)
    Tstart = datetime.strptime(opts.start,'%Y-%m-%d-%H:%M:%S')
    Tend = Tstart + timedelta(minutes=opts.dur)
    TstartMJD = (Time(Tstart,format='datetime')).mjd
    TendMJD = (Time(Tend,format='datetime')).mjd
    startandbeam = []
    startandbeam = np.asarray(startandbeam)
    try:
        os.mkdir(opts.outdir)
    except OSError as error:
        print error
    #for B in range(opts.nbeams):
    # make an (Nbeam,3) array with [beam No., median/mean, std] to keep running
    # statistics for padding. making the noise should then just be 
    # np.random.normal(mean/median[beam No], std[beam No], numsamps).astype(np.uint8)
    for B in range(opts.stbeam,opts.stbeam+opts.nbeams):
        print "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"
        print "working with beam number %s" % B
        print "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$"
        thefiles = find_files(opts.d,B) 
        if len(thefiles) == 0:
            print "there are no files with this beam number in this directory"
            continue
        a = opts.outname + "_beam_" + str(B+1) + ".fil"
        outguy = os.path.join(opts.outdir,a)
        print "number of filterbanks for this beam: %s " % len(thefiles)
        start_file, where_is_it = find_start_file(TstartMJD,TendMJD,thefiles)
        if start_file == 0:
            print "No files within requested time range"
            continue
        thelog = opts.outname + "_beam_" + str(B+1) + "_log.txt"
        LOG = os.path.join(opts.outdir,thelog)
        with open(LOG,'w') as L:
            L.write("in-between files padding:\n")
        create_new_file(start_file,outguy,False,0)
        print "start filterbank: %s " % start_file
        print "starting counter (list index for this beam's filterbanks): %s " % where_is_it
        bigfilelen, samptime, beam_mean, beam_std = add_to_new(thefiles,where_is_it,TendMJD,outguy,LOG)
        bigboy = (bigfilelen/512) * samptime / 60
        print "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&"
        print "number of samples %s " % (bigfilelen/512.)
        print "the length of the big file is %s minutes" % bigboy
        print "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&"
        bbbeam = int(get_info(start_file)[0])
        with open(start_file) as F:
            hhhead = header(F)
        ssstart = get_headparam(hhhead,['tstart'])[0]
        #ssstart = get_info(start_file)[1]
        startandbeam = np.append(startandbeam,[bbbeam,ssstart,ssstart + bigboy*60/86400,
            beam_mean, beam_std],axis=0)
        #startandbeam.append((bbbeam,ssstart,ssstart + bigboy*60/86400))
    #startandbeam = np.asarray(startandbeam)
    startandbeam = np.reshape(startandbeam,(len(startandbeam)/5,5))
    print startandbeam
    print "earliest start"
    print min(startandbeam[:,1])
    earliest_start = min(startandbeam[:,1])
    #index_of_first = int(np.where(startandbeam[:,1] == earliest_start))
    #print max(startandbeam[:,2])
    latest_end = max(startandbeam[:,2])
    for B in range(len(startandbeam[:,0])):
        beam = int(startandbeam[B,0])
        print "beam %s" % beam
        earlydiff = startandbeam[B,1] - earliest_start
        print "diff from latest start in MJD: %s" % earlydiff
        earlydiff_samples = earlydiff *86400 / samptime
        print "diff from latest start in samples: %s" % earlydiff_samples
        latediff = latest_end - startandbeam[B,2]
        print "diff from earliest end in MJD: %s" % latediff
        latediff_samples = latediff * 86400 /samptime
        print "diff from earliest end in samples: %s" % latediff_samples
    for B in range(len(startandbeam[:,0])):
        beam = int(startandbeam[B,0])
        thelog = opts.outname + "_beam_" + str(beam+1) + "_log.txt"
        LOG = os.path.join(opts.outdir,thelog)
        with open(LOG,'a') as L:
            L.write("pre and post padding:\n")
        print "#$#$#$$#$#$#$#$$#$#$"
        print "pre- and post-padding beam number %s" % beam
        print "#$#$#$$#$#$#$#$$#$#$"
        origi = opts.outname + "_beam_" + str(beam + 1) + ".fil"
        paddy = opts.outname + "_beam_" + str(beam + 1) + "_padded.fil"
        origuy = os.path.join(opts.outdir,origi)
        padguy = os.path.join(opts.outdir,paddy)
        print "TESTING SHIT"
        print origuy
        print os.path.isfile(origuy)
        if not os.path.isfile(origuy):
            print "I'm here for some reason"
            print "it's because you're trying to pad a file that doesn't exist buddy"
            continue
        create_new_file(origuy,padguy,True,earliest_start)
        earlydiff = startandbeam[B,1] - earliest_start
        earlydiff_samples = round(earlydiff * 86400 / samptime)
        latediff = latest_end - startandbeam[B,2]
        latediff_samples = round(latediff * 86400 / samptime)
        beam_mean = startandbeam[B,3]
        beam_std = startandbeam[B,4]
        thepad(origuy,padguy,earlydiff_samples,latediff_samples,beam_mean,beam_std,LOG)
        subprocess.check_call(["rm",origuy])
    print startandbeam
         
 
    
        
