import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import optparse
import glob
import subprocess
import os

def dispdelay(DM,LOFREQ,HIFREQ):
    dconst = 4.15e+06
    delay = DM * dconst * (1.0 / LOFREQ**2 - 1.0 / HIFREQ**2) # in ms
    return delay

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
            print "hi there, I'm assuming the DADA header is 4096 bits"
        else:
            print "hi there, I'm reading a filterbank header now"
            thehead = header(F)
            headlen = len(thehead)
        #F.seek(4096+NCHAN*STARTSAMP)
        F.seek(headlen+NCHAN*STARTSAMP)
        #data = np.fromfile(F,dtype=np.uint8,count=int(NCHAN*NUMSAMP))
        data = np.fromfile(F,dtype=DTYPE,count=int(NCHAN*NUMSAMP))
        data = np.reshape(data,(-1,NCHAN)).T
    return data

def read_mask(MASK):
    with open(MASK,'r') as x:
        np.fromfile(x,dtype=np.float64,count=6)
        np.fromfile(x,dtype=np.int32,count=3)
        nzap = np.fromfile(x,dtype=np.int32,count=1)[0]
        mask_zap_chans = np.fromfile(x,dtype=np.int32,count=nzap)
    return mask_zap_chans

def plot_dedisp_ts(NCHAN,FTOP,FCHAN,DM,SAMPTIME,SAMPUSE,DATA,PLOTNAME):
    # Dedispersed time series
    ts_DD = np.zeros((1,SAMPUSE))
    for chan in np.arange(NCHAN):
        chandata = DATA[chan,:]
        chanfreq = FTOP - chan * FCHAN
        dmdelay = dispdelay(DM,chanfreq,FTOP)
        dmdelay_samp = int(np.round(dmdelay / SAMPTIME))
        #print "dmdelay {} ms, {} samples".format(dmdelay,dmdelay_samp)
        sampuse = SAMPUSE
        ts_DD = np.add(ts_DD,chandata[dmdelay_samp:dmdelay_samp + sampuse])
    print "ts DD shape {}".format(np.shape(ts_DD))
    plt.plot(ts_DD[0,:])
    plt.xlabel('Time [samples]')
    plt.ylabel('Amplitude')
    plt.title('De-dispersed time series')
    plt.savefig(PLOTNAME)
    plt.close()
    #plt.show()

def plot_bandpass(DATA,PLOTNAME):
    # summed data for each channel
    plt.plot(np.sum(DATA,axis=1))
    plt.xlabel('Channel')
    plt.ylabel('Amplitude')
    plt.title('Bandpass')
    plt.savefig(PLOTNAME)
    plt.close()

def plot_dynspec_raw(DATA,PLOTNAME):
    # raw dynamic spectra
    print "data shape {}".format(np.shape(DATA))
    plt.imshow(DATA,aspect='auto',interpolation='none',cmap='binary')
    plt.xlabel('Time [samples]')
    plt.ylabel('Frequency [channels]')
    plt.title('Dynamic spectrum')
    plt.savefig(PLOTNAME)
    plt.close()
    #plt.show()

def plot_tscrunch(AVGFACT,SAMPUSE,NCHAN,DATA,PLOTNAME): # sampuse_noDD, which datagrab
    # timescrunch 16 samples together
    nchan = NCHAN
    avgfactor = AVGFACT
    avgsampuse = int(np.floor(SAMPUSE / avgfactor))
    timmy = np.zeros((NCHAN,avgsampuse))
    for tavg in np.arange(avgsampuse):
        subdata = DATA[:, tavg * avgfactor : (tavg + 1) * avgfactor]
        subavg = np.reshape(np.sum(subdata,axis=1), (nchan,1)) / avgfactor
        timmy[:, tavg] = subavg[:,0] 
    print "timmy shape {}".format(np.shape(timmy))
    plt.imshow(timmy,aspect='auto',interpolation='none',cmap='binary')
    plt.xlabel('Time [samples]')
    plt.ylabel('Frequency [channels]')
    plt.title('Dynamic spectrum, time scrunched')
    plt.savefig(PLOTNAME)
    plt.close()
    #plt.show()

def plot_freqscrunch(OUTBANDS,NCHAN,SAMPUSE,DATA,FTOP,FCHAN,DM,SAMPTIME,PLOTNAME): # datagrab_noDD,sampuse_noDD
    # freqscrunch to 64 channels dyn spec
    outbands = OUTBANDS
    perband = int(NCHAN / outbands)
    downsamped = np.zeros((outbands,SAMPUSE))
    for band in np.arange(outbands):
        subdata = DATA[band * perband : (band + 1) * perband, :]
        bandtop = FTOP - band * perband * FCHAN
        for chan in np.arange(perband):
            chanfreq = bandtop - chan * FCHAN
            dmdelay = dispdelay(DM,chanfreq,bandtop)
            dmdelay_samp = int(np.round(dmdelay / SAMPTIME))
            #downsamped[band,:] = np.add(downsamped[band,:],subdata[chan,dmdelay_samp : dmdelay_samp + sampuse_noDD])
            downsamped[band,:] = np.add(downsamped[band,:],subdata[chan,dmdelay_samp : dmdelay_samp + SAMPUSE])
    print "downsamped shape {}".format(np.shape(downsamped))
    plt.imshow(downsamped, aspect='auto',interpolation='none',cmap='binary')
    plt.xlabel('Time [samples]')
    plt.ylabel('Frequency [channels]')
    plt.title('Dynamic spectrum, freq. scrunched')
    plt.savefig(PLOTNAME)
    plt.close()
    #plt.show()
    return downsamped

def plot_tfscrunch(SAMPUSE,AVGFACT,OUTBANDS,DATA,PLOTNAME): # needs data from freqscrunch, sampuse_noDD, avgfactor, outbands
    # freq and time scrunched
    avgfactor = AVGFACT
    outbands = OUTBANDS
    avgsampuse = int(np.floor(SAMPUSE / avgfactor))
    timmy = np.zeros((outbands,avgsampuse))
    for tavg in np.arange(avgsampuse):
        subdata = DATA[:, tavg * avgfactor : (tavg + 1) * avgfactor]
        subavg = np.reshape(np.sum(subdata,axis=1), (outbands,1)) / avgfactor
        timmy[:,tavg] = subavg[:,0]
    print "timmy shape {}".format(np.shape(timmy))
    plt.imshow(timmy,aspect='auto',interpolation='none',cmap='binary')
    plt.xlabel('Time [samples]')
    plt.ylabel('Frequency [channels]')
    plt.title('Dynamic spectrum, time and freq. scrunched')
    plt.savefig(PLOTNAME)
    plt.close()
    #plt.show()

def plot_DDF(OUTBANDS,SAMPUSE,DATA,FTOP,FCHAN,DM,NCHAN,SAMPTIME,PLOTNAME): # uses downsamped form fscrunched, sampuse, 
    perband = int(NCHAN / OUTBANDS)
    outbands = OUTBANDS
    sampuse = SAMPUSE
    samptime = SAMPTIME
    ftop = FTOP
    fchan = FCHAN
    dedispF = np.zeros((outbands,sampuse))
    for chan in np.arange(outbands):
        chandata = DATA[chan,:]
        chanfreq = ftop - chan * perband * fchan
        dmdelay = dispdelay(DM,chanfreq,ftop)
        dmdelay_samp = int(np.round(dmdelay / samptime))
        dedispF[chan,:] = np.add(dedispF[chan,:],chandata[dmdelay_samp : dmdelay_samp + sampuse]) 
    print "dedispF shape {}".format(np.shape(dedispF))
    plt.imshow(dedispF,aspect='auto',interpolation='none',cmap='binary')
    plt.xlabel('Time [samples]')
    plt.ylabel('Frequency [channels]')
    plt.title('De-dispersed dynamic spectrum, freq. scrunched')
    plt.savefig(PLOTNAME)
    plt.close()
    #plt.show()
    return dedispF

def plot_DDFT(SAMPUSE,AVGFACT,OUTBANDS,DATA,PLOTNAME): # uses dedispF data, 
    # dedisperse freq and time scrunch
    sampuse = SAMPUSE; avgfactor = AVGFACT; outbands = OUTBANDS
    avgsampuse = int(np.floor(sampuse / avgfactor))
    print "DDFT avgsampuse {}".format(avgsampuse)
    print "DDFT sampuse {}".format(sampuse)
    dedispFT = np.zeros((outbands,avgsampuse))
    for tavg in np.arange(avgsampuse):
        subdata = DATA[:, tavg * avgfactor : (tavg + 1) * avgfactor]
        subavg = np.reshape(np.sum(subdata,axis=1), (outbands,1)) / avgfactor
        dedispFT[:,tavg] = subavg[:,0]
    print "dedispFT shape {}".format(np.shape(dedispFT))
    plt.imshow(dedispFT,aspect='auto',interpolation='none',cmap='binary')
    plt.xlabel('Time [samples]')
    plt.ylabel('Frequency [channels]')
    plt.title('De-dispersed dynamic spectrum, freq. and time scrunched')
    plt.savefig(PLOTNAME)
    plt.close()
    #plt.show()

def plot_DDT(SAMPUSE,AVGFACT,NCHAN,DATA,FTOP,FCHAN,SAMPTIME,DM,PLOTNAME): # datagrab_noDD
    # dedisperse time scrunch
    sampuse = SAMPUSE; avgfactor = AVGFACT; nchan = NCHAN; ftop = FTOP
    fchan = FCHAN; samptime = SAMPTIME
    avgsampuse = int(np.floor(sampuse / avgfactor))
    print "DDT avgsampuse {}".format(avgsampuse)
    print "DDT sampuse {}".format(sampuse)
    dedispT = np.zeros((nchan,avgsampuse))
    dedispy = np.zeros((nchan,sampuse))
    for chan in np.arange(nchan):
        chandata = DATA[chan,:]
        chanfreq = ftop - chan * fchan
        dmdelay = dispdelay(DM,chanfreq,ftop)
        dmdelay_samp = int(np.round(dmdelay / samptime))
        dedispy[chan,:] = chandata[dmdelay_samp : dmdelay_samp + sampuse]     
    for tavg in np.arange(avgsampuse):
        subdata = dedispy[:, tavg * avgfactor : (tavg + 1) * avgfactor]
        subavg = np.reshape(np.sum(subdata,axis=1), (nchan,1)) / avgfactor
        dedispT[:, tavg] = subavg[:,0] 
    print "dedispy shape {}".format(np.shape(dedispy))
    print "dedispT shape {}".format(np.shape(dedispT))
    plt.imshow(dedispT,aspect='auto',interpolation='none',cmap='binary')
    plt.xlabel('Time [samples]')
    plt.ylabel('Frequency [channels]')
    plt.title('De-dispersed dynamic spectrum, time scrunched')
    plt.savefig(PLOTNAME)
    plt.close()
    #plt.show()

if __name__=="__main__":
    desc = """Candidate plotter. Give me a file and info and I'll make nice plots. 
        I can see whether it's .dada or .fil so you don't have to worry about that. 
        What I need from you is \n 
        --ftop <top of band in MHz> \n 
        --fchan <channel bandwidth in MHz> \n
        --nchan <number of channels> \n
        --samptime <sample time in ms> \n
        --data <path to .fil/.dada> \n
        --out <path to where plots go> (dir is created if it doesn't exist) \n
        --dm <DM of candidate> \n
        --samp <sample of candidate in the file> \n
        If you don't like the default values you can also specify: \n
        --fscrunch <freq scrunch to X channels> \n
        --tscrunch <time scrunch X samples together> \n
        --dtype <data type> \n
        -w <width of candidate in samples> \n
        --beamno <beam numer> (only used for plot name) \n
        What you get are eight plots: i) dedispersed time series \n
        ii) dynamic spectrum \n iii) time scrunched (downsampled) dynamic spectrum \n
        iv) frequency scrunched dynamic spectrum \n v) time and freq scrunched dynamic spectrum \n
        vi) de-dispersed dynamic spectrum, freq scrunched \n 
        vii) de-dispersed dynamic spectrum, time scrunched (downsampled) \n
        viii) de-dispersed dynamic spectrum, freq and time scrunched \n
        You can clone the code from https://github.com/ghenning/PAFcode.git and play with it as you like.
        NEW: added manual zapping and zapping from rfifind mask: \n
        --zap <startchan:endchan-startchan:endchan> i.e. zap as many ranges as you want but separate ranges
        with a dash (-) and start/stop of range with a colon (:), NO SPACES! \n
        --mask <path to rfifind mask> reads the zapped channels from the mask and applies it to plots"""
    parser = optparse.OptionParser(description=desc)
    parser.add_option('--ftop',dest='ftop',type='float',
        help="Top part of band in MHz (Default=1492.203704 [PAF])", default=1492.203704)
    parser.add_option('--fchan',dest='fchan',type='float',
        help="channel bandwidth in MHz (Default=0.592593 [PAF])", default=0.592593)
    parser.add_option('--nchan',dest='nchan',type='int',
        help="number of channels (Default=512 [PAF])", default=512)
    parser.add_option('--samptime',dest='samptime',type='float',
        help="sample time in ms (Default=54e-03 [PAF])", default=54e-03)
    parser.add_option('--data',dest='data',type='string', 
        help="data file")
    parser.add_option('--fscrunch',dest='fscrunch',type='int',
        help="frequency scrunch to X channels (Default=64)", default=64)
    parser.add_option('--tscrunch',dest='tscrunch',type='int',
        help="time scrunch X samples together (Default=16)", default=16)
    parser.add_option('--dtype',dest='dtype',type='str',
        help="data type (Default=np.uint8)", default=np.uint8)
    parser.add_option('--out',dest='out',type='string',
        help="where to put your plots")
    parser.add_option('--dm',dest='DM',type='float',
        help="DM [Default=0]",default=0)
    parser.add_option('--samp',dest='samp',type='int',
        help="Sample number of candidate")
    parser.add_option('-w',dest='w',type='int',
        help="Width of candidate in samples (Default=10)",default=10)
    parser.add_option('--beamno',dest='beamno',type='int',
        help="beam number (only used for name of plot) (Default=1)",default=1)
    parser.add_option('--mask',dest='mask',type='str',
        help="rfifind mask to apply",default='')
    parser.add_option('--zap',dest='zap',type='str',
        help="zap specific channels (<x:y-x:y-x:y>) (no spaces!)",default='')
    (opts,args) = parser.parse_args()
    if opts.mask:
        print "Using mask {}".format(os.path.basename(opts.mask))
        MASK_ME = True
    if not opts.mask:
        print "no mask applied"
        MASK_ME = False
    if opts.zap:
        print "manual zapping"
        #print opts.zap
        zapsplit = opts.zap.split('-')
        print "zap these bastards {}".format(zapsplit)
        ZAP_ME = True
    if not opts.zap:
        print "no manual zapping"
        ZAP_ME = False
    fbot = opts.ftop - opts.nchan * opts.fchan
    ddelay = dispdelay(opts.DM,fbot,opts.ftop)
    ddelay_samp = int(np.round(ddelay/opts.samptime))
    padding = int(np.round(20 / opts.samptime))
    startsamp = int(opts.samp - padding)
    sampuse = 2*padding + opts.w
    to_grab = 2*padding + opts.w + ddelay_samp
    sampuse_noDD = 2*padding + opts.w + ddelay_samp
    to_grab_noDD = 2*padding + opts.w + 2*ddelay_samp
    timeint = round(opts.samp*opts.samptime*10**-3,2) # only used for name of plot
    base = os.path.splitext(os.path.basename(opts.data))[0]
    ext = os.path.splitext(os.path.basename(opts.data))[-1]
    print "ext {}".format(ext)
    isfil = True
    if ext=='.dada':
        print "I see you're using .dada"
        isfil = False
    if ext=='.fil':
        print "I see you're using .fil"
        isfil = True
    if not os.path.exists(opts.out):
        try:
            subprocess.check_call(["mkdir",opts.out])
        except OSError as error:
            print error
    dat = grab_data(opts.data, startsamp, to_grab, opts.nchan,opts.dtype,isfil)
    if MASK_ME:
        maskzaps = read_mask(opts.mask)
        dat[maskzaps,:] = 0
    if ZAP_ME:
        for i in range(len(zapsplit)):
            print "zapping {}".format(zapsplit[i].strip())
            zapchunk = zapsplit[i].strip()
            zapstart = int(zapchunk.split(':')[0])
            zapend = int(zapchunk.split(':')[1])
            zaprange = np.arange(zapstart,zapend)
            dat[zaprange,:] = 0
    print "data shape {}".format(np.shape(dat))
    plotn = "{}_time{}_DM{}_beam{}_bp.png".format(base,timeint,int(opts.DM),opts.beamno)
    plotname = os.path.join(opts.out,plotn)
    print plotname 
    plot_bandpass(dat,plotname)
    plotn = "{}_time{}_DM{}_beam{}_tsDD.png".format(base,timeint,int(opts.DM),opts.beamno)
    plotname = os.path.join(opts.out,plotn)
    print plotname 
    plot_dedisp_ts(opts.nchan,opts.ftop,opts.fchan,opts.DM,opts.samptime,sampuse,dat,plotname)
    plotn = "{}_time{}_DM{}_beam{}_dynspec.png".format(base,timeint,int(opts.DM),opts.beamno)
    plotname = os.path.join(opts.out,plotn)
    print plotname 
    plot_dynspec_raw(dat,plotname)
    plotn = "{}_time{}_DM{}_beam{}_dynspecT.png".format(base,timeint,int(opts.DM),opts.beamno)
    plotname = os.path.join(opts.out,plotn)
    print plotname 
    plot_tscrunch(opts.tscrunch,sampuse_noDD,opts.nchan,dat,plotname)
    dat = grab_data(opts.data,startsamp,to_grab_noDD,opts.nchan,opts.dtype,isfil)
    if MASK_ME:
        dat[maskzaps,:] = 0
    if ZAP_ME:
        for i in range(len(zapsplit)):
            print "zapping {}".format(zapsplit[i].strip())
            zapchunk = zapsplit[i].strip()
            zapstart = int(zapchunk.split(':')[0])
            zapend = int(zapchunk.split(':')[1])
            zaprange = np.arange(zapstart,zapend)
            dat[zaprange,:] = 0
    print "data shape {}".format(np.shape(dat))
    plotn = "{}_time{}_DM{}_beam{}_dynspecF.png".format(base,timeint,int(opts.DM),opts.beamno)
    plotname = os.path.join(opts.out,plotn)
    print plotname 
    binned_data = plot_freqscrunch(opts.fscrunch,opts.nchan,sampuse_noDD,dat,opts.ftop,opts.fchan,opts.DM,opts.samptime,plotname)
    plotn = "{}_time{}_DM{}_beam{}_dynspecFT.png".format(base,timeint,int(opts.DM),opts.beamno)
    plotname = os.path.join(opts.out,plotn)
    print plotname 
    plot_tfscrunch(sampuse_noDD,opts.tscrunch,opts.fscrunch,binned_data,plotname)
    plotn = "{}_time{}_DM{}_beam{}_dynspecDDF.png".format(base,timeint,int(opts.DM),opts.beamno)
    plotname = os.path.join(opts.out,plotn)
    print plotname 
    binned_data = plot_DDF(opts.fscrunch,sampuse,binned_data,opts.ftop,opts.fchan,opts.DM,opts.nchan,opts.samptime,plotname)
    plotn = "{}_time{}_DM{}_beam{}_dynspecDDFT.png".format(base,timeint,int(opts.DM),opts.beamno)
    plotname = os.path.join(opts.out,plotn)
    print plotname 
    plot_DDFT(sampuse,opts.tscrunch,opts.fscrunch,binned_data,plotname)
    plotn = "{}_time{}_DM{}_beam{}_dynspecDDT.png".format(base,timeint,int(opts.DM),opts.beamno)
    plotname = os.path.join(opts.out,plotn)
    print plotname 
    plot_DDT(sampuse,opts.tscrunch,opts.nchan,dat,opts.ftop,opts.fchan,opts.samptime,opts.DM,plotname)
    
