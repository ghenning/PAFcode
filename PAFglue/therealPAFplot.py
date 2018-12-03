import numpy as np
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

def dispdelay(DM,LOFREQ,HIFREQ):
    dconst = 4.15e+06
    delay = DM * dconst * (1.0 / LOFREQ**2 - 1.0 / HIFREQ**2) # in ms
    return delay

def grab_data(FILE,STARTSAMP,NUMSAMP,NCHAN,DTYPE):
    with open(FILE,'r') as F:
        thehead = header(F)
        headlen = len(thehead)
        F.seek(headlen + NCHAN * STARTSAMP)
        data = np.fromfile(F,dtype=DTYPE,count=NCHAN*NUMSAMP) # make dtype option as well?
        data = np.reshape(data, (-1,NCHAN)).T 
    return data

def siftcands(CANDFILE):
    cand = np.genfromtxt(CANDFILE)
    filtcand = cand[(cand[:,0] >= 10.0) & (cand[:,5] >= 25.0) & (cand[:,10] <= 65536)]
    # some other sifting stuff in the future hopefully
    return filtcand

def plot_dedisp_ts(NCHAN,FTOP,FCHAN,DM,SAMPTIME,SAMPUSE,DATA,PLOTNAME):
    # Dedispersed time series
    ts_DD = np.zeros((1,SAMPUSE))
    for chan in np.arange(NCHAN):
        chandata = DATA[chan,:]
        chanfreq = FTOP - chan * FCHAN
        dmdelay = dispdelay(DM,chanfreq,FTOP)
        dmdelay_samp = int(np.round(dmdelay / SAMPTIME))
        #print "dmdelay {} ms, {} samples".format(dmdelay,dmdelay_samp)
        ts_DD = np.add(ts_DD,chandata[dmdelay_samp:dmdelay_samp + sampuse])
    print "ts DD shape {}".format(np.shape(ts_DD))
    plt.plot(ts_DD[0,:])
    plt.xlabel('Time [samples]')
    plt.ylabel('Amplitude')
    plt.title('De-dispersed time series')
    plt.savefig(PLOTNAME)
    plt.close()
    #plt.show()

def plot_dynspec_raw(DATA,PLOTNAME):
    # raw dynamic spectra
    plt.imshow(DATA,aspect='auto',interpolation='none',cmap='binary')
    plt.xlabel('Time [samples]')
    plt.ylabel('Frequency [channels]')
    plt.title('Dynamic spectrum')
    plt.savefig(PLOTNAME)
    plt.close()
    #plt.show()

def plot_tscrunch(AVGFACT,SAMPUSE,NCHAN,DATA,PLOTNAME): # sampuse_noDD, which datagrab
    # timescrunch 16 samples together
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

def plot_freqscrunch(OUTBANDS,NCHAN,SAMPUSE,DATA,FTOP,FCHAN,DM,PLOTNAME): # datagrab_noDD,sampuse_noDD
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
            dmdelay_samp = int(np.round(dmdelay / samptime))
            downsamped[band,:] = np.add(downsamped[band,:],subdata[chan,dmdelay_samp : dmdelay_samp + sampuse_noDD])
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

def plot_DDF(OUTBANDS,SAMPUSE,DATA,FTOP,FCHAN,DM,NCHAN,PLOTNAME): # uses downsamped form fscrunched, sampuse, 
    perband = int(NCHAN / OUTBANDS)
    outbands = OUTBANDS
    sampuse = SAMPUSE
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

def plot_DDT(SAMPUSE,AVGFACT,NCHAN,DATA,FTOP,FCHAN,SAMPTIME,PLOTNAME): # datagrab_noDD
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
    desc = "Hi I'm a description"
    parser = optparse.OptionParser(description=desc)
    parser.add_option('--ftop',dest='ftop',type='float'
        ,help="Top part of band (Default is for PAF)", default=1492.203704)
    parser.add_option('--fchan',dest='fchan',type='float',
        help="channel bandwidth (Default is for PAF)", default=0.592593)
    parser.add_option('--nchan',dest='nchan',type='int',
        help="number of channels (Default is for PAF)", default=512)
    parser.add_option('--samptime',dest='samptime',type='float',
        help="sample time in ms (Default is for PAF)", default=54e-03)
    parser.add_option('--candfile',dest='candfile',type='string',
        help="candidate file from coincidencer")
    parser.add_option('--datadir',dest='datadir',type='string',
        help="directory of filterbanks")
    parser.add_option('--fscrunch',dest='fscrunch',type='int',
        help="frequency scrunch to X channels (Default: 64)",default=64)
    parser.add_option('--tscrunch',dest='tscrunch',type='int',
        help="time scrunch X samples together (Default: 16)",default=16)
    parser.add_option('--dtype',dest='dtype',type='str',
        help="data type (Default: np.uint8)",default=np.uint8)
    # maybe just hardcode this if it is being troublesome, seems to work fine
    (opts,args) = parser.parse_args()
    dir1 = os.path.dirname(opts.candfile)
    plotpath = os.path.join(dir1,"plots")
    if not os.path.exists(plotpath):
        try:
            subprocess.check_call(["mkdir",plotpath])
        except OSError as error:
            print error    
    ftop = opts.ftop
    fchan = opts.fchan
    nchan = opts.nchan
    fbot = ftop - nchan * fchan
    samptime = opts.samptime
    #print "ftop {}".format(ftop)
    #print "fchan {}".format(fchan)
    #print "nchan {}".format(nchan)
    #print "fbot {}".format(fbot)
    #print "samptime {}".format(samptime)
    #print "bla {}".format(opts.dtype)
    # make dir for plots
    # read cand file
    # sift candidates
    # loop through candidates
    # find corresponding beam
    # candidate info
    # make plots
    candidates = siftcands(opts.candfile)
    print "siftcand shape {}".format(np.shape(candidates))
    fils = glob.glob(os.path.join(opts.datadir,"*.fil"))
    print "fils {}".format(fils)
    for cand in candidates:
        if cand[0] != cand[12]:
            print "this event is stronger in another beam, skip me please"
            continue
        DM = cand[5]
        timeint = int(cand[2])
        ddelay = dispdelay(DM,fbot,ftop)
        ddelay_samp = int(np.round(ddelay / samptime))
        candwidth = 2 ** cand[3]
        padding = int(np.round(20 / samptime)) 
        startsamp = int(cand[1] - padding)
        W_samp = int(np.round(candwidth / samptime))
        sampuse = 2*padding + W_samp
        sampuse_noDD = 2*padding + W_samp + ddelay_samp
        to_grab = 2*padding + W_samp + ddelay_samp
        to_grab_noDD = 2*padding + W_samp + 2*ddelay_samp
        beamNO = int(cand[-1])
        print "beam number {}".format(beamNO)
        for fil in fils:
            findme = 'beam_{}'.format(beamNO)
            name = os.path.basename(fil)
            if name.find(findme) > 0:
                thefil = fil
                break
        datagrab = grab_data(thefil,startsamp,to_grab,nchan,opts.dtype)
        #datagrab = grab_data(fil,startsamp,to_grab_noDD,nchan)
        plotn = "beam{}_time{}_DM{}_tsDD.png".format(beamNO,timeint,int(DM))
        plotname = os.path.join(plotpath,plotn)
        print plotname
        plot_dedisp_ts(nchan,ftop,fchan,DM,samptime,sampuse,datagrab,plotname)
        plotn = "beam{}_time{}_DM{}_dynspec.png".format(beamNO,timeint,int(DM))
        plotname = os.path.join(plotpath,plotn)
        print plotname
        plot_dynspec_raw(datagrab,plotname)
        plotn = "beam{}_time{}_DM{}_dynspecT.png".format(beamNO,timeint,int(DM))
        plotname = os.path.join(plotpath,plotn)
        print plotname
        plot_tscrunch(opts.tscrunch,sampuse_noDD,nchan,datagrab,plotname)
        datagrab = grab_data(thefil,startsamp,to_grab_noDD,nchan,opts.dtype)
        plotn = "beam{}_time{}_DM{}_dynspecF.png".format(beamNO,timeint,int(DM))
        plotname = os.path.join(plotpath,plotn)
        print plotname
        binned_data = plot_freqscrunch(opts.fscrunch,nchan,sampuse_noDD,datagrab,ftop,fchan,DM,plotname)
        plotn = "beam{}_time{}_DM{}_dynspecFT.png".format(beamNO,timeint,int(DM))
        plotname = os.path.join(plotpath,plotn)
        print plotname
        plot_tfscrunch(sampuse_noDD,opts.tscrunch,opts.fscrunch,binned_data,plotname)
        plotn = "beam{}_time{}_DM{}_dynspecDDF.png".format(beamNO,timeint,int(DM))
        plotname = os.path.join(plotpath,plotn)
        print plotname
        binned_data = plot_DDF(opts.fscrunch,sampuse,binned_data,ftop,fchan,DM,nchan,plotname)
        plotn = "beam{}_time{}_DM{}_dynspecDDFT.png".format(beamNO,timeint,int(DM))
        plotname = os.path.join(plotpath,plotn)
        print plotname
        plot_DDFT(sampuse,opts.tscrunch,opts.fscrunch,binned_data,plotname)
        plotn = "beam{}_time{}_DM{}_dynspecDDT.png".format(beamNO,timeint,int(DM))
        plotname = os.path.join(plotpath,plotn)
        print plotname
        plot_DDT(sampuse,opts.tscrunch,nchan,datagrab,ftop,fchan,samptime,plotname)









