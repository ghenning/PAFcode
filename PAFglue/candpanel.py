import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import argparse
import os
import subprocess

def dispdelay(DM,LOFREQ,HIFREQ):
    dconst = 4.15e+06
    delay = DM * dconst * (1. / LOFREQ**2 - 1. / HIFREQ**2) # in ms
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
            #print "hi there, I'm assuming the DADA header is 4096 bits"
        else:
            #print "hi there, I'm reading a filterbank header now"
            thehead = header(F)
            headlen = len(thehead)
        F.seek(headlen+NCHAN*STARTSAMP)
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

def get_ts(NCHAN,FTOP,FCHAN,DM,SAMPTIME,SAMPUSE,DATA):
    ts_DD = np.zeros((1,SAMPUSE))
    for chan in np.arange(NCHAN):
        chandata = DATA[chan,:]
        chanfreq = FTOP - chan * FCHAN
        dmdelay = dispdelay(DM,chanfreq,FTOP)
        dmdelay_samp = int(np.round(dmdelay / SAMPTIME))
        sampuse = SAMPUSE
        ts_DD = np.add(ts_DD,chandata[dmdelay_samp:dmdelay_samp + sampuse])
    return ts_DD

def DM_T_plot(TS,DM,TOGRAB,TSAMP,COL):
    ax[1,0].imshow(TS,aspect='auto',interpolation='hamming',cmap=COL)
    ax[1,0].set_ylabel('DM')
    ax[1,0].set_xlabel('Time (ms)')
    ms_range = int(TOGRAB * TSAMP)
    ran = ax[1,0].get_xlim()
    ax[1,0].set_xticks(np.linspace(ran[0],ran[-1]),11)
    center = np.sum(ran)/2
    ticran = np.linspace(ran[0],ran[-1],11)
    ax[1,0].set_xticks([])
    ax[1,0].set_xticks(ticran)
    tic = ax[1,0].get_xticks()
    ticlab = np.round(TSAMP * tic)
    ticlab -= np.max(ticlab)/2
    ax[1,0].set_xticklabels(np.round(ticlab))
    ax[1,0].tick_params(axis='x',which='minor',bottom=False)
    ytic = ax[1,0].get_yticks()
    newy = ytic[1:-1]
    yt = np.linspace(DM[0],DM[-1],len(newy))
    ax[1,0].set_yticks(newy)
    ax[1,0].set_yticklabels(yt)
    ax[1,0].set_title('DM(t)')

def TS_plot(TS):
    ax[0,0].tick_params(axis='y',which='both',left=False,labelleft=False)
    ax[0,0].plot(TS,color='k')
    ax[0,0].set_ylabel('Intensity')
    ax[0,0].set_xlabel('Sample')
    ax[0,0].set_title('Dedispersed time series')
    ax[0,0].set_xlim(0,np.shape(TS)[0])

def DS_plot(OUTBANDS,AVGFAC,DM,NCHAN,DATA,FTOP,FCHAN,TSAMP,SAMPUSE,COL):
    perband = int(NCHAN/OUTBANDS)
    ds = np.zeros((OUTBANDS,SAMPUSE))
    ds2 = np.zeros((OUTBANDS,2*SAMPUSE))
    for band in np.arange(OUTBANDS):
        subdata = DATA[band * perband : (band + 1) * perband, :]
        bandtop = FTOP - band * perband * FCHAN
        for chan in np.arange(perband):
            chanfreq = bandtop - chan * FCHAN
            dmdelay = dispdelay(DM,chanfreq,bandtop)
            dmdelay_samp = int(np.round(dmdelay/TSAMP))
            ds[band,:] = np.add(ds[band,:],subdata[chan,dmdelay_samp : dmdelay_samp + SAMPUSE])
            ds2[band,:] = np.add(ds2[band,:],subdata[chan,dmdelay_samp : dmdelay_samp + SAMPUSE*2])
    avgsampuse = int(np.floor(SAMPUSE / AVGFAC))
    avgsampuse2 = int(np.floor(SAMPUSE*2 / AVGFAC))
    timmy = np.zeros((OUTBANDS,avgsampuse))
    timmy2 = np.zeros((OUTBANDS,avgsampuse2))
    for tavg in np.arange(avgsampuse):
        subdata = ds[:, tavg * AVGFAC : (tavg + 1) * AVGFAC]
        subavg = np.reshape(np.sum(subdata,axis=1),(OUTBANDS,1))/AVGFAC
        timmy[:,tavg] = subavg[:,0]
    for tavg in np.arange(avgsampuse2):
        subdata = ds[:, tavg * AVGFAC : (tavg + 1) * AVGFAC]
        subavg = np.reshape(np.sum(subdata,axis=1),(OUTBANDS,1))/AVGFAC
        timmy2[:,tavg] = subavg[:,0]
    ddd = np.zeros((OUTBANDS,avgsampuse))
    for chan in np.arange(OUTBANDS):
        chandata = timmy2[chan,:]
        chanfreq = FTOP - chan * perband * FCHAN
        dmdelay = dispdelay(DM,chanfreq,FTOP)
        dmdelay_samp = int(np.round(dmdelay/TSAMP/AVGFAC))
        ddd[chan,:] = np.add(ddd[chan,:],chandata[dmdelay_samp:dmdelay_samp+avgsampuse])
    ax[0,1].imshow(timmy,aspect='auto',interpolation='hamming',cmap=COL)
    ax[1,1].imshow(ddd,aspect='auto',interpolation='hamming',cmap=COL)
    ends = np.where(ddd[-1]==0.)[0][0]
    ax[1,1].set_xlim(0,ends-2)
    ax[0,1].set_ylabel("Subband ({} channels per)".format(perband))
    ax[1,1].set_ylabel("Subband ({} channels per)".format(perband))
    ax[0,1].set_xlabel("Sample (downsampled by {})".format(AVGFAC))
    ax[1,1].set_xlabel("Sample (downsampled by {})".format(AVGFAC))
    ax[0,1].set_title('Dynamic spectrum')
    ax[1,1].set_title('Dedispersed dynamic spectrum')

def testing():
    f = "/home/henning/work/fakedata/7bfake.fil"
    nchan = 512
    ftop = 1510
    fchan = .585938
    fbot = ftop - nchan * fchan
    samptime = 54.61333e-3
    thesamp = 91553
    halfsec = int(.5e3/samptime)
    tograb = halfsec
    startsamp = thesamp - int(.5*tograb)
    dt = np.uint8
    dms = np.arange(450,551)
    max_dispdelay = dispdelay(dms[-1],fbot,ftop)
    max_dispdelay_samp = int(max_dispdelay/samptime)
    dat = grab_data(f,startsamp,tograb+max_dispdelay_samp,nchan,dt,True)
    fluffer = int(50/samptime)
    dat2 = grab_data(f,thesamp-fluffer,2*max_dispdelay_samp + 10*fluffer,nchan,dt,True)
    tees = np.zeros((len(dms),tograb))
    
    for i in range(len(dms)):
        tmp_dm = dms[i]
        tmp_ts = get_ts(nchan,ftop,fchan,tmp_dm,samptime,tograb,dat)
        tees[i,:] = tmp_ts
    
    fig,ax = plt.subplots(2,2,figsize=(12,8))
    
    DM_T_plot(tees,dms)
    TS_plot(tees[len(tees)/2])
    DS_plot(128,8,500,nchan,dat2,ftop,fchan,samptime,max_dispdelay_samp+fluffer)
    plt.tight_layout()
    plt.show()
    print np.shape(tees)

if __name__=="__main__":
    desc = "Candidate plotter. Plots four panels for a candidate from \
            fil or dada: Dedispersed time series, DM(t), and dynamic spectra\
            (dedispersed/non-dedispersed). Too see standard options, use the\
            --print_defaults flag."
    parser = argparse.ArgumentParser(description=desc)
    optional = parser._action_groups.pop()
    required = parser.add_argument_group('required args')
    required.add_argument('--data',type=str,
            help="Filterbank file")
    required.add_argument('--out',type=str,
            help="Plot output directory")
    required.add_argument('--dm',type=float,
            help="Candidate DM")#,required=True)
    required.add_argument('--samp',type=int,
            help="Sample number of candidate")#,required=True)
    optional.add_argument('--cx',action='store_true',
            help="C+/CX receiver default settings",default=False)
    optional.add_argument('--PAF',action='store_true',
            help="PAF default settings",default=False)
    optional.add_argument('--seven',action='store_true',
            help="7beam default settings",default=False)
    optional.add_argument('--ftop',type=float,
            help="Top of the band in MHz, default is PAF")
    optional.add_argument('--fchan',type=float,
            help="Channel bandwidth in MHz, default is PAF")
    optional.add_argument('--nchan',type=int,
            help="Number of channels, default is PAF")
    optional.add_argument('--samptime',type=float,
            help="Sample time in ms, default is PAF")
    optional.add_argument('--downsamp',type=int,
            help="Downsampling factor (power of 2), default is 8",
            default=8)
    optional.add_argument('--subbands',type=int,
            help="Number of subbands (power of 2), default is 128",
            default=128)
    optional.add_argument('--dtype',type=str,
            help="Data type, default is np.uint8",default=np.uint8)
    optional.add_argument('--beamno',type=int,
            help="Beam number (used for outfile name), default is 0",
            default=0)
    optional.add_argument('--print_defaults',action='store_true',
            help="Print default settings for receivers",default=False)
    optional.add_argument('--window',type=float,
            help="Width of time series plot in s, default is 0.5 seconds",
            default=0.5)
    optional.add_argument('--dm_range',type=int,
            help="DM range around candidate DM for DM-t plot, default is 100",
            default=100)
    optional.add_argument('--cmap',type=str,
            help="Cmap for plots, default is jet. For black and white, choose 'binary'",
            default='jet')
    optional.add_argument('--mask',type=str,
            help="Path to rfifind mask to apply")
    optional.add_argument('--zap',type=str,
            help="Manually zap channels (--zap x:y-x:y-x:y) (no spaces)")
    #optional.add_argument()
    parser._action_groups.append(optional)
    args = parser.parse_args()

    # standards for seven beam, paf, cx
    sb = [1510,.585938,512,54.61333e-3]
    pp = [1451.962963,0.44907407421875,512,216e-3]
    cc = [8000,0.976562,4096,131.072e-3]

    # Printing default settings
    if args.print_defaults:
        print "7-beam{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}"\
            .format("\n","FTOP (MHz)","\t",sb[0],"\n",
                    "FCHAN (MHz)","\t",sb[1],"\n",
                    "NCHAN","\t\t",sb[2],"\n",
                    "TSAMP (ms)","\t",sb[3])
        print "\n"
        print "PAF{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}"\
            .format("\n","FTOP (MHz)","\t",pp[0],"\n",
                    "FCHAN (MHz)","\t",pp[1],"\n",
                    "NCHAN","\t\t",pp[2],"\n",
                    "TSAMP (ms)","\t",pp[3])
        print "\n"
        print "C+/CX{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}{}"\
            .format("\n","FTOP (MHz)","\t",cc[0],"\n",
                    "FCHAN (MHz)","\t",cc[1],"\n",
                    "NCHAN","\t\t",cc[2],"\n",
                    "TSAMP (ms)","\t",cc[3])
        exit()

    # Checking if inputs are missing
    # I could use 'required=True' when adding these arguments,
    # but then they would need to be present when using the
    # --print_defaults flag
    if args.data==None:
        print "\nLooks like you forgot to feed me data..."
        print "use --data to get rid of me\n"
        parser.print_help()
        exit()
    if args.out==None:
        print "\nLooks like you forgot to set an output directory..."
        print "use --out to get rid of me\n"
        parser.print_help()
        exit()
    if args.dm==None:
        print "\nLooks like you hate dedispersion..."
        print "use --dm to get rid of me\n"
        parser.print_help()
        exit()
    if args.samp==None:
        print "\nI don't know where in the data your candidate is..."
        print "use --samp to get rid of me\n"
        parser.print_help()
        exit()

    # Setting basic inputs
    if args.cx: # CX presets
        print "Using CX presets"
        ftop = cc[0]
        fchan = cc[1]
        nchan = cc[2]
        tsamp = cc[3]
    if args.seven: # 7beam presets
        print "Using 7beam presets"
        ftop = sb[0]
        fchan = sb[1]
        nchan = sb[2]
        tsamp = sb[3]
    if (args.PAF or args.ftop==None) and (not args.cx and not args.seven): # PAF presets
        print "Using PAF presets"
        ftop = pp[0]
        fchan = pp[1]
        nchan = pp[2]
        tsamp = pp[3]
    if not args.PAF and not args.cx and not args.seven: # manual inputs
        ftop = args.ftop
        fchan = args.fchan
        nchan = args.nchan
        tsamp = args.samptime
    fbot = ftop - nchan * fchan # bottom of the band
    f = args.data
    thesamp = args.samp
    win = args.window * 1e3 # plot window in ms
    dt = args.dtype
    dm = args.dm
    dm_ran = args.dm_range
    dms = np.arange(int(dm)-dm_ran/2,int(dm)+dm_ran/2+1) # dm range 
    subbands = args.subbands
    downsamp = args.downsamp
    beamno = args.beamno
    color = args.cmap
    winsamp = int(win/tsamp) # plot window in samples
    tograb = winsamp
    startsamp = thesamp - int(.5*tograb)
    max_dispdelay = dispdelay(dms[-1],fbot,ftop)
    max_dispdelay_samp = int(max_dispdelay/tsamp)

    # Status
    print "Working on {}".format(f)

    # Check if output dir exists, otherwise try to make it
    if not os.path.exists(args.out):
        try:
            subprocess.check_call(["mkdir",args.out])
        except OSError as error:
            print error

    # Check extension of data (fil vs dada)
    base = os.path.splitext(os.path.basename(f))[0]
    ext = os.path.splitext(os.path.basename(f))[-1]
    if ext=='.dada':
        is_fil = False
    if ext=='.fil':
        is_fil = True

    # Checking if mask or zapping is requested
    MASK_ME = False
    if not args.mask==None:
        print "Using mask {}".format(os.path.basename(args.mask))
        MASK_ME = True
    ZAP_ME = False
    if not args.zap==None:
        zapsplit = args.zap.split('-')
        print "Zapping channels {}".format(zapsplit)
        ZAP_ME = True

    # Grabbing data
    dat = grab_data(f,startsamp,tograb+max_dispdelay_samp,nchan,dt,is_fil)
    fluffer = int(50/tsamp)
    dat2 = grab_data(f,thesamp-fluffer,2*max_dispdelay_samp + 10*fluffer,nchan,dt,is_fil)

    # Initiate time series
    tees = np.zeros((len(dms),tograb))

    # Apply mask/zapping
    if MASK_ME:
        maskzaps = read_mask(args.mask)
        dat[maskzaps,:] = 0
        dat2[maskzaps,:] = 0
    if ZAP_ME:
        for i in range(len(zapsplit)):
            zapchunk = zapsplit[i].strip()
            print "zapping {}".format(zapchunk)
            zapstart = int(zapchunk.split(':')[0])
            zapend = int(zapchunk.split(':')[1])
            zaprange = np.arange(zapstart,zapend+1)
            print "zaprange {}".format(zaprange)
            dat[zaprange,:] = 0
            dat2[zaprange,:] = 0
    
    # getting the time series
    for i in range(len(dms)):
        tmp_dm = dms[i]
        tmp_ts = get_ts(nchan,ftop,fchan,tmp_dm,tsamp,tograb,dat)
        tees[i,:] = tmp_ts
    
    # creating panel for plots
    fig,ax = plt.subplots(2,2,figsize=(12,8))
    
    # DM vs T plot
    DM_T_plot(tees,dms,tograb,tsamp,color)

    # Dedispersed time series
    TS_plot(tees[len(tees)/2])

    # Dynamic spectra plots
    DS_plot(subbands,downsamp,dm,nchan,dat2,ftop,fchan,tsamp,max_dispdelay_samp+fluffer,color)

    # Creating plot name
    time = thesamp * tsamp * 1e-3
    time = np.round(time,2)
    dm = np.round(dm,2)
    plotn = "{}_time{}_DM{}_beam{}.png".format(base,time,dm,beamno) 
    plotname = os.path.join(args.out,plotn)

    # Save plot
    plt.tight_layout()
    plt.savefig(plotname,dpi=90)
    print "Done"
    #plt.show()


