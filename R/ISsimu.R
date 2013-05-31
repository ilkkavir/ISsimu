# 
# Incoherent scatter radar simulator 
# 
# 
# I. Virtanen 2010, 2012
#

writeTimestampsFile.gdf <- function(prefix='data-',extension='.gdf',nstart=1,nend=1,times=c(0),fname='timestamps.log'){


  fid <- file(fname,'w')
  cat('',file=fid)
  times <- rep(times,length.out=(nend-nstart+1))

  for(k in seq(nstart,nend)){

    cnum <- sprintf('%.0f',k)
    nc   <- nchar(cnum)
    stmpstr <- paste(prefix,substr('000000',1,6-nc),cnum,'.gdf ',sprintf('%20.8f',times[k-nstart+1]),sep='')
    cat(stmpstr,'\n',file=fid,append=TRUE)
    print(stmpstr)
  }

  close(fid)

}

TXenv <- function(exp){
# 
# TX envelope of the given experiment
#
# exp is a list with entries "IPP", "code", and "baudLength"
#  IPP is a vector of inter-pulse periods
#  code is a list of complex vectors giving the phase codes of each pulse
#  baudLength is the length of a single baud in the phase-coding; each element of each code is repeated baudLength times in the envelope
# 
# 
# I. Virtanen 2010, 2012
#

  # mach the lengths of IPP and code cycles
  nipp  <- length(exp$IPP)
  ncode <- length(exp$code)
  nmax  <- max(nipp,ncode)
  nmin  <- min(nipp,ncode)
  k <- 1
  while((k*nmax)%%nmin) k <- k+1
  exp$IPP <- rep(exp$IPP,(k*nmax/nipp))
  exp$code <- rep(exp$code,(k*nmax/ncode))

  env <- rep((0+0i),sum(exp$IPP))

  env[1:(length(exp$code[[1]])*exp$baudLength)] <- rep(exp$code[[1]],each=exp$baudLength)
  if(length(exp$IPP)>1){
    for(k in seq(2,length(exp$IPP))){
      env[(sum(exp$IPP[1:(k-1)])+1):(sum(exp$IPP[1:(k-1)])+(length(exp$code[[k]])*exp$baudLength))] <- rep(exp$code[[k]],each=exp$baudLength)
    }
  }

  return(env)

} # TXenv





ISsimu.iri.old <- function(time=c(2000,1,1,11,0,0),latitude=69.5864,longitude=19.2272,hmin=50,hmax=1.0e3,sampFreq=1.0e5,experiment=list(code=list(c(1)),IPP=c(10000),baudLength=c(1000)),radarFreq=233e6,flen=1000000,spectrumScale=1e14){
# 
# simulated incoherent scatter radar signal with plasma parameters taken from ionospheric models
# 
# time = c(year,month,day,hour,minute,second)  UT
# latitude and longitude in degrees
# hmax = the maximum height in km
# sampFreq = sampling frequency in Hz
# experiment = list(list(code),c(ipps),c(baudlengths)) # in us
# radarFreq in Hz
# flen  = number of complex samples in a single data file
# spectrumScale = normalisation factor of the spectrum to keep the signal level inside the dynamic range of the data format
# 

  # save the simulation parameters
  simuParam   <- list(time=time,latitude=latitude,longitude=longitude,hmin=hmin,hmax=hmax,sampFreq=sampFreq,
                     experiment=experiment,radarFreq=radarFreq,p_m0=c(30.5,16.0,1.0),flen=flen)
  save(simuParam,file='simuParam.Rdata')


  # the transmission envelope
  experiment$IPP <- floor(experiment$IPP*sampFreq/1e6)
  experiment$baudLength <- max(floor(experiment$baudLength*sampFreq/1e6),1)
  tx          <- TXenv(experiment)
  txlen       <- length(tx)

  # create the frequency axis for IS spectrum calculation
  fmax        <- 2*radarFreq/100000
  nf          <- 2*round(fmax/10)
  freqs       <- c(seq(0,(nf/2)),seq((-nf/2+1),-1))*10.0
  heights     <- seq(hmin,hmax,by=(1.0e6*.1498962)/sampFreq)

  # the spectrum 
  spectrum    <- ISspectrum.iri( time=time , latitude=latitude , longitude=longitude , heights=heights , freq=freqs , fradar=radarFreq )

  # fill in the missing values (they really are zeros!), and normalize the spectrum
  spectrum[is.na(spectrum)] <- 0
  spectrum <- spectrum*spectrumScale

  # zero-padding and square root of the spectrum
  fdiff       <- floor(sampFreq/10)-nf
  nh          <- length(heights)
  nc          <- nf/2
  nf2         <- nf + fdiff
  spectrumsqr <- matrix(0,nrow=nh,ncol=nf2)
  spectrumsqr[,1:nc] <- spectrum[,1:nc]
  spectrumsqr[,(nc+fdiff+1):nf2] <- spectrum[,(nc+1):nf]

  # range in time units
  ranges      <- floor(heights/(299792.458/2)*sampFreq)

  # smallest and largest range
  rmin        <- min(ranges)
  rmax        <- max(ranges)
  for(k in seq(nh)){
    spectrumsqr[k,] <- sqrt(spectrumsqr[k,]) / ranges[k] /  sampFreq
  }

  # matrix for the non-correlating random signal
  sigm        <- matrix(nrow=nh,ncol=nf2)
  # initial random signal
  sigm[,]     <- rnorm(nf2*nh) + 1i*rnorm(nf2*nh)

  # matrix for the correlating signal
  sigcm       <- sigm[,]*0

  # overlap in simulation windows (longest correlation), 10 ms
  overlap     <- floor(sampFreq/1e4)

  # timestamps file
  fid         <- file('timestamps.log','w')
  unixtime    <- makeUnixTime(time) 
  cat(paste('simudata-000001.gdf ',as.character(unixtime),'.0000000',sep=''),'\n',file=fid)
  close(fid)

  # save the parameters to a file before continuing
  save(heights,freqs,spectrum,file='modelspectrum.Rdata')

  # vectors for the radar signal and tx bits
  rsig        <- rep((0+0i),flen)
  rtx         <- rep(F,flen)

  # an infinite loop, break with ctrl-c
  fnum        <- 1 # current file
  snum        <- 1 # current sample index in file
  tnum        <- 1 # current index in the transmission signal (tx) vector

  tx <- tx/2**14
  while(T){

    # create proper correlating signals at all ranges
    for(k in seq(1,nh)) sigcm[k,] <- fft( ( fft(sigm[k,]) * spectrumsqr[k,] ) , inverse=T ) / nf2

    # signal values
    for(k in seq((nf2-overlap))){

      # if the radar is transmitting, simply copy the TX signal
      if(tx[tnum]!=0){
        rsig[snum] <- tx[tnum]
        rtx[snum]  <- TRUE
      }else{
        # proper part of the envelope
        curenv     <- currentEnvelope(tx,tnum,txlen,rmin,rmax) 
        rsig[snum] <- sum(sigcm[,k]*curenv)
        rtx[snum]  <- FALSE
      }

      # increase counters
      snum <- snum + 1
      tnum <- tnum + 1

      # if the data vector is full, write a new data file
      if(snum>flen){
        writeSimuDataFile( fnum=fnum , rsig=rsig , rtx=rtx , flen=flen , fileType='gdf' )
        fnum <- fnum + 1
        snum <- 1
      }

      # tx envelope is cycled
      if(tnum>txlen) tnum <- 1

    }

    # shift the overlapping part of sigm and generate new random signals
    if(overlap>0) sigm[,1:overlap] <- sigm[,(nf2-overlap+1):nf2]
    sigm[,(overlap+1):nf2]         <- rnorm((nf2-overlap)*nh) + 1i*rnorm((nf2-overlap)*nh)

  }

} #ISSimu.iri

Rdata2gdf.dir <- function(ddir='.',prefix='simudata',echoScale=1,txScale=1){

  f <- dir(ddir,pattern=prefix,full.names=TRUE)
  f <- f[grep('Rdata',f)]

  for(fn in f) Rdata2gdf(fn,echoScale,txScale)
  
}

Rdata2gdf <- function(fname,echoScale=1,txScale){
  load(fname)
  nc <- nchar(fname)
  gdfprefix <- substr(fname,1,nc-13)
  gdfnum <- as.numeric(substr(fname,nc-11,nc-6))
  rsig[rtx] <- rsig[rtx]*txScale
  rsig[!rtx] <- rsig[!rtx]*echoScale
  writeSimuDataFile.gdf(fnum=gdfnum,rsig=rsig,rtx=rtx,flen=flen,prefix=gdfprefix,nameadd=NULL)
  invisible()
}

writeSimuDataFile <- function(fnum,rsig,rtx,flen,prefix='simudata',nameadd=NULL,fileType='gdf'){
  if(tolower(fileType[1])=='gdf'){
    writeSimuDataFile.gdf(fnum=fnum,rsig=rsig,rtx=rtx,flen=flen,prefix=prefix,nameadd=nameadd)
  }else if(tolower(fileType[1])=='rdata'){
    writeSimuDataFile.Rdata(fnum=fnum,rsig=rsig,rtx=rtx,flen=flen,prefix=prefix,nameadd=nameadd)
  }else{
    stop(paste('Unkonwn data file type',fileType[1]))
  }

  invisible()
}

writeSimuDataFile.gdf <- function(fnum,rsig,rtx,flen,prefix='simudata',nameadd=NULL){

  # data file name
  fchar <- as.character(fnum)
  nf <- nchar(fchar)
  fname <- paste(prefix,nameadd,'-',substr('000000',1,(6-nf)),fchar,'.gdf',sep='')

  # open the file
  f <- file(fname,'wb')

  # conversion to integer format with tx bits in the lowest bit of the imaginary part
  resig <-as.integer( trunc( Re(rsig)*(2**14) )*2 )
  imsig <-as.integer( trunc( Im(rsig)*(2**14) )*2 ) 

  if(any(abs(resig)>2**14)|any(abs(imsig)>2**14)) warning('Found data values larger than 2**14, the signal will be cut')

  imsig[rtx] <- imsig[rtx]+1

  # proper integer array of the whole data
  risig <- vector(length(2*flen),mode='integer')
  risig[seq(1,(2*flen),by=2)] <- resig
  risig[seq(2,(2*flen),by=2)] <- imsig

  # write the binary data to the file
  writeBin(as.integer(risig),con=f,size=2,endian='big')

  # close the file
  close(f)

  print(fname)

  return()

} # writeSimuDataFile.gdf

writeSimuDataFile.Rdata <- function(fnum,rsig,rtx,flen,prefix='simudata',nameadd=NULL){

  # data file name
  fchar <- as.character(fnum)
  nf <- nchar(fchar)
  fname <- paste(prefix,nameadd,'-',substr('000000',1,(6-nf)),fchar,'.Rdata',sep='')

  # make sure that vector lengths match with flen
  rsig <- rsig[1:flen]
  rtx  <- rtx[1:flen]
  
  save(rsig=rsig,rtx=rtx,flen=flen,file=fname)

  print(fname)

  return()

} # writeSimuDataFile.Rdata

currentEnvelope <- function(tx,tnum,txlen,rmin,rmax){

  si <- tnum - rmin
  li <- tnum - rmax
  ind <- seq(si,li,by=-1)
  while(any(ind<1)) ind[ind<1] <- ind[ind<1]+txlen
  return(tx[ind])

} # currentEnvelope

dateSec <- function(time=c(2000,1,1,0,0,0)){
# 
# conversion from c(year,month,day,hour,minute,seconds) to c(year,seconds)
# 
# 
  # number of days in each month
  mdays <- c(31,28,31,30,31,30,31,31,30,31,30,31)
  mdays[2] <- ifelse(time[1]%%4==0,ifelse(time[1]%%100==0,ifelse(time[1]%%400==0,29,28),29),28)

  # year and seconds from the beginning of the year
  return(c(time[1], (((( sum(mdays[1:(time[2]-1)]) + time[3] - 1 )*24 + time[4] )*60 + time[5])*60 + time[6]) ))

} # dateSec


makeUnixTime <- function(time){
# time = c(year,monh,day,hour,minute,seconds)
# convert a time string into seconds counted from 01.01.1970
#

  ysecs <- 0
  for (k in seq(1970,(time[1]-1))){
      ysecs <- ysecs + dateSec(c(k,12,31,24,0,0))[2]
  }

  return((ysecs + dateSec(time)[2]))

} #makeUnixTime

addRandomNoise <- function(ddir='.',prefix='simudata-',std=1,firstfile=1,lastfile=Inf,maxmiss=10,echoScale=1,txScale=1,fileType=c('Rdata','gdf')){
# 
# Add random noise to simulated radar data
# 
# 
  

  # read the data and add noise
  k <- firstfile
  nmiss <- 0
  repeat{
    if(k>lastfile) break
    if(nmiss>maxmiss) break
    fchar <- as.character(k)
    nf <- nchar(fchar)
    if(tolower(fileType[1])=='rdata'){
      fname <- file.path(ddir,paste('simudata-',substr('000000',1,(6-nf)),fchar,'.Rdata',sep=''))
    }else if(tolower(fileType[1])=='gdf'){
      fname <- file.path(ddir,paste('simudata-',substr('000000',1,(6-nf)),fchar,'.gdf',sep=''))
    }else{
      stop(paste('Unknown file type:',fileType[1]))
    }
    if(file.exists(fname)){
      nmiss <- 0
      if(tolower(fileType[1])=='rdata'){
        load(fname)
        ntx        <- sum(rtx)
        necho      <- flen - ntx
        rsig[!rtx] <- ( rsig[!rtx] + rnorm(necho,0,std) + 1i*rnorm(necho,0,std) ) * echoScale
        rsig[rtx]  <- rsig[rtx] * txScale
        writeSimuDataFile(k,rsig,rtx,flen,nameadd=paste('_noise_',as.character(std),sep=''),fileType='Rdata')
      }else{
        flen  <- file.info(fname)$size/4
        rdata <- readData.gdf(fname,flen)
        necho <- sum(!rdata$idatai)
        rdata$cdata[!rdata$idatai] <- ( rdata$cdata[!rdata$idatai] +  rnorm(necho,0,std) + 1i*rnorm(necho,0,std) ) * echoScale
        writeSimuDataFile(k,rdata$cdata/2**14,rdata$idatai,flen,nameadd=paste('_noise_',as.character(std),sep=''),fileType='gdf')
      }
    }else{
      nmiss <- nmiss + 1
    }
    k <- k+1
  } 
  
  
} # addRandomNoise


DABpowerSpectrum <- function( transmissionMode=1 , f=seq(-2e6,2e6,by=10) , centreFreq=0)
  {
    #
    # Power spectral density of T-DAB transmission according to standard ETS 300 401 (RE/JTC-00DAB-4)
    #
    # INPUT:
    #    transmissionMode  one of the transmission modes (1...4) in the standard
    #    f                 frequency axis
    #
    # OUTPUT:
    #   powerSpectrum      power spectral density at each frequency point
    #
    #

    if(!any(transmissionMode==c(1,2,3,4))) stop('Transmission mode must be 1, 2, 3, or 4.')


    # list of all tx modes in the standard (From Table 1, page 31, in RE/JTC-00DAB-4)
    txmodes <- list(
                    c('nCarriers'=1536,'carrierSpacing'=1000,'symbolDuration'=1246e-6,'guardInterval'=246e-6),
                    c('nCarriers'=384 ,'carrierSpacing'=4000,'symbolDuration'=312e-6 ,'guardInterval'=62e-6 ),
                    c('nCarriers'=192 ,'carrierSpacing'=8000,'symbolDuration'=156e-6 ,'guardInterval'=31e-6 ),
                    c('nCarriers'=768 ,'carrierSpacing'=2000,'symbolDuration'=623e-6 ,'guardInterval'=123e-6)
                    )
    # select the appropriate mode from the list
    txmode <- txmodes[[transmissionMode]]

    # zero-vector for the spectrum
    powerSpectrum <- f[]*0

    # go through all carrier frequencies and add their contribution to the power spectral density
    # (According to the equation in section 15.3, page 190, in RE/JTC-00DAB-4)
    for( k in seq( (-txmode[['nCarriers']]/2) , (txmode[['nCarriers']]/2)) ){
      if(k!=0){
        fk            <- centreFreq +  k*txmode[['carrierSpacing']]
        x             <- pi*(f-fk)*txmode[['symbolDuration']]
        y             <- (sin(x)/x)**2
        y[x==0]       <- 1
        powerSpectrum <- powerSpectrum + y
      }
    }
    return(powerSpectrum)
  }



ISsimu.general <- function(ISspectra,rmin=1,TXenvelope,flen=1000000,fileType=c('Rdata','gdf'),time0=0,timestep=flen/1e6){
# 
# simulated incoherent scatter radar signal with power spectral densities given in ISspectra
#
#
# ISspectra  a matrix with power spectral densities as row-vectors,
#            first row in ISspectra is spectrum at range rmin. The spectra must be at even range intervals, i.e.
#            ISspectrum[2,] contains spectrum from range rmin+1, ISspectrum[3,] from rmin+2, etc.
# rmin       range from which the spectrum in ISspectra[1,] is assumed to be from
# TXenvelope transmission envelope of the experiment
# flen       number of samples in each output file
# fileType   output file format, either 'gdf' (16-bit signed integers) or 'Rdata' (64-bit floats)
# 

  # number of frequencies and number of ranges
  sdims <- dim(ISspectra)
  nh <- sdims[1]
  nf <- sdims[2]
  rmax <- rmin + nh - 1

  # squareroot of the spectra
  ISspectraSqr <- t(sqrt(ISspectra))

  # matrix for the non-correlating random signal
  sigm        <- matrix(ncol=nh,nrow=nf)
  
  # initial random signal
  sigm[,]     <- rnorm(nf*nh) + 1i*rnorm(nf*nh)

  # matrix for the correlating signal
  sigcm       <- sigm[,]*0

  # overlap in simulation windows (longest correlation), 10 ms
  overlap     <- floor(nf/1e4)

  # save the parameters to a file before continuing
  save(ISspectra,rmin,TXenvelope,flen,file='ISsimuParam.Rdata')

  # vectors for the radar signal and tx bits
  rsig        <- rep((0+0i),flen)
  rtx         <- rep(F,flen)

  # an infinite loop, break with ctrl-c
  fnum        <- 1 # current file
  snum        <- 1 # current sample index in file
  tnum        <- 1 # current index in the transmission signal (tx) vector

  tx <- TXenvelope/2**14
  txlen <- length(tx)
  repeat{


    gc()
    
    # create proper correlating signals at all ranges
#    for(k in seq(1,nh)) sigcm[k,] <- fft( ( fft(sigm[k,]) * ISspectraSqr[k,] ) , inverse=T ) / nf
    sigcm <- mvfft( ( mvfft(sigm) * ISspectraSqr ) ,inverse=TRUE ) / nf

    # signal values
    for(k in seq((nf-overlap))){

      # if the radar is transmitting, simply copy the TX signal
      if(tx[tnum]!=0){
        rsig[snum] <- tx[tnum]
        rtx[snum]  <- TRUE
      }else{
        # proper part of the envelope
        curenv     <- currentEnvelope(tx,tnum,txlen,rmin,rmax) 
        rsig[snum] <- sum(sigcm[k,]*curenv)
        rtx[snum]  <- FALSE
      }

      # increase counters
      snum <- snum + 1
      tnum <- tnum + 1

      # if the data vector is full, write a new data file
      if(snum>flen){
        writeSimuDataFile( fnum=fnum , rsig=rsig , rtx=rtx , flen=flen , fileType=fileType[1] )
        writeTimestampsFile.gdf(prefix='simudata-',extension='.gdf',nstart=1,nend=fnum,times=seq(fnum)*timestep+time0,fname='timestamps.log')
        fnum <- fnum + 1
        snum <- 1
      }

      # tx envelope is cycled
      if(tnum>txlen) tnum <- 1

    }

    # shift the overlapping part of sigm and generate new random signals
    if(overlap>0) sigm[1:overlap,] <- sigm[(nf-overlap+1):nf,]
    sigm[(max(overlap,0)+1):nf,]          <- rnorm((nf-max(overlap,0))*nh) + 1i*rnorm((nf-max(overlap,0))*nh)

  }

} #ISsimu.general






ISsimu.iri <- function(time=c(2000,1,1,11,0,0),latitude=69.5864,longitude=19.2272,hmin=50,hmax=1.0e3,sampFreq=1.0e5,experiment=list(code=list(c(1)),IPP=c(10000),baudLength=c(1000)),radarFreq=233e6,flen=1000000,spectrumScale=1e30,fileType=c('Rdata','gdf')){
# 
# simulated incoherent scatter radar signal with plasma parameters taken from ionospheric models
# 
# time = c(year,month,day,hour,minute,second)  UT
# latitude and longitude in degrees
# hmax = the maximum height in km
# sampFreq = sampling frequency in Hz
# experiment = list(list(code),c(ipps),c(baudlengths)) # in us
# radarFreq in Hz
# flen  = number of complex samples in a single data file
# spectrumScale = normalisation factor of the spectrum to keep the signal level inside the dynamic range of the data format
# 

  # save the simulation parameters
  simuParam   <- list(time=time,latitude=latitude,longitude=longitude,hmin=hmin,hmax=hmax,sampFreq=sampFreq,
                     experiment=experiment,radarFreq=radarFreq,p_m0=c(30.5,16.0,1.0),flen=flen,spectrumScale=spectrumScale)
  save(simuParam,file='ISsimu.iriParam.Rdata')


  # the transmission envelope
  experiment$IPP <- floor(experiment$IPP*sampFreq/1e6)
  experiment$baudLength <- max(floor(experiment$baudLength*sampFreq/1e6),1)
  tx          <- TXenv(experiment)
  txlen       <- length(tx)

  # create the frequency axis for IS spectrum calculation
  fmax        <- min(8*radarFreq/100000,sampFreq/2)
  nf          <- 2*round(fmax/10)
  freqs       <- c(seq(0,(nf/2)),seq((-nf/2+1),-1))*10.0
  heights     <- seq(hmin,hmax,by=(1.0e6*.149896229)/sampFreq)

  # range in time units
  ranges      <- floor(heights/(299792.458/2)*sampFreq)

  # remove 0-heights
  heights     <- heights[ranges>0]
  ranges      <- ranges[ranges>0]

  # the spectrum 
  spectrum    <- ISspectrum.iri( time=time , latitude=latitude , longitude=longitude , heights=heights , freq=freqs , fradar=radarFreq , savePlasmaParams=TRUE)

  # fill in the missing values (they really are zeros!), and normalize the spectrum
  spectrum[is.na(spectrum)] <- 0
  spectrum <- spectrum*spectrumScale
  # there will be some negative values due to numerical inaccuracies, set them to zero
  spectrum[spectrum<0] <- 0

  # smallest and largest range
  rmin        <- min(ranges)
  rmax        <- max(ranges)

  # zero-padding and shifting the zero-frequency to the first column
  fdiff       <- floor(sampFreq/10)-nf
  nh          <- length(heights)
  nc          <- nf/2
  nf2         <- nf + fdiff
  spectrump <- matrix(0,nrow=nh,ncol=nf2)
  spectrump[,1:nc] <- spectrum[,1:nc]
  spectrump[,(nc+fdiff+1):nf2] <- spectrum[,(nc+1):nf]

  # scale with range squared
  for(k in seq(nh)){
    spectrump[k,] <- spectrump[k,] / (ranges[k] *  sampFreq )**2
  }


#  # timestamps file
#  fid         <- file('timestamps.log','w')
#  unixtime    <- makeUnixTime(time) 
#  cat(paste('simudata-000001.gdf ',as.character(unixtime),'.0000000',sep=''),'\n',file=fid)
#  close(fid)


  ISsimu.general( ISspectra=spectrump , rmin=rmin , TXenvelope=tx , flen=flen ,fileType=fileType[1] , time0=unixtime , timestep=flen/sampFreq)

} #ISsimu.iri



IQsample <- function(rsig,t=seq(length(rsig)),cfreq,nfilter){

  # IQ-sampling of the real signal rsig, with centre frequency cfreq and filter length (plus decimation)
  # of nfilter samples

  iqsig <- rsig * exp(1i*2*pi*cfreq*t)
  boxcar <- t[]*0
  boxcar[1:nfilter]<-1/sqrt(nfilter)
  iqfilter <- fft( ( fft(iqsig) * fft(boxcar) )[1:floor(length(t)/nfilter)] , inverse=TRUE ) / length(rsig)

  return(iqfilter)
  
}

IQ2real <- function(iqsig,t=seq(length(iqsig)),cfreq,nfilter){

  # an imperfect reconstruction of the original real input signal from an IQ signal
  #
  # t is  a time vector for the *output* data vector!!
  #

  ns <- length(iqsig)
  
  boxcar <- rfft <- rep(0,ns*nfilter)
  
  boxcar[1:nfilter] <- 1*sqrt(nfilter)
  
  iqfft <- fft(iqsig)
  
  rfft[1:ns] <- iqfft

  boxfft <- fft(boxcar)
  # replace zeros with ones, which is possible because we assumed a narrow-band signal
  boxfft[abs(boxfft)<1e-10] <- 1

  # a signal with one spectral peak at correct place
  rsig <- fft(  rfft / boxfft , inverse=TRUE )   * exp( -1i*2*pi*cfreq*t )

  # add the mirrored spectrum peak and return
  frsig <- fft(rsig)
  frsig[2:(ns*nfilter)] <- frsig[2:(ns*nfilter)] + Conj(rev(frsig[2:(ns*nfilter)]))
  rsig <- fft( frsig , inverse=TRUE) / ns**2

  return(Re(rsig))
  

}


correlatedNoise <- function(spectrum,overlap,n,rnorm0=(rnorm(length(spectrum))+1i*rnorm(length(spectrum)))/sqrt(2)){
  # Random noise signal with power spectral density given in spectrum
  #
  # INPUT:
  #   spectrum   a real vector of the power spectral density
  #   overlap    hard to explain, see the code... overlap < length(spectrum)!
  #   n          number of points to simulate
  #   rnorm0     a random non-correlated signal. This can be picked from the output of a previous function call
  #              to guarantee that outputs from consequtive calls have proper correlation properties.
  #
  # I. Virtanen 2012
  #

  ns <- length(spectrum)
  if(length(rnorm0)!=ns){
    warning("Lengths of spectrum and rnorm0 do not match, replacing rnorm0 with a random sequence")
    rnorm0 <- (rnorm(ns)+1i*rnorm(ns))/sqrt(2)
  }

  outlist <- list()

  outlist$signal <- rep((0+0i),n)

  specsqr <- sqrt(spectrum)
  
  k <- 1
  while(k<n){
    sigcor <- fft( fft(rnorm0) * specsqr , inverse=TRUE) /sqrt(ns)
    m <- min((n-k),(ns-overlap))
    outlist$signal[k:(k+m-1)] <- sigcor[1:m]
    rnorm0[1:(ns-m)] <- rnorm0[(m+1):ns]
    rnorm0[(ns-m+1):ns] <- (rnorm(m)+1i*rnorm(m))/sqrt(2)
    k <- k+ns-overlap
  }

  outlist$rnorm0=rnorm0

  return(outlist)
  
}
