# 
# Incoherent scatter radar simulator 
# 
# 
# I. Virtanen 2010, 2012
#


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
# I. Virtanen 2010
# 

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

addRandomNoise <- function(ddir='.',prefix='simudata-',std=1,firstfile=1,lastfile=Inf,maxmiss=10,echoScale=1){
# 
# Add random noise to simulated radar data
# 
# 
  

  # read the data and add noise
  k <- firstfile
  nmiss <- 0
  while(T){
    if(k>lastfile) break
    if(nmiss>maxmiss) break
    fchar <- as.character(k)
    nf <- nchar(fchar)
    fname <- file.path(ddir,paste('simudata-',substr('000000',1,(6-nf)),fchar,'.gdf',sep=''))
    if(file.exists(fname)){
      nmiss <- 0
      flen  <- file.info(fname)$size/4
      rdata <- readRawData(fname,flen)
      necho <- sum(!rdata$TX)
      rdata$data[!rdata$TX] <- (rdata$data[!rdata$TX] +  rnorm(necho,0,std) + 1i*rnorm(necho,0,std))/echoScale
      writeSimuDataFile(k,rdata$data/2**14,rdata$TX,flen,nameadd=paste('_noise_',as.character(std),sep=''))
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



ISsimu.general <- function(ISspectra,rmin=1,TXenvelope,flen=1000000,fileType=c('Rdata','gdf')){
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
  ISspectraSqr <- sqrt(ISspectra)

  # matrix for the non-correlating random signal
  sigm        <- matrix(nrow=nh,ncol=nf)
  
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

    # create proper correlating signals at all ranges
    for(k in seq(1,nh)) sigcm[k,] <- fft( ( fft(sigm[k,]) * ISspectraSqr[k,] ) , inverse=T ) / nf

    # signal values
    for(k in seq((nf-overlap))){

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
        writeSimuDataFile( fnum=fnum , rsig=rsig , rtx=rtx , flen=flen , fileType=fileType[1] )
        fnum <- fnum + 1
        snum <- 1
      }

      # tx envelope is cycled
      if(tnum>txlen) tnum <- 1

    }

    # shift the overlapping part of sigm and generate new random signals
    if(overlap>0) sigm[,1:overlap] <- sigm[,(nf-overlap+1):nf]
    sigm[,(overlap+1):nf]          <- rnorm((nf-overlap)*nh) + 1i*rnorm((nf-overlap)*nh)

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
                     experiment=experiment,radarFreq=radarFreq,p_m0=c(30.5,16.0,1.0),flen=flen)
  save(simuParam,file='ISsimu.iri2Param.Rdata')


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


  # timestamps file
  fid         <- file('timestamps.log','w')
  unixtime    <- makeUnixTime(time) 
  cat(paste('simudata-000001.gdf ',as.character(unixtime),'.0000000',sep=''),'\n',file=fid)
  close(fid)


  ISsimu.general( ISspectra=spectrump , rmin=rmin , TXenvelope=tx , flen=flen ,fileType=fileType[1] )

} #ISsimu.iri



