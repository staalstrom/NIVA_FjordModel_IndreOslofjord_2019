#================================================================
# Who:  Andre Staalstrom (ans@niva.no)
# When: 3.10.2019
# What: Tar modellresulatater fra NIVA Fjordmodell som er lagret
#       som CSV fil, og legger dette inn i NetCDF filer
#================================================================
# Disse pakkene må installeres (gjøres bare en gang)
#----------------------------------------------------------------
#install.packages("chron")
#install.packages("ncdf4")
#install.packages("RColorBrewer")
# Og så må bibliotekene lastes inn (gjøres hver gang)
library(chron)
library(RColorBrewer)
library(lattice)
library(ncdf4)
# Set working directory
workdir <- "C:\\Users\\ANS\\OneDrive - NIVA\\WinProgProjects\\NIVA_FjordModel_IndreOslofjord_2019\\Compile_and_Run\\"
setwd(workdir)
getwd()
#================================================================
# Setup
#================================================================
# input
# xls <- "EUTRO.xlsx"  # not used
# output file names
nc    <- "nivafm_b1_BF.nc"
nc[2] <- "nivafm_b2_LY.nc"
nc[3] <- "nivafm_b3_VF.nc"
nc[4] <- "nivafm_b4_BB.nc"
nc[5] <- "nivafm_b5_BL.nc"
nc[6] <- "nivafm_b6_OH.nc"
nc[7] <- "nivafm_b7_BR.nc"
nc[8] <- "nivafm_b8_HF.nc"
nc[9] <- "nivafm_b9_SN.nc"
# # CTD stations in each basin
# lat    <- 59.72577   # Gp1
# lat[2] <- 59.84048   # Cj1
# lat[3] <- 59.89880   # Aq2
# lat[4] <- 59.82716   # Dm1
# lat[5] <- 59.88083   # Bl4
# lat[6] <- 59.87285   # Cq1
# lat[7] <- 59.88047   # Bn1
# lat[8] <- 59.78630   # Ep1
# lat[9] <- 59.81500   # Dk1
 
# lon    <- 10.72645    # Gp1
# lon[2] <- 10.50698    # Cj1
# lon[3] <- 10.74487    # Aq2
# lon[4] <- 10.61727    # Dm1
# lon[5] <- 10.56883    # Bl4
# lon[6] <- 10.73667    # Cq1
# lon[7] <- 10.64673    # Bn1
# lon[8] <- 10.72378    # Ep1
# lon[9] <- 10.56938    # Dk1
#================================================================
# Read CSV into R
#================================================================
NUM <- read.csv(file="EUTRO.csv", header=FALSE, sep=";", dec = ".")

FILETIME <- NUM[1,2]
iBASSENG <- NUM[,3]
iLAG     <- NUM[,4]
iOBSNUM  <- NUM[,5]
iT_AAR   <- NUM[,6]
iSESONG  <- NUM[,7]
iDYP_M   <- NUM[,8]
iTEMP    <- NUM[,9]
iSAL     <- NUM[,10]
iDENSI   <- NUM[,11]
iOXYG    <- NUM[,12]
iTOTN    <- NUM[,13]
iNO3     <- NUM[,14]
iNH4     <- NUM[,15]
iTOTP    <- NUM[,16]
iPO4     <- NUM[,17]
iSIO2    <- NUM[,18]
iCZOO    <- NUM[,19]
iBACT    <- NUM[,20]
iCFYT1   <- NUM[,21]
iNFYT1   <- NUM[,22]
iPFYT1   <- NUM[,23]
iSFYT    <- NUM[,24]
iCFYT2   <- NUM[,25]
iNFYT2   <- NUM[,26]
iPFYT2   <- NUM[,27]
iODM     <- NUM[,28]
iDOC     <- NUM[,29]
iTOTC    <- NUM[,30]
iPARTC   <- NUM[,31]
iPARTN   <- NUM[,32]
iPARTP   <- NUM[,33]
iCDFLUX  <- NUM[,34]
iNDFLUX  <- NUM[,35]
iPDFLUX  <- NUM[,36]
iSDFLUX  <- NUM[,37]
iCSED    <- NUM[,38]
iNSED    <- NUM[,39]
iPSED    <- NUM[,40]
iSSED    <- NUM[,41]
iASED    <- NUM[,42]
iPADS    <- NUM[,43]
iSPP     <- NUM[,44]
iSPPSED  <- NUM[,45]
iSPPDV   <- NUM[,46]
iSPPSEDDV<- NUM[,47]
iCHL1    <- NUM[,48]
iCHL2    <- NUM[,49]
iC1       <- NUM[,50]
iTOTN       <- NUM[,51]
iTOTP       <- NUM[,52]



I <- length(iTEMP)

maxib <- max(iBASSENG)
#================================================================
# Prep data for NetCDF file
#================================================================
# LOOP THROUGH BASINS
for(ib in 1:maxib) {
  cat("BASIN: ",ib,"\n")


#ib<-1
#cat("BASIN: ",ib,"\n")


# Find unique times
T_AAR <- iT_AAR[iBASSENG==ib & iLAG==2]
# T_AAR = T_AAR(2:end); 
N <- length(T_AAR)
# Find number of layers in basin
# KK=max(iLAG(iBASSENG==ib));
# KK=max(iLAG(iBASSENG==ib & iOBSNUM==1));
# Unique LAGs
# LAG = 1:KK;
LAG <- iLAG[iBASSENG==ib & iOBSNUM==1]
KK <- length(LAG)

depth<-iDYP_M[iBASSENG==ib & iOBSNUM==1]

#================================================================
# Write time vector to empty NetCDF file (this expands the vars)
#================================================================
T_DAY <- T_AAR*365
# time = T_DAY + datenum(2000,1,1);
time <- T_DAY*24*3600

# set path and filename
#ncpath <- workdir
#ncname <- nc[ib]  
#ncfname <- paste(ncpath, ncname, sep="")

#================================================================
# Make dimensions
#================================================================
zdim <- ncdim_def( 'depth', 'm', -depth)
tdim <- ncdim_def( 'time', 'seconds since 2000-01-01', time )

#================================================================
# Make variable
#================================================================
mv <- 1.e36 # missing value
vTEMP <- ncvar_def( 'TEMP', 'Celsius', list(zdim,tdim), mv ,longname='Temperature')
vSAL  <- ncvar_def( 'SAL', 'PSU', list(zdim,tdim), mv ,longname='Salinity')
vDENS <- ncvar_def( 'DENS', 'kg/m3', list(zdim,tdim), mv ,longname='Density (sigma-t)')
vDOC  <- ncvar_def( 'DOC', 'micro g/L', list(zdim,tdim), mv ,longname='Dissolved organic matter')
vOXYG  <- ncvar_def( 'OXYG', 'ml O2/L', list(zdim,tdim), mv ,longname='Oxygen volume concentration')
vTOTP  <- ncvar_def( 'TOTP', 'micro g P/L', list(zdim,tdim), mv ,longname='Totale phosphorous')
vPO4  <- ncvar_def( 'PO4', 'micro g P/L', list(zdim,tdim), mv ,longname='Phosphate')
vTOTN  <- ncvar_def( 'TOTN', 'micro g N/L', list(zdim,tdim), mv ,longname='Totale nitrogen')
vNO3  <- ncvar_def( 'NO3', 'micro g N/L', list(zdim,tdim), mv ,longname='Sum of nitrate and nitrite')
vNH4  <- ncvar_def( 'NH4', 'micro g N/L', list(zdim,tdim), mv ,longname='Ammonium')
vSIO2  <- ncvar_def( 'SIO2', 'micro g SIO2/L', list(zdim,tdim), mv ,longname='Silicate')
vCDET  <- ncvar_def( 'CDET', 'micro g C/L', list(zdim,tdim), mv ,longname='Particulate carbon')
vNDET  <- ncvar_def( 'NDET', 'micro g N/L', list(zdim,tdim), mv ,longname='Particulate nitrogen')
vPDET  <- ncvar_def( 'PDET', 'micro g P/L', list(zdim,tdim), mv ,longname='Particulate phosphorous')
vCHL1  <- ncvar_def( 'CHL1', 'micro g/L', list(zdim,tdim), mv ,longname='Chlorophyll a in diatoms')
vCFYT1  <- ncvar_def( 'CFYT1', 'micro g/L', list(zdim,tdim), mv ,longname='Carbon in diatoms')
vNFYT1  <- ncvar_def( 'NFYT1', 'micro g/L', list(zdim,tdim), mv ,longname='Nitrogen in diatoms')
vPFYT1  <- ncvar_def( 'PFYT1', 'micro g/L', list(zdim,tdim), mv ,longname='Phosphourous in diatoms')
vSFYT1  <- ncvar_def( 'SFYT1', 'micro g/L', list(zdim,tdim), mv ,longname='Silicate in diatoms')
vCHL2  <- ncvar_def( 'CHL2', 'micro g/L', list(zdim,tdim), mv ,longname='Chlorophyll a in other phytoplankon')
vCFYT2  <- ncvar_def( 'CFYT2', 'micro g/L', list(zdim,tdim), mv ,longname='Carbon in other phytoplankon')
vNFYT2  <- ncvar_def( 'NFYT2', 'micro g/L', list(zdim,tdim), mv ,longname='Nitrogen in other phytoplankon')
vPFYT2  <- ncvar_def( 'PFYT2', 'micro g/L', list(zdim,tdim), mv ,longname='Phosphourous in other phytoplankon')

vC1  <- ncvar_def( 'C1', '-', list(zdim,tdim), mv ,longname='Passive dye')
vCFYT  <- ncvar_def( 'CFYT', 'micro g/L', list(zdim,tdim), mv ,longname='Carbon in plankton')
vCHL  <- ncvar_def( 'CHL', 'micro g/L', list(zdim,tdim), mv ,longname='Chlorophyll a in plankton')

vCZOO  <- ncvar_def( 'CZOO', 'micro g C/L', list(zdim,tdim), mv ,longname='Zooplankton')
vBACT  <- ncvar_def( 'BACT', 'micro g C/L', list(zdim,tdim), mv ,longname='Marine bacteria')


#================================================================
# Create a new NetCDF file
#================================================================
ncin <- nc_create( nc[ib], list(vTEMP,vSAL,vDENS,vDOC,vOXYG,vTOTP,vPO4,vTOTN,vNO3,vNH4,vSIO2,vCDET,vNDET,vPDET,vCHL1,vCFYT1,vNFYT1,vPFYT1,vSFYT1,vCHL2,vCFYT2,vNFYT2,vPFYT2,vC1,vCFYT,vCHL,vCZOO,vBACT))
print(ncin)

#================================================================
# Define empty vars
#================================================================
TEMP  <-matrix(NA,ncol=N,nrow=KK)
SAL   <-matrix(NA,ncol=N,nrow=KK) 
DENS  <-matrix(NA,ncol=N,nrow=KK)
DOC   <-matrix(NA,ncol=N,nrow=KK)
OXYG  <-matrix(NA,ncol=N,nrow=KK) 
TOTP  <-matrix(NA,ncol=N,nrow=KK)
PO4   <-matrix(NA,ncol=N,nrow=KK) 
TOTN  <-matrix(NA,ncol=N,nrow=KK)
NO3   <-matrix(NA,ncol=N,nrow=KK)
NH4   <-matrix(NA,ncol=N,nrow=KK)
SIO2  <-matrix(NA,ncol=N,nrow=KK)

CDET  <-matrix(NA,ncol=N,nrow=KK)
NDET  <-matrix(NA,ncol=N,nrow=KK)
PDET  <-matrix(NA,ncol=N,nrow=KK)

CHL1  <-matrix(NA,ncol=N,nrow=KK)
CFYT1 <-matrix(NA,ncol=N,nrow=KK)
NFYT1 <-matrix(NA,ncol=N,nrow=KK)
PFYT1 <-matrix(NA,ncol=N,nrow=KK)
SFYT1 <-matrix(NA,ncol=N,nrow=KK)
CHL2  <-matrix(NA,ncol=N,nrow=KK)
CFYT2 <-matrix(NA,ncol=N,nrow=KK)
NFYT2 <-matrix(NA,ncol=N,nrow=KK)
PFYT2 <-matrix(NA,ncol=N,nrow=KK)

C1    <-matrix(NA,ncol=N,nrow=KK)

CHL  <-matrix(NA,ncol=N,nrow=KK)
CFYT <-matrix(NA,ncol=N,nrow=KK)

CZOO  <-matrix(NA,ncol=N,nrow=KK)
BACT <-matrix(NA,ncol=N,nrow=KK)

#================================================================
# Loop through it all
#================================================================
for(i in 1:I) {
    if (iBASSENG[i]==ib){
# Find the place to put the data
		nput=which(iT_AAR[i]==T_AAR)
        kput=which(LAG==iLAG[i]);
		cat(nput," ",kput,"\n")
		
        TEMP[kput,nput]=iTEMP[i]
        SAL[kput,nput]=iSAL[i]
        DENS[kput,nput]=iDENSI[i]
        DOC[kput,nput]=iDOC[i]
        OXYG[kput,nput]=iOXYG[i]
        TOTP[kput,nput]=iTOTP[i]
        PO4[kput,nput]=iPO4[i]
        TOTN[kput,nput]=iTOTN[i]
        NO3[kput,nput]=iNO3[i]
        NH4[kput,nput]=iNH4[i]
        SIO2[kput,nput]=iSIO2[i]
        CDET[kput,nput]=iPARTC[i]
        NDET[kput,nput]=iPARTN[i]
        PDET[kput,nput]=iPARTP[i]
        CHL1[kput,nput]=iCHL1[i]
        CFYT1[kput,nput]=iCFYT1[i]
        NFYT1[kput,nput]=iNFYT1[i]
        PFYT1[kput,nput]=iPFYT1[i]
        SFYT1[kput,nput]=iSFYT[i]
        CHL2[kput,nput]=iCHL2[i]
        CFYT2[kput,nput]=iCFYT2[i]
        NFYT2[kput,nput]=iNFYT2[i]
        PFYT2[kput,nput]=iPFYT2[i]
        C1[kput,nput]=iC1[i]
		CHL[kput,nput]=iCHL1[i]+iCHL2[i]
        CFYT[kput,nput]=iCFYT1[i]+iCFYT2[i]
		CZOO[kput,nput]=iCZOO[i]
        BACT[kput,nput]=iBACT[i]
        TOTN[kput,nput]=iTOTN[i]
        TOTP[kput,nput]=iTOTP[i]
    }
}
#TEMP(TEMP>1e35)=NaN;
#================================================================
# Write to NetCDF file 
#================================================================
ncvar_put( ncin, vTEMP, TEMP,  start=c(1,1), count=c(KK,N))
ncvar_put( ncin, vSAL,  SAL,   start=c(1,1), count=c(KK,N))
ncvar_put( ncin, vDENS, DENS,  start=c(1,1), count=c(KK,N)) 
ncvar_put( ncin, vDOC,  DOC,   start=c(1,1), count=c(KK,N))
ncvar_put( ncin, vOXYG, OXYG,  start=c(1,1), count=c(KK,N))
ncvar_put( ncin, vTOTP, TOTP,  start=c(1,1), count=c(KK,N))
ncvar_put( ncin, vPO4,  PO4,   start=c(1,1), count=c(KK,N))
ncvar_put( ncin, vTOTN, TOTN,  start=c(1,1), count=c(KK,N))
ncvar_put( ncin, vNO3,  NO3,   start=c(1,1), count=c(KK,N))
ncvar_put( ncin, vNH4,  NH4,   start=c(1,1), count=c(KK,N))
ncvar_put( ncin, vSIO2, SIO2,  start=c(1,1), count=c(KK,N))
ncvar_put( ncin, vCDET, CDET,  start=c(1,1), count=c(KK,N))
ncvar_put( ncin, vNDET, NDET,  start=c(1,1), count=c(KK,N))
ncvar_put( ncin, vPDET, PDET,  start=c(1,1), count=c(KK,N))
ncvar_put( ncin, vCHL1, CHL1,  start=c(1,1), count=c(KK,N))
ncvar_put( ncin, vCFYT1,CFYT1, start=c(1,1), count=c(KK,N))
ncvar_put( ncin, vNFYT1,NFYT1, start=c(1,1), count=c(KK,N))
ncvar_put( ncin, vPFYT1,PFYT1, start=c(1,1), count=c(KK,N))
ncvar_put( ncin, vSFYT1,SFYT1, start=c(1,1), count=c(KK,N))
ncvar_put( ncin, vCHL2, CHL2,  start=c(1,1), count=c(KK,N))
ncvar_put( ncin, vCFYT2,CFYT2, start=c(1,1), count=c(KK,N))
ncvar_put( ncin, vNFYT2,NFYT2, start=c(1,1), count=c(KK,N))
ncvar_put( ncin, vPFYT2,PFYT2, start=c(1,1), count=c(KK,N))
ncvar_put( ncin, vC1,   C1,    start=c(1,1), count=c(KK,N))
ncvar_put( ncin, vCFYT, CFYT,  start=c(1,1), count=c(KK,N))
ncvar_put( ncin, vCHL,  CHL,   start=c(1,1), count=c(KK,N))
ncvar_put( ncin, vCZOO, CZOO,  start=c(1,1), count=c(KK,N))
ncvar_put( ncin, vBACT, BACT,  start=c(1,1), count=c(KK,N))

#================================================================
# Close NetCDF file 
#================================================================
nc_close( ncin )

# quick map
#image(time,depth,TEMP, col=rev(brewer.pal(10,"RdBu")))
}