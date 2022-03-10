# DxO PureRAW. Noise reduction based on neural networks
# www.overfitting.net
# https://www.overfitting.net/2022/03/dxo-pureraw-reduccion-de-ruido-con.html

library(tiff)

# READ RAW DATA

# RAW integer extraction using DCRAW:
# dcraw -v -r 1 1 1 1 -o 0 -k 512 -S 16383 -4 -T iso25600.dng
img1=readTIFF("iso25600_trunc.tiff", native=F, convert=F)

# dcraw -v -r 1 1 1 1 -o 0 -k 0 -S 16383 -4 -T iso25600.dng
img1=readTIFF("iso25600_notrunc.tiff", native=F, convert=F)
BLACK=512
SAT=16383
MAX=max(img1)
img1=img1-BLACK/SAT  # linearize to 0..1 preserving negative values
img1=img1*MAX/max(img1)

# dcraw -v -r 1 1 1 1 -o 0 -4 -T iso25600dxo.dng
img2=readTIFF("iso25600dxo.tiff", native=F, convert=F)

# dcraw -v -r 1 1 1 1 -o 0 -4 -T iso25600topaz.dng
# img2=readTIFF("iso25600topaz.tiff", native=F, convert=F)

img1=img1[13:4012, 13:6012,]  # crop area (13,13)-(6012,4012)
img1=img1[1800:2119, 886:4996, 2]  # crop G patches
img2=img2[1800:2119, 886:4996, 2]  # crop G patches


# S/N values
NPATCHES=24
S1=array(0,NPATCHES)
N1=S1
S2=S1
N2=S1

# Loop 24 patches
ALTO=nrow(img1)
ANCHO=ncol(img1)/NPATCHES
OFFX=20
OFFY=5

par(mfrow=c(4,6))
BREAKS=90
for (j in 1:NPATCHES) {
    i=which(row(img1)>=OFFY & row(img1)<=ALTO-OFFY &
            col(img1)>=ANCHO*(j-1)+OFFX & col(img1)<=ANCHO*j-OFFX)

    # Plot histograms before/after
    xmin=min(min(img1[i]), min(img2[i]))
    xmax=max(max(img1[i]), max(img2[i]))
    xrange=seq(from=xmin, to=xmax, length.out=BREAKS)
    
    h1=hist(img1[i], breaks=xrange, plot=FALSE)
    h2=hist(img2[i], breaks=xrange, plot=FALSE)
    
    ymax=max(h1$counts, h2$counts)
    h1$counts=h1$counts/ymax  # normalize histograms to 0..1
    h2$counts=h2$counts/ymax
    
    plot(h1, main=paste0('patch ',j), breaks=xrange, ylim=c(0,1),
         xlab='', ylab='', axes=FALSE,
         col=rgb(0.8, 0.8, 0.8, 1/1), border=rgb(0.8, 0.8, 0.8,0))
    plot(h2, breaks=xrange,
         xlab='', ylab='',
         col=rgb(1,0,0,0.3), border=rgb(1,0,0,0), add=TRUE)
    axis(side=1)  # only x axis
    abline(v=0)
    
    # Signal/noise calculations
    S1[j]=mean(img1[i])  # S=mean
    S2[j]=mean(img2[i])
    N1[j]=var(img1[i])^0.5  # N=stdev
    N2[j]=var(img2[i])^0.5
}
dev.off()


# Check S1 vs S2: S2 slightly higher than S1 (we assume S1 correct)
plot(log2(S1), log2(S2), xlim=c(-7,0), ylim=c(-7,0), col='red',
     main='S1 vs S2',
     xlab='S1 RAW exposure (EV)', ylab='S2 RAW exposure (EV)')
lines(c(-7,0), c(-7,0), col='gray')


# SNR cuves in dB
plot(log2(S1), 20*log10(S2/N2), xlim=c(-6,0), ylim=c(0,30),
     main='DxO PureRAW SNR enhancement',
     xlab='RAW exposure (EV)', ylab='SNR (dB)')
lines(log2(S1), 20*log10(S1/N1), col='red')
abline(h=12, v=0, lty=2)

# SNR curves in EV
plot(log2(S1), log2(S2/N2), xlim=c(-6,0), ylim=c(0,5),
     main='DxO PureRAW SNR enhancement',
     xlab='RAW exposure (EV)', ylab='SNR (EV)')
lines(log2(S1), log2(S1/N1), col='red')
# abline(h=0:5, v=-9:0, col='gray', lty=2)
abline(h=2, lty=2)


# SNR gain in dB
plot(log2(S1), 20*log10((S2/N2)/(S1/N1)), xlim=c(-6,0), ylim=c(0,10),
     main='DxO PureRAW SNR enhancement',
     xlab='RAW exposure (EV)', ylab='DR gain (dB)', col='red')
abline(h=mean(20*log10((S2/N2)/(S1/N1))), lty=2)

# SNR gain in EV
plot(log2(S1), log2((S2/N2)/(S1/N1)), xlim=c(-6,0), ylim=c(0,1.5),
     main='DxO PureRAW SNR enhancement',
     xlab='RAW exposure (EV)', ylab='DR gain (EV)', col='red')
abline(h=mean(log2((S2/N2)/(S1/N1))), lty=2)
