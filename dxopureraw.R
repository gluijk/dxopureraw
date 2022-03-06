# DxO PureRAW. Noise reduction based on Neural Networks
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
img1=img1-BLACK/SAT
img1=img1/max(img1)

# dcraw -v -r 1 1 1 1 -o 0 -4 -T iso25600dxo.dng
img2=readTIFF("iso25600dxo.tiff", native=F, convert=F)

# dcraw -v -r 1 1 1 1 -o 0 -4 -T iso25600topaz.dng
# img2=readTIFF("iso25600topaz.tiff", native=F, convert=F)


img1=img1[13:4012, 13:6012,]  # crop area (13,13)-(6012,4012)
img1=img1[1800:2119, 886:4996, 2]  # crop G patches
img2=img2[1800:2119, 886:4996, 2]  # crop G patches

ALTO=nrow(img1)
ANCHO=ncol(img1)/24
OFFX=20
OFFY=5

# mark 24 patches
S1=array(0,24)
N1=S1
S2=S1
N2=S1
Ratio1=S1
Ratio2=S1

# Loop 24 patches
for (j in 1:24) {
    i=which(row(img1)>=OFFY & row(img1)<=ALTO-OFFY &
            col(img1)>=ANCHO*(j-1)+OFFX & col(img1)<=ANCHO*j-OFFX)
    hist(img1[i], breaks=100)
    S1[j]=mean(img1[i])  # S=mean
    S2[j]=mean(img2[i])
    N1[j]=var(img1[i])^0.5  # N=stdev
    N2[j]=var(img2[i])^0.5
    tmp=img1[img1==0]
    Ratio1[j]=length(tmp)/length(img1[i])  # check black clipping
    tmp=img2[img2==0]
    Ratio2[j]=length(tmp)/length(img2[i])  # check black clipping
}


# Check S1 vs S2
plot(S1, S2, type='b', col='red')

plot(S1/S2, ylim=c(0,2), col='red')
abline(h=1, lty=2)


# SNR gain in dB
plot(log(S2,2), 20*log10(S2/N2), xlim=c(-6,0), ylim=c(0,30),
     main='DxO PureRAW SNR enhacement',
     xlab='RAW exposure (EV)', ylab='SNR (dB)')
lines(log(S2,2), 20*log10(S1/N1), col='red')

# SNR gain in EV
plot(log2(S2), log2(S2/N2), xlim=c(-6,0), ylim=c(0,5),
     main='DxO PureRAW SNR enhacement',
     xlab='RAW exposure (EV)', ylab='SNR (EV)')
lines(log2(S2), log2(S1/N1), col='red')


# DR gain in dB
plot(log2(S2), 20*log10(S2/N2)-20*log10(S1/N1), xlim=c(-6,0), ylim=c(0,10),
     main='DxO PureRAW DR enhacement',
     xlab='RAW exposure (EV)', ylab='DR gain (dB)', col='red')
abline(h=mean(20*log10(S2/N2)-20*log10(S1/N1)), lty=2)

# DR gain in EV
plot(log2(S2), log2((S2/N2)/(S1/N1)), xlim=c(-6,0), ylim=c(0,1.5),
     main='DxO PureRAW DR enhacement',
     xlab='RAW exposure (EV)', ylab='DR gain (EV)', col='red')
abline(h=mean(log2((S2/N2)/(S1/N1))), lty=2)
