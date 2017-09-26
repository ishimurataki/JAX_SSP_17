plot.mediation_letmesleep <- function(med, 
                                      col="firebrick4", 
                                      chrlen=mouse.chrlen, ...){
  library(RColorBrewer)
  
  names(med) = toupper(names(med))
  
  #Check input
  stopifnot(c("CHR", "POS","LOD") %in% names(med))
  
  if (!("GMB" %in% names(med)))
    med$GMB <- gmb.coordinates(med$CHR, med$POS, chrlen=chrlen)
  
  # reorganize chr-lengths as in gmb.coordinates
  max.chr <- max(as.numeric(med$CHR[grep("[0-9]+", med$CHR)]))
  unique.chr <- levels(factor(factor(med$CHR, levels=c(1:max.chr, "X", "Y", "M"))))
  chrlen = chrlen[unique.chr]
  
  # if pos > MAXPOS, then pos is NOT in Mb and needs conversion
  MAXPOS = 3000
  if (max(chrlen) > MAXPOS) chrlen <- chrlen / 10^6
  
  ### Create X-axis gmb values for mediators (e.g. proteins/RNA) genome positions
  gmb.mediation.extended = c(min(med$GMB)-30,med$GMB,max(med$GMB)+30) #Used to extend the x-axis of the plot a bit
  chrlen = c(0, cumsum(chrlen))
  chrmid = chrlen[-length(chrlen)] + (diff(chrlen) * 0.5)
  
  # creating color distribution spread from the average
  med$strayfromLOD <- log(abs(med$LOD - mean(med$LOD))) - min(log(abs(med$LOD - mean(med$LOD))), na.rm = TRUE)
  breaks <- seq(from = min(med$strayfromLOD), to = max(med$strayfromLOD), length.out = 18)
  for(i in 1:21771){
    print(i)
    if(med$CHR[i] %in% c("1", "3", "5", "7", "9", "11", "13", "15", "17", "19", "Y")){
      med$col[i] <- "#2d7c00"
    }
    if(med$CHR[i] %in% c("2", "4", "6", "8", "10", "12", "14", "16", "18", "X")){
      med$col[i] <- "#1DB409"
    }
  }
  
  # 5E4FA2, 378EBA, 75C8A4, BEE5A0, F1F9A9, FEEDA2, FDBE6E, FDBE6E, F67B49, D8434D, 9E0142
  
  #### Create the plot
  par(font = 2, font.lab = 2, font.axis = 2, xaxs="i",las = 1, mar=c(3, 4, 3, 1) + 0.1)
  plot(gmb.mediation.extended, c(0,med$LOD,0), col = 0, ylim = c(0, max(med$LOD)*1.05),ylab = "Conditioned LOD", xaxt = "n", xlab = "", ...)
  usr = par("usr")
  rect(chrlen[2 * 1:(length(chrlen) * 0.5)-1], usr[3],
       chrlen[2 * 1:(length(chrlen) * 0.5)], usr[4],
       border = "grey60", lty=3, col = rgb(0.96,0.96,0.96))
  rect(usr[1], usr[3], usr[2], usr[4], border = "black")
  points(med$GMB, med$LOD, type = "p", pch=21,lty = 1,lwd=1, cex=0.9, col= med$col) #dark red points
  text(chrmid, 0.97 * usr[4], names(chrlen)[-1], cex=1)
  colnumbers <- which(med$LOD %in% sort(med$LOD)[1:10])
}