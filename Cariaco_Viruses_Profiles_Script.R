# For plotting biogeochemical profiles from CAR216 for Mara et al virus paper
# Started 3/20/19

# rm(list = ls())

#install.packages("oce")
#install.packages("gridGraphics")
library(oce)
library(dplyr)
library(tidyverse)
library(readxl)
library(cowplot)
library(gridGraphics)

# First try cast 2, from which 3 of the virus samples came
C216_2_2 =read.oce("C216_2/C216_2/C216_2_2.cnv")
plot(C216_2_2)
plotProfile(C216_2_2, ytype = "depth", "oxygen5") # Note that both the down and upcasts are plotted

### Subset to remove downcast and equilibration data
# To do this, plot indices of consecutive datapoints against pressure to determine where downcast ends and upcast begins
plotScan(C216_2_2)

#locator() 
# this is interactive. Click on plot at beginning and end of upcast. Then hit Escape, and it will tell you 
# values where you clicked. Then you know where to trim the dataset by scan number.

# Use plotScan to check if ctdTrim is cutting correct portions
plotScan(ctdTrim(C216_2_2, "range", parameters = list(item = "scan", from= 930, to = 1810)))
C216_2_2_trim <- ctdTrim(C216_2_2, "range", parameters = list(item = "scan", from= 930, to = 1810))
plotProfile(C216_2_2_trim, ytype = "depth", "oxygen5")

# import cast 4 (where other 3 samples were taken) and do same
C216_2_4 =read.oce("C216_2/C216_2/C216_2_4.cnv")
plot(C216_2_4)
plotProfile(C216_2_4, ytype = "depth", "oxygen5")
plotScan(C216_2_4)
#locator() 
plotScan(ctdTrim(C216_2_4, "range", parameters = list(item = "scan", from= 915, to = 1810)))
C216_2_4_trim <- ctdTrim(C216_2_4, "range", parameters = list(item = "scan", from= 930, to = 1810))
plotProfile(C216_2_4_trim, ytype = "depth", "oxygen5")

# Plot both oxygen data from both casts on same plot
plotProfile(C216_2_2_trim, ytype = "depth", "oxygen5")
lines(C216_2_4_trim[["oxygen5"]], C216_2_4_trim[["depth"]], lty = "dashed")


# Worked but both of these profiles need to be corrected by the sensor offset
# Based on O2 values in anoxic zone, I previously determined this offset to be an average of 2.30 for cast 2 and 2.29umol/kg for cast 4
C216_2_2_trim[["oxygen5"]] <- C216_2_2_trim[["oxygen5"]]-2.3
C216_2_4_trim[["oxygen5"]] <- C216_2_4_trim[["oxygen5"]]-2.29

# Plot corrected oxygen data from both casts on same plot again
plotProfile(C216_2_2_trim, ytype = "depth", "oxygen5")
lines(C216_2_4_trim[["oxygen5"]], C216_2_4_trim[["depth"]], lty = "dashed")





# This is fine but using this package, you can't plot 2 variables in same profile (except temp and salinity) so I am going to pull 
# data from the CTD dataset item and plot using ggplot or something
C216_2_2 <- tibble("depth" = C216_2_2_trim[["depth"]], "oxygen" = C216_2_2_trim[["oxygen5"]], "BAT" = C216_2_2_trim[["beamAttenuation"]])
C216_2_4 <- tibble("depth" = C216_2_4_trim[["depth"]], "oxygen" = C216_2_4_trim[["oxygen5"]], "BAT" = C216_2_4_trim[["beamAttenuation"]])

# Also import microbial count and geochemical data from microscopy from C216_2
# NOTE that bacterial counts were made during C216_2 and C216_3 while VLP and nutrients, sulfide were measured during C216_3.
# For now just plot C216_2 data
x <- read_excel("C216_2/C216_2_MicroscopyData.xlsx")
#colnames(x) <- NULL
C216_2_VLP <- tibble("depth" = x$`Depth (m)`[1:6], "VLPx10_8_L-1" = x$`x10^8 VLP L-1`[1:6], "VLP_SD" = x$`SD x10^8 L-1..4`[1:6])
C216_2_prok <- tibble("depth" = x$Depth, "Proksx10_8_L-1" = x$`x10^8 Proks L-1`, "Prok_SD" = x$`SD x10^8 L-1..8`)


# Combine VLP and Prok abundance data into same tibble as C216_2_2 so that can plot them in same figure
temp <- full_join(C216_2_2, C216_2_VLP)
C216_2_2 <- full_join(temp,C216_2_prok)



# Plot and save 
svg(filename = "Figures/C216_2_oxygen_microabund.svg")
par(mar=c(5,6,5,11)+.1, xaxs="i", yaxs="i") # needed to adjust margins to fit plot
with(C216_2_2, plot(C216_2_2$oxygen, C216_2_2$depth, "l", lwd = 3, ylab = "depth [m]", xlab = "", ylim = rev(c(0,1000)), axes=FALSE)) # Do not plot any axes
lines(C216_2_4$oxygen,C216_2_4$depth,col = "grey", lwd = 3)
axis(3)   # Draw the x-axis above the plot area
axis(2)   # Draw the y-axis to the left of the plot area
mtext("oxygen [µM]", side=3, line=3) # O2 label on top
box() # puts box back around plot
# Add in VLP 
par(new = T)
with(C216_2_2, plot(C216_2_2$`VLPx10_8_L-1`, C216_2_2$depth, pch=16, axes=F, xlab=NA, ylab=NA, cex=1.2, ylim = rev(c(0,1000)), xlim = c(0,12)))
lines(C216_2_2$`VLPx10_8_L-1`,C216_2_2$depth,lty="dashed")
# And plot the standard deviation
# This is tricky because for depth profiles, they are x-axis error bars
# Do this by drawing arrows between the +/- SD and with modified "arrowheads"
arrows(C216_2_2$`VLPx10_8_L-1`-C216_2_2$VLP_SD, C216_2_2$depth, C216_2_2$`VLPx10_8_L-1`+C216_2_2$VLP_SD, C216_2_2$depth, length=0.05, angle=0, code=3)
# and Prokaryote abundance
par(new = T)
with(C216_2_2, plot(C216_2_2$`Proksx10_8_L-1`, C216_2_2$depth, pch=17, axes=F, xlab=NA, ylab=NA, cex=1.2, ylim = rev(c(0,1000)), xlim = c(0,12)))
lines((C216_2_2[882:891,]%>%arrange(depth))$`Proksx10_8_L-1`,(C216_2_2[882:891,]%>%arrange(depth))$depth,lty="dotted") # have to use "arrange" because depths were not in descending order in tibble and lines were plotted in zig zaggy way
arrows(C216_2_2$`Proksx10_8_L-1`-C216_2_2$Prok_SD, C216_2_2$depth, C216_2_2$`Proksx10_8_L-1`+C216_2_2$Prok_SD, C216_2_2$depth, length=0.05, angle=0, code=3)
axis(side = 1)
mtext(side = 1, line = 3, parse(text=paste("x10","^8*","L", "^-1", sep=" ", collapse = NULL)))
# Legend
legend(5, 800, legend=c("Oxygen, Cast 2", "Oxygen, Cast 4", "VLP Abundance", "Prokaryote Abundance"), col=c("black", "grey", "black", "black"), lty=c(1,1,2,3), pch = c(NA_integer_, NA_integer_, 16,17), cex=0.8, lwd = c(3,3,1,1))
dev.off()



## Next make a figure from C216_3 with more biogeochemical data (O2, BAT, H2S, nutrients, micro abundance)
# import cast 3 (full depth profile)
C216_3_3 =read.oce("C216_3/C216_3_transmissometerdatanotgood/C216_3_3.cnv")
plot(C216_3_3)
plotProfile(C216_3_3, ytype = "depth", "oxygen5")
plotScan(C216_3_3)
#locator() 
plotScan(ctdTrim(C216_3_3, "range", parameters = list(item = "scan", from= 1325, to = 2625)))
C216_3_3_trim <- ctdTrim(C216_3_3, "range", parameters = list(item = "scan", from= 1325, to = 2625))
plotProfile(C216_3_3_trim, ytype = "depth", "oxygen5")

# Worked but both of these profiles need to be corrected by the sensor offset
# Based on O2 values in anoxic zone, I previously determined the oxygen offset to be an average of 2.36
C216_3_3_trim[["oxygen5"]] <- C216_3_3_trim[["oxygen5"]]-2.36

# Plot corrected oxygen data
plotProfile(C216_3_3_trim, ytype = "depth", "oxygen5")




# Pull data from the CTD dataset item and plot in base R
C216_3_3 <- tibble("depth" = C216_3_3_trim[["depth"]], "oxygen" = C216_3_3_trim[["oxygen5"]], "BAT" = C216_3_3_trim[["beamAttenuation"]])


# Import microbial count and geochemical data from C216_3
x <- read_excel("C216_3/C216_3_metadata.xlsx")

# Combine metadata and CTD data in same tibble
C216_3 <- full_join(C216_3_3,x)


# Plot and save 

# First plot- oxygen, sulfide, micro abundance
svg(filename = "Figures/C216_3_oxygen_h2s_microabund.svg")
par(mar=c(5,6,5,11)+.1, xaxs="i", yaxs="i") # needed to adjust margins to fit plot
with(C216_3, plot(C216_3$oxygen, C216_3$depth, "l", lwd = 3, ylab = "depth [m]", xlab = "", ylim = rev(c(0,1000)), axes=FALSE, xlim = c(0,150))) # Do not plot any axes
lines(C216_3$H2S_µM,C216_3$depth,col = "grey", lwd = 3)
points(C216_3$H2S_µM,C216_3$depth, pch = 16, col = "grey", cex=0.8) # also add points for H2S because they are grab samples (not CTD data)
arrows(C216_3$H2S_µM-C216_3$H2S_SD, C216_3$depth, C216_3$H2S_µM+C216_3$H2S_SD, C216_3$depth, length=0.05, angle=0, code=3, col = "grey")
axis(3)   # Draw the x-axis above the plot area
axis(2)   # Draw the y-axis to the left of the plot area
mtext("oxygen or sulfide [µM]", side=3, line=3) # O2 label on top
box() # puts box back around plot
# ADD Prokaryote abundance
par(new = T)
with(C216_3, plot(C216_3$`Proksx10_8_L_-1`, C216_3$depth, pch=17, axes=F, xlab=NA, ylab=NA, cex=1.2, ylim = rev(c(0,1000)), xlim = c(0,6)))
lines((C216_3[1302:1319,]%>%arrange(depth))$`Proksx10_8_L_-1`,(C216_3[1302:1319,]%>%arrange(depth))$depth,lty="dotted") 
arrows(C216_3$`Proksx10_8_L_-1`-C216_3$Prok_SD, C216_3$depth, C216_3$`Proksx10_8_L_-1`+C216_3$Prok_SD, C216_3$depth, length=0.05, angle=0, code=3)
axis(side = 1)
mtext(side = 1, line = 3, parse(text=paste("x10","^8*","L", "^-1", sep=" ", collapse = NULL)))
# Legend
legend(2.3, 800, legend=c("CTD Oxygen", "Sulfide", "Prokaryote Abundance"), col=c("black", "grey", "black"), lty=c(1,1, 3), pch = c(NA_integer_, 16, 17), cex=0.8, lwd = c(3,3, 1))
dev.off()

# Second plot- oxygen and sulfide with NH4 and PO4
svg(filename = "Figures/C216_3_oxygen_h2s_nh4_po4.svg")
par(mar=c(5,6,5,11)+.1, xaxs="i", yaxs="i") # needed to adjust margins to fit plot
with(C216_3, plot(C216_3$oxygen, C216_3$depth, "l", lwd = 3, ylab = "depth [m]", xlab = "", ylim = rev(c(0,1000)), axes=FALSE, xlim = c(0,150))) # Do not plot any axes
lines(C216_3$H2S_µM,C216_3$depth,col = "grey", lwd = 3)
points(C216_3$H2S_µM,C216_3$depth, pch = 16, col = "grey", cex=0.8) # also add points for H2S because they are grab samples (not CTD data)
arrows(C216_3$H2S_µM-C216_3$H2S_SD, C216_3$depth, C216_3$H2S_µM+C216_3$H2S_SD, C216_3$depth, length=0.05, angle=0, code=3, col = "grey")
axis(3)   # Draw the x-axis above the plot area
axis(2)   # Draw the y-axis to the left of the plot area
mtext("oxygen or sulfide [µM]", side=3, line=3) # O2 label on top
box() # puts box back around plot
# Add Nutrients
# NH4
par(new = T)
with(C216_3, plot(C216_3$NH4_µM, C216_3$depth, pch=15, axes=F, xlab=NA, ylab=NA, cex=1.2, ylim = rev(c(0,1000)), xlim = c(0,25)))
lines((C216_3[1302:1319,]%>%arrange(depth))$NH4_µM,(C216_3[1302:1319,]%>%arrange(depth))$depth,lty="dotdash") 
axis(side = 1)
# PO4
par(new = T)
with(C216_3, plot(C216_3$PO4_µM, C216_3$depth, pch=16, axes=F, xlab=NA, ylab=NA, cex=1.2, ylim = rev(c(0,1000)), xlim = c(0,25)))
lines((C216_3[1302:1319,]%>%arrange(depth))$PO4_µM,(C216_3[1302:1319,]%>%arrange(depth))$depth,lty="dashed") 
axis(side = 1)
mtext("ammonia or phosphate [µM]", side = 1, line = 3)
# Legend
legend(7, 775, legend=c("CTD Oxygen", "Sulfide", "Ammonia", "Phosphate"), col=c("black", "grey", "black", "black"), lty=c(1,1,3, 3), pch = c(NA_integer_, 16,15, 16), cex=0.8, lwd = c(3,3,1,1))
dev.off()


# Third plot- oxygen and sulfide with NO3
svg(filename = "Figures/C216_3_oxygen_h2s_no3.svg")
par(mar=c(5,6,5,11)+.1, xaxs="i", yaxs="i") # needed to adjust margins to fit plot
with(C216_3, plot(C216_3$oxygen, C216_3$depth, "l", lwd = 3, ylab = "depth [m]", xlab = "", ylim = rev(c(0,1000)), axes=FALSE, xlim = c(0,150))) # Do not plot any axes
lines(C216_3$H2S_µM,C216_3$depth,col = "grey", lwd = 3)
points(C216_3$H2S_µM,C216_3$depth, pch = 16, col = "grey", cex=0.8) # also add points for H2S because they are grab samples (not CTD data)
arrows(C216_3$H2S_µM-C216_3$H2S_SD, C216_3$depth, C216_3$H2S_µM+C216_3$H2S_SD, C216_3$depth, length=0.05, angle=0, code=3, col = "grey")
axis(3)   # Draw the x-axis above the plot area
axis(2)   # Draw the y-axis to the left of the plot area
mtext("oxygen or sulfide [µM]", side=3, line=3) # O2 label on top
box() # puts box back around plot
# Add Nutrients
# NO3
par(new = T)
with(C216_3, plot(C216_3$NO3_µM, C216_3$depth, pch=15, axes=F, xlab=NA, ylab=NA, cex=1.2, ylim = rev(c(0,1000)), xlim = c(0,13)))
lines((C216_3[1302:1319,]%>%arrange(depth))$NO3_µM,(C216_3[1302:1319,]%>%arrange(depth))$depth,lty="dotdash") 
axis(side = 1)
mtext("nitrate [µM]", side = 1, line = 3)
# Legend
legend(1, 775, legend=c("CTD Oxygen", "Sulfide", "Nitrate"), col=c("black", "grey", "black"), lty=c(1,1,3), pch = c(NA_integer_, 16,15), cex=0.8, lwd = c(3,3,1))
dev.off()


# Fourth plot- oxygen and sulfide with NO2
svg(filename = "Figures/C216_3_oxygen_h2s_no2.svg")
par(mar=c(5,6,5,11)+.1, xaxs="i", yaxs="i") # needed to adjust margins to fit plot
with(C216_3, plot(C216_3$oxygen, C216_3$depth, "l", lwd = 3, ylab = "depth [m]", xlab = "", ylim = rev(c(0,1000)), axes=FALSE, xlim = c(0,150))) # Do not plot any axes
lines(C216_3$H2S_µM,C216_3$depth,col = "grey", lwd = 3)
points(C216_3$H2S_µM,C216_3$depth, pch = 16, col = "grey", cex=0.8) # also add points for H2S because they are grab samples (not CTD data)
arrows(C216_3$H2S_µM-C216_3$H2S_SD, C216_3$depth, C216_3$H2S_µM+C216_3$H2S_SD, C216_3$depth, length=0.05, angle=0, code=3, col = "grey")
axis(3)   # Draw the x-axis above the plot area
axis(2)   # Draw the y-axis to the left of the plot area
mtext("oxygen or sulfide [µM]", side=3, line=3) # O2 label on top
box() # puts box back around plot
# Add Nutrients
# NO2
par(new = T)
with(C216_3, plot(C216_3$NO2_µM, C216_3$depth, pch=16, axes=F, xlab=NA, ylab=NA, cex=1.2, ylim = rev(c(0,1000)), xlim = c(0,0.5)))
lines((C216_3[1302:1319,]%>%arrange(depth))$NO2_µM,(C216_3[1302:1319,]%>%arrange(depth))$depth,lty="dashed") 
axis(side = 1)
mtext("nitrite [µM]", side = 1, line = 3)
# Legend
legend(.05, 800, legend=c("CTD Oxygen", "Sulfide", "Nitrite"), col=c("black", "grey", "black"), lty=c(1,1,2), pch = c(NA_integer_, 16,16), cex=0.8, lwd = c(3,3,1))
dev.off()



### 3/28/19
# After talk with Vivian and Ginny, they just want the plot from C216_2 and C216_3 plot with oxygen sulfide, NO3, and NH4
svg(filename = "Figures/Hydrography_Fig_a.svg")
par(mar=c(5,6,5,11)+.1, xaxs="i", yaxs="i") # needed to adjust margins to fit plot
plot(C216_2_2$oxygen, C216_2_2$depth, "l", lwd = 3, ylab = "depth [m]", xlab = "", ylim = rev(c(0,1000)), axes=FALSE) # Do not plot any axes
  lines(C216_2_4$oxygen,C216_2_4$depth,col = "grey", lwd = 3)
  axis(3)   # Draw the x-axis above the plot area
  axis(2)   # Draw the y-axis to the left of the plot area
  mtext("oxygen [µM]", side=3, line=3) # O2 label on top
  box() # puts box back around plot
  # Add in VLP 
  par(new = T)
  with(C216_2_2, plot(C216_2_2$`VLPx10_8_L-1`, C216_2_2$depth, pch=16, axes=F, xlab=NA, ylab=NA, cex=1.5, ylim = rev(c(0,1000)), xlim = c(0,12)))
  lines(C216_2_2$`VLPx10_8_L-1`,C216_2_2$depth,lty="dashed")
  # And plot the standard deviation
  # This is tricky because for depth profiles, they are x-axis error bars
  # Do this by drawing arrows between the +/- SD and with modified "arrowheads"
  arrows(C216_2_2$`VLPx10_8_L-1`-C216_2_2$VLP_SD, C216_2_2$depth, C216_2_2$`VLPx10_8_L-1`+C216_2_2$VLP_SD, C216_2_2$depth, length=0.05, angle=0, code=3)
  # and Prokaryote abundance
  par(new = T)
  with(C216_2_2, plot(C216_2_2$`Proksx10_8_L-1`, C216_2_2$depth, pch=17, axes=F, xlab=NA, ylab=NA, cex=1.5, ylim = rev(c(0,1000)), xlim = c(0,12)))
  lines((C216_2_2[882:891,]%>%arrange(depth))$`Proksx10_8_L-1`,(C216_2_2[882:891,]%>%arrange(depth))$depth,lty="dotted") # have to use "arrange" because depths were not in descending order in tibble and lines were plotted in zig zaggy way
  arrows(C216_2_2$`Proksx10_8_L-1`-C216_2_2$Prok_SD, C216_2_2$depth, C216_2_2$`Proksx10_8_L-1`+C216_2_2$Prok_SD, C216_2_2$depth, length=0.05, angle=0, code=3)
  axis(side = 1)
  mtext(side = 1, line = 3, parse(text=paste("x10","^8*","L", "^-1", sep=" ", collapse = NULL)))
  # Legend
  legend(3.8, 725, legend=c("CTD Oxygen, Cast 2", "CTD Oxygen, Cast 4", "VLP Abundance", "Prokaryote Abundance"), col=c("black", "grey", "black", "black"), lty=c(1,1,2,3), pch = c(NA_integer_, NA_integer_, 16,17), cex=0.9, lwd = c(3,3,1,1))
dev.off()

# C216_3- oxygen and sulfide with NH4 and NO3- leave out depth axis
svg(filename = "Figures/Hydrography_Fig_b.svg")
par(mar=c(5,6,5,11)+.1, xaxs="i", yaxs="i") # needed to adjust margins to fit plot
plot(C216_3$oxygen, C216_3$depth, "l", lwd = 3, ylab = "", xlab = "", ylim = rev(c(0,1000)), axes=FALSE, xlim = c(0,150))
  lines(C216_3$H2S_µM,C216_3$depth,col = "grey", lwd = 3)
  points(C216_3$H2S_µM,C216_3$depth, pch = 16, col = "grey", cex=1.5) # also add points for H2S because they are grab samples (not CTD data)
  arrows(C216_3$H2S_µM-C216_3$H2S_SD, C216_3$depth, C216_3$H2S_µM+C216_3$H2S_SD, C216_3$depth, length=0.05, angle=0, code=3, col = "grey")
  axis(3)   # Draw the x-axis above the plot area
  axis(2, labels = FALSE)   # Draw the y-axis to the left of the plot area
  mtext("oxygen or sulfide [µM]", side=3, line=3) # O2 label on top
  box() # puts box back around plot
  # Add Nutrients
  # NH4
  par(new = T)
  with(C216_3, plot(C216_3$NH4_µM, C216_3$depth, pch=15, axes=F, xlab=NA, ylab=NA, cex=1.2, ylim = rev(c(0,1000)), xlim = c(0,25)), cex=1.5)
  lines((C216_3[1302:1319,]%>%arrange(depth))$NH4_µM,(C216_3[1302:1319,]%>%arrange(depth))$depth,lty="dotted") 
  axis(side = 1)
  # NO3
  par(new = T)
  with(C216_3, plot(C216_3$NO3_µM, C216_3$depth, pch=16, axes=F, xlab=NA, ylab=NA, cex=1.2, ylim = rev(c(0,1000)), xlim = c(0,25)), cex=1.5)
  lines((C216_3[1302:1319,]%>%arrange(depth))$NO3_µM,(C216_3[1302:1319,]%>%arrange(depth))$depth,lty="dashed") 
  axis(side = 1)
  mtext("ammonia or nitrate [µM]", side = 1, line = 3)
  # Legend
  legend(2, 725, legend=c("CTD Oxygen", "Sulfide", "Ammonia", "Nitrate"), col=c("black", "grey", "black", "black"), lty=c(1,1,3, 2), pch = c(NA_integer_, 16,15, 16), cex=0.9, lwd = c(3,3,1,1))
dev.off()



# Put together as left and right panels in Inkscape

