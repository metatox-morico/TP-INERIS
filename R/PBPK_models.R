# first install library deSolve
# load library deSolve
library(deSolve)

# one compartment model
model1 <- function(t, Conc, k) { 
  dC  <- -k * Conc  
  return(list(dC=dC))
}

# times at which the concentration is estimated
times  <- seq(from=0, to=5, by=0.25) #h

# initial state
Conc<-200 # microM/L
# parameter
k<-2 #/h

# simulation
out1 <- ode(times=times, func=model1, y=Conc, parms=k)
# graph
plot(out1, type="b", xlab="time (h)", ylab=expression("Concentration ("*mu*M.L^{-1}*")"))







# same model with a test on parameters, and allowing for several parameters (only k in fact)

model1bis <- function(t, C, parameters) { 
  with (as.list(parameters), {
    if ((C>0)){
      if (is.na(k))
        stop("missing k")
    }
    dC  <- -k * C  
    return(list(dC))
  })
}
parms<-c("k"=2)
out1 <- ode(times=times, func=model1bis, y=C, parms=parms)
plot(out1, type="b", xlab="time (h)", ylab=expression("Concentration ("*mu*M.L^{-1}*")"))






# Two-compartment model with diffusion between compartments, 
# 4 parameters (2 compound-specific, 2 model-specific), 1st order metabolism 
model2 <- function(t, Y, parameters) { 
  with(as.list(Y), {
    with (as.list(parameters), {
      CCentral <- QCentral / VCentral
      CPeriph  <- QPeriph / VPeriph
      
      dQCentral <- kdiff * (CPeriph-CCentral) - CL * CCentral
      dQPeriph  <- -kdiff * (CPeriph-CCentral) 

      return(list(c(dQCentral, dQPeriph), # derivatives
           c(CCentral = CCentral, CPeriph = CPeriph)))    # outputs
    })
  })
}
parms2<-c("kdiff"=1, "CL"=10, "VCentral"=10, "VPeriph"=2)
# initial state
Y2<-c("QCentral"=0, "QPeriph"=100)
# simulation
out2 <- ode(times=times, func=model2, y=Y2, parms=parms2)
# graph
plot(out2, type="b", xlab="time (h)", ylab=expression("Concentration ("*mu*M.L^{-1}*")"), 
     ylim=list(c(0,30), c(0,100), c(0,3), c(0,50)))

#Calulate AUC
sum(out2[,"CCentral"]) # beware: the concentration does not reach 0

# calculate half-life
CMax_CP<-max(out2[,"CPeriph"])
abs_diff<-abs(out2[,"CPeriph"]-CMax_CP/2)
Half_life_CPeriph<-times[which(abs_diff==min(abs_diff))]

CMax_CC<-max(out2[,"CCentral"])
CC_aftermax<-out2[,"CCentral"][-c(1:(which(out2[,"CCentral"]==CMax_CC)-1))]
abs_diff<-abs(CC_aftermax-CMax_CC/2)
Half_life_CCentral<-times[which(abs_diff==min(abs_diff))+which(out2[,"CCentral"]==CMax_CC)-1]




# full PBPK model for acetaminophene. oral absorption, 1st order metabolism, 1st order renal elimination

model <- function(t, Y, parameters) { 
  with(as.list(Y), {
    with (as.list(parameters), {
      if ((Qadmin>0)){
        if (is.na(kGut))
          stop("missing kGut")
        if (is.na(Frac))
          stop("missing Frac")
      }
      
      CFat  <- QFat / (BW*scVFat)
      CPoor <- QPoor / (BW*scVPoor)
      CRich <- QRich / (BW*scVRich)
      CLiver <- QLiver / (BW*scVLiver)
      CSkin <- QSkin / (BW*scVSkin)
      CArt <- QArt / (scVBlood * BW/3) 
      CVen <- QVen / (scVBlood * BW*2/3) 
      
      FFat   <- FBlood * (BW * scFFat)
      FPoor <- FBlood * (BW * scFPoor)
      FRich <- FBlood * (BW * scFRich)
      FLiver <- FBlood * (BW * scFLiver)
      FSkin <- FBlood * (BW * scFSkin)
      
      dQFat  <- FFat * (CArt - CFat / PCFat) 
      dQPoor <- FPoor * (CArt - CPoor / PCPoor) 
      dQRich <- FRich * (CArt - CRich / PCRich)
      dQGut  <- -kGut * QGut  # NA if kGut=NA (no oral).Should be 0 and Qgut=0.  replace kGut by 0 if NA and QGut==NA and Qadmin==NA? If Qadmin=0 OK, if QGut=0 OK.
      dQLiver <- sum(FLiver * (CArt - CLiver / PCLiver), 
                     kGut * QGut, 
                     - CLH * fup * CLiver / PCLiver, na.rm=TRUE)  
      dQSkin <- FSkin  * (CArt - CSkin  / PCSkin)
      
      dQArt <- FBlood * (CVen - CArt)  * BW 
      
      dQVen <-  FFat * CFat / PCFat + 
        FPoor * CPoor / PCPoor +
        FRich * CRich / PCRich +
        FLiver * CLiver / PCLiver +
        FSkin * CSkin / PCSkin - 
        FBlood * CVen * BW -
        Gr * fup * CVen
      
      list(c(dQVen, dQArt, dQFat, dQPoor, dQRich, dQGut, dQLiver, dQSkin), # derivatives
           c(CVen = CVen,     # outputs
             CArt = CArt, 
             CFat = CFat, 
             CPoor = CPoor,
             CRich = CRich, 
             CLiver = CLiver,
             CSkin = CSkin))
    })
  })
}
# parameters
# already set:
parms_physio <- c("scFFat"=0.05, "scFPoor"=0.17,   "scFRich"=0.48, "scFLiver"=0.25, "scFSkin"=0.05,  # no unit
                  "FBlood"=0.089,     # litres / minute / kg
                  "scVFat"=0.25, "scVPoor"=0.46,   "scVRich"=0.08, "scVLiver"=0.025, "scVSkin"=0.11, "scVBlood"=0.075, # no unit
                  "BW"=73, # kg
                  "Gr"=0.035)  # litres / minute recalculated based on estimated renal elimination
# chemical-dependant. Non relevant parameters can be set to NA
parms_chemical <- c("PCFat"=0.25, "PCPoor"=0.66,   "PCRich"=0.77, "PCLiver"=0.77, "PCSkin"=0.66,  # no unit
                    "CLH"=10, # hepatic clearance
                    "fup"=0.88,
                    "Frac"=0.95,
                    "kGut"=0.025)   #/min  
# exposure parameters for single exposure
parms_exposure <- c("Qadmin"= 500*1000*73/151) # oral

parms<-c(parms_physio, parms_chemical, parms_exposure)

# initial state
Y <- c("QVen"=0, "QArt"=0, "QFat"=0, "QPoor"=0, "QRich"=0, "QGut"=parms["Frac"] * parms["Qadmin"], "QLiver"=0, "QSkin"=0 )
names(Y)<-c("QVen", "QArt", "QFat", "QPoor", "QRich", "QGut", "QLiver", "QSkin")    

# times at which the function should be evaluated
times  <- seq(from=0, to=1000, by=10)





# calculations
# single dose
out <- ode(times=times, func=model, y=Y, parms=parms)
# plot:
jpeg("out.jpg", width=20, height=20, res=300, unit="cm", pointsize=12)
par(mar=c(3,3.5, 3, 0.5), mgp=c(1.8, 0.8, 0))
plot(out, which=c("CVen", "CArt", "CFat", "CPoor", "CRich", "QGut", "CLiver", "CSkin"), type="b", xlab="time (min)", ylab=expression("Concentration ("*mu*M.L^{-1}*")"))
dev.off()
