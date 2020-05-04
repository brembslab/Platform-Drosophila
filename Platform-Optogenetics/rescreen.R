rm(list=ls()) #clean memory
gc()          #collect garbage

############################################# --- begin user data --- #############################################
datadir="D:/data/optogenetics"   #where is the data located
codedir="B:/GitHub/screen-analysis-Joystick" #location of other R or Rmd files used in this script, normally location of this script
htmlname="rescreen.html" #filename for HTML evaluation sheet
groupfilename="rescreen_all.txt" #filename for text file with datafiles assigned to experimental groups
############################################## --- end user data --- ##############################################

#load libraries
library(gridExtra)
library(reshape2)
library(ggplot2)


#################################################################### Functions ####################################################

# Check if the light is on or of depending of the hysteresis and the previous trace
check_switch <- function(trace_point, switch_off, switch_on) {
  
  if (trace_point == FALSE && switch_on == TRUE) {trace_point <- TRUE}
  
  if (trace_point == TRUE && switch_off == TRUE) {trace_point <- FALSE}
  
  return(trace_point)	# return the last timestamp in the dataframe
}

# make a vector of light on or of for the trace and hysteresis
light_state <- function(trace, Hysteresis) {
  
  switch_on <- (trace > Hysteresis) # potential signals to turn on the light
  switch_off <- (trace < -Hysteresis) # potential switch-off signals
  state <- vector("logical",length= lengthExp) # vector allocation for light ON/OFF state
  
  if (trace[1] > Hysteresis) {state[1]<-TRUE} 
  
  for (i in 2:lengthExp){
    state[i] <- state[i-1]
    state[i] <- check_switch(state[i], switch_off[i], switch_on[i])
  }
  
  return(state)
}

# Calculate a PI from a boolean light state vector
calculate_PI <- function(state_vector) {
  PI <- (sum(state_vector)-sum(!state_vector))/length(state_vector)
  return(PI)
}

all_false <- function(thres,trace,window=4800) {
  
  lag20 <- abs(diff(trace, lag = 20))
  over_thres <- lag20 > thres
  flat <- vector("logical",length(trace)-window)
  for (i in 1:(length(trace)-window))
  {
    flat[i]<-any(over_thres[i:(i+window)])
  }
  
  keep<-all(flat)
  
  return(keep)
  
}

#samplesizes annotations
samplesizes.annotate <- function(boxes, samplesizes)
{
  annotate("text",
           x=boxes,
           y=-Inf, vjust = -0.5,
           label=paste("N=", samplesizes[boxes]))
}


########################################## Initialize some parameters  ###########################################
tested_flies <- read.table(paste(datadir,"/",groupfilename, sep = ""), quote="\"", comment.char="#") #load group file
setwd(datadir)
all_screens <- unique(tested_flies$V2)
no_of_screens <- length(all_screens)
skip<-37
flat_thres <- 0.9
PI_thres <- 1 #flies rejected for pretest PI=1

groupedPIs <- list()

on_wiggle1<- NA
on_wiggle2<- NA
on_wiggle3 <- NA


off_wiggle1<- NA
off_wiggle2<- NA
off_wiggle3 <- NA


Nexp_total <- 60 # there won't be more experiments than these #round(length(tested_flies$V1)/no_of_screens)*4

effectsize_mat <- matrix(NA, 3*Nexp_total, no_of_screens)
just_rein_mat <- matrix(NA, 3*Nexp_total, no_of_screens)
used_traces <- as.numeric(rep(NA,no_of_screens))

on_wiggle_mat <- matrix(NA, 3*Nexp_total, no_of_screens)
off_wiggle_mat <- matrix(NA, 3*Nexp_total, no_of_screens)
all_PIs <- as.numeric(rep(NA,10))

for(j in 1:no_of_screens)
  
{
  count <- 0 # a variable that keeps a count of the number of graphs taken 
  
  Nexp <- length(tested_flies$V1[tested_flies$V2==all_screens[j]])
  data_files <- as.character(tested_flies$V1[tested_flies$V2==all_screens[j]])

  
  
  PI_platform <- matrix(NA, 3*Nexp, 10)    # Variable where PIs are saved
  combo_PI_platform <- matrix(NA, 3*Nexp, 5) # variable where the mean of PI's are stored for 2-2 segments
  
  
  on_wiggle <- as.numeric(rep(NA,3*Nexp))
  off_wiggle <- as.numeric(rep(NA,3*Nexp))
  diff_wiggle <- as.numeric(rep(NA,3*Nexp))
  
  effectsize <- as.numeric(rep(NA,3*Nexp)) # variable defining effectsize score
  effectsize_reinf <- as.numeric(rep(NA,3*Nexp)) # variable defining effectsize score

  group_name <- all_screens[j]
  
  # Start a for loop for the number of files to analyze
  for(i in 1:Nexp){

    
    ############################## Import in a dataframe just the values. ################################
    
     data <- read.table(data_files[i], header = FALSE, sep = "\t", quote = "\"" , dec = ".", fill = TRUE, skip = skip , comment.char = "", nrows = 24037-skip,col.names=c("n","t.s.","pos1","pos2","pos3"))
    
    
    
    ################################### Import in a dataframe the information of the experiments. ############################################
    
    
    info <-read.table(data_files[i], header = FALSE, sep = "", 
                      col.names = paste0("V",seq_len(20)), fill = TRUE)
    info  <- info[1:20,]
    
    
    ######################################## Extracting some parameters from the meta data #########################################
    
    lengthExp <- length(data$t.s.)     # Number of data points
    #Hysteresis <- as.numeric(as.character(info$V3[19]))  # Hysteresis of the experiment
    Hysteresis <- 0
    
    # Side of light for each platform
    light_side1 <- c(as.character(info$V3[6]),as.character(info$V4[6]),as.character(info$V5[6]),as.character(info$V6[6]),as.character(info$V7[6]),as.character(info$V8[6]),as.character(info$V9[6]),as.character(info$V10[6]),as.character(info$V11[6]),as.character(info$V12[6])) 
    light_side2 <- c(as.character(info$V3[10]),as.character(info$V4[10]),as.character(info$V5[10]),as.character(info$V6[10]),as.character(info$V7[10]),as.character(info$V8[10]),as.character(info$V9[10]),as.character(info$V10[10]),as.character(info$V11[10]),as.character(info$V12[10])) 
    light_side3 <- c(as.character(info$V3[14]),as.character(info$V4[14]),as.character(info$V5[14]),as.character(info$V6[14]),as.character(info$V7[14]),as.character(info$V8[14]),as.character(info$V9[14]),as.character(info$V10[14]),as.character(info$V11[14]),as.character(info$V12[14])) 
    
    right_platform1 <- all(light_side1=="right")
    right_platform2 <- all(light_side2=="right")
    right_platform3 <- all(light_side3=="right")
    
    TimeExp <- data$t.s.[lengthExp]   # The total time it took for the experiment to complete
    data$Sampling<-c(0,diff(data$t.s., lag = 1)) # Calculating Inter Sample intervals (ISI)
    MaxSample<-max(data$Sampling)  # Checking what it the maximal ISI
    # Segment rechnen und plotten mit dem trace
    
    segment<- seq(from = 0,to = lengthExp, lengthExp/10)
 
    ############################################################### PI calculation ##############################################
    
    # Save in state variable if the light is on or off taking care of hysteresis
    data$state1 <- light_state(data$pos1,Hysteresis)
    data$state2 <- light_state(data$pos2,Hysteresis)
    data$state3 <- light_state(data$pos3,Hysteresis)
    
    # Change the ON or OFF state if the platform were set to reinforce left
    if(right_platform1==FALSE){ data$state1 <- !data$state1}
    if(right_platform2==FALSE){ data$state2 <- !data$state2}
    if(right_platform3==FALSE){ data$state3 <- !data$state3}
    
    # Condition 1 for the flat line
    keep1<-all_false(flat_thres,data$pos1)
    keep2<-all_false(flat_thres,data$pos2)
    keep3<-all_false(flat_thres,data$pos3)
    
    # Condition 2: Calculate PIs and sort the ones with good pretest
    
    PI_platform1 <- vector("numeric", length = 10)
    #combo_PI_platform1 <- vector("numeric", length = 5)
    PI_platform2 <- vector("numeric", length = 10)
    #combo_PI_platform2 <- vector("numeric", length = 5)
    PI_platform3 <- vector("numeric", length = 10)
    #combo_PI_platform3 <- vector("numeric", length = 5)
    
    
    effectsize1 <- vector("numeric", 1)
    effectsize2 <- vector("numeric", 1)
    effectsize3 <- vector("numeric", 1)
    
    wiggle1<-abs(diff(data$pos1))
    wiggle2<-abs(diff(data$pos2))
    wiggle3<-abs(diff(data$pos3))
    
    on_wiggle1 <- mean(wiggle1[data$state1[-24000]], na.rm = TRUE)
    off_wiggle1 <- mean(wiggle1[!data$state1[-24000]], na.rm = TRUE)
    on_wiggle2 <- mean(wiggle2[data$state2[-24000]], na.rm = TRUE)
    off_wiggle2 <- mean(wiggle2[!data$state2[-24000]], na.rm = TRUE)
    on_wiggle3 <- mean(wiggle3[data$state3[-24000]], na.rm = TRUE)
    off_wiggle3 <- mean(wiggle3[!data$state3[-24000]], na.rm = TRUE)
    
     for(oo in 1:10){
      PI_platform1[oo] <- calculate_PI(data$state1[segment[oo]:segment[oo+1]])
      PI_platform2[oo] <- calculate_PI(data$state2[segment[oo]:segment[oo+1]])
      PI_platform3[oo] <- calculate_PI(data$state3[segment[oo]:segment[oo+1]])
    }
    

    #subtract pretest
    effectsize1 <- mean(PI_platform1[c(2,3,4,5,6,7,8,9)])-PI_platform1[1]
    effectsize2 <- mean(PI_platform2[c(2,3,4,5,6,7,8,9)])-PI_platform2[1]
    effectsize3 <- mean(PI_platform3[c(2,3,4,5,6,7,8,9)])-PI_platform3[1]
    
    just_reinf1 <- mean(PI_platform1[c(2,3,4,5,6,7,8,9)])
    just_reinf2 <- mean(PI_platform2[c(2,3,4,5,6,7,8,9)])
    just_reinf3 <- mean(PI_platform3[c(2,3,4,5,6,7,8,9)])
    
    # Condition 3: Tally light encounters in the first training period
    light_encounter1 <- sum(abs(diff(data$state1[(lengthExp/10):((lengthExp*2)/10)])))
    light_encounter2 <- sum(abs(diff(data$state2[(lengthExp/10):((lengthExp*2)/10)])))
    light_encounter3 <- sum(abs(diff(data$state3[(lengthExp/10):((lengthExp*2)/10)])))
    
    
    
    if(keep1 && abs(PI_platform1[1])<PI_thres & light_encounter1>1){
      count<-count+1
      #pos <- pos + data$pos1
      PI_platform[((3*i)-2),] <- PI_platform1
      #combo_PI_platform[((3*i)-2),] <- combo_PI_platform1
      effectsize[((3*i)-2)] <- effectsize1
      effectsize_reinf[((3*i)-2)] <- just_reinf1
      on_wiggle[((3*i)-2)] <- on_wiggle1
      off_wiggle[((3*i)-2)] <- off_wiggle1
    }
    
    if(keep2 && abs(PI_platform2[1])<PI_thres & light_encounter2>1){
      count<-count+1
      #pos <- pos + data$pos2
      PI_platform[((3*i)-1),] <- PI_platform2
      #combo_PI_platform[((3*i)-1),] <- combo_PI_platform2 
      effectsize[((3*i)-1)] <- effectsize2
      effectsize_reinf[((3*i)-1)] <- just_reinf2
      on_wiggle[((3*i)-1)] <- on_wiggle2
      off_wiggle[((3*i)-1)] <- off_wiggle2
    }
    
    if(keep3 && abs(PI_platform3[1])<PI_thres & light_encounter3>1){
      count<-count+1
      #pos <- pos + data$pos3
      PI_platform[((3*i)),] <- PI_platform3
      #combo_PI_platform[(3*i),] <- combo_PI_platform3
      effectsize[(3*i)] <- effectsize3
      effectsize_reinf[((3*i))] <- just_reinf3
      on_wiggle[((3*i))] <- on_wiggle3
      off_wiggle[((3*i))] <- off_wiggle3
    }
  }
  
  effectsize_mat[1:length(effectsize),j] <- effectsize
  on_wiggle_mat[1:length(on_wiggle),j] <- on_wiggle
  off_wiggle_mat[1:length(off_wiggle),j] <- off_wiggle
  just_rein_mat[1:length(effectsize_reinf),j] <- effectsize_reinf
  all_PIs <- rbind(all_PIs,PI_platform)
  used_traces[j]<-count

  # Boxplot of the PIs
  boxplot(PI_platform, col="grey",xlab="",ylab="PI",main=group_name, ylim = c(-1, 1),names=c("Pretest","Training","Training","Training","Training","Training","Training","Training","Training","Test"), cex.lab=1.5, cex.axis = 1.2)
  abline(h = 0, untf = FALSE, col="black",lwd=3)
  group <- NULL
  for(s in 1:10){
    a <- rep(s,Nexp*3)
    group <- append(group,a)
  }
  stripchart(as.vector(PI_platform)~group,vertical = TRUE, method = "jitter",pch = 21, col = "maroon", bg = "bisque",add = TRUE) 
  
  accepted_flies <- which(!is.na(PI_platform[,1]))
  
  matplot(t(PI_platform[accepted_flies,]), type = c("b"),xlab=group_name,main=group_name, pch=1,col = 1:4)

colnames(PI_platform)<-c("pretest","training1","training2","training3","training4","training5","training6","training7","training8","test")
groupedPIs[[j]]=as.data.frame(na.omit(PI_platform)) #store the PIs without the excluded flies
    
}  #for nofScreens

names(groupedPIs) <- as.character(all_screens) #name the list elements in the list if PIs

colnames(effectsize_mat) <- as.character(all_screens)
colnames(just_rein_mat) <- as.character(all_screens)

boxplot(effectsize_mat,col="yellow",xlab="",ylab="reinforcement - pretest", names = as.character(all_screens), cex.lab=1.0, cex.axis = 1, las=2)
abline(h = 0, untf = FALSE, col="black",lwd=3)

boxplot(just_rein_mat,col="yellow",xlab="",ylab="just_reinforcement", names = as.character(all_screens), cex.lab=1.0, cex.axis = 1, las=2)
abline(h = 0, untf = FALSE, col="black",lwd=3)

diff_wiggle_mat=on_wiggle_mat-off_wiggle_mat
boxplot(diff_wiggle_mat,col="yellow",xlab="",ylab="wiggle difference(on - off)", names = as.character(all_screens),cex.lab=1.0, cex.axis = 1, las = 2)
abline(h = 0, untf = FALSE, col="black",lwd=3)

## Reinforcement barplot

mean_effectsize <- apply(effectsize_mat,2,function(effectsize_mat){mean(effectsize_mat, na.rm = TRUE)})
std_dev_effectsize <- apply(effectsize_mat,2,function(effectsize_mat){sd(effectsize_mat, na.rm = TRUE)})

std_err_effectsize <- std_dev_effectsize/count**0.5

barCenters <- barplot(height = mean_effectsize,
                      beside = true, las = 2,ylim=c(-1,1),
                      cex.names = 1,
                      main = "mixed Reinforcement with yellow light",
                      ylab = "Reinforcement",
                      border = "black", axes = TRUE)


segments(barCenters, mean_effectsize - std_err_effectsize, barCenters,
         mean_effectsize + std_err_effectsize, lwd = 1.5)

## Wiggle barplot

mean_diff_wiggle <- apply(diff_wiggle_mat,2,function(diff_wiggle_mat){mean(diff_wiggle_mat, na.rm = TRUE)})
std_dev_diff_wiggle <- apply(diff_wiggle_mat,2,function(diff_wiggle_mat){sd(diff_wiggle_mat, na.rm = TRUE)})
std_err_diff_wiggle <- std_dev_diff_wiggle/count**0.5

barCenters <- barplot(height = mean_diff_wiggle,
                      beside = true, las = 2,ylim=c(-0.3,0.3),
                      cex.names = 1,
                      main = "mixed Wiggle with yellow light",
                      ylab = "wiggle difference(on - off)",
                      border = "black", axes = TRUE)


segments(barCenters, mean_diff_wiggle - std_err_diff_wiggle, barCenters,
         mean_diff_wiggle + std_err_diff_wiggle, lwd = 1.5)


## Total PI overall and number of experiments per line

boxplot(all_PIs,main="overall PI for the whole screen - mixed")

names(used_traces)<-all_screens
barplot(used_traces,las=2,main = "mixed experiments")


####################################

#### call RMarkdown to generate HTML file with evaluations #################################
rmarkdown::render(paste(codedir,"/rescreen.Rmd", sep=""), 
                  output_file = paste(htmlname), 
                  output_dir = datadir)
#### end RMarkdown for project evaluations #################################################
