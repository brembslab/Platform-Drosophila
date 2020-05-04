<<<<<<< Updated upstream:saurabh_modifiedChristian.R


#################################################################### Functions ####################################################

# Check if the light is on or of depending of the hysteresis and the previous trace
check_switch <- function(trace_point, switch_off, switch_on) {
  
  if (trace_point == FALSE && switch_on == TRUE) { trace_point <- TRUE}
  
  if (trace_point == TRUE && switch_off == TRUE) {trace_point <- FALSE}
  
  return(trace_point)	# return the last timestamp in the dataframe
}

# make a vector of light on or of for the trace and hysteresis
ligth_state <- function(trace, Hysteresis) {
  
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

########################################## Initialize some parameters  ###########################################

## Set initial parameters
data_path<-"C:/Users/chise/Desktop/mixed_all/"
flat_thres <- 0.8
PI_thres <- 2 #no flies rejected on this basis
#Hysteresis <- as.numeric(as.character(info$V3[19]))  # Hysteresis of the experiment
Hysteresis <- 0
filename_data<-"/mixed_all.txt"
# setwd("C:/Users/LocalAdmin/Desktop/new_screens/mixed/mixed_all")


#tested_flies <- read.table("C:/Users/LocalAdmin/Desktop/new_screens/all/all.txt", quote="\"", comment.char="")
#tested_flies <- read.table("C:/Users/LocalAdmin/Desktop/new_screens/AMANDA/amanda_all/amanda_all.txt", quote="\"", comment.char="")
#tested_flies <- read.table("C:/Users/LocalAdmin/Desktop/new_screens/SAURABH/saurabh_all/saurabh_all.txt", quote="\"", comment.char="")
tested_flies <- read.table(paste0(data_path,filename_data), quote="\"", comment.char="")
all_screens <- unique(tested_flies$V2)
no_of_screens <- length(all_screens)
skip<-37



########################################################### Allocate variables  ###################################################################
on_wiggle1<- NA
on_wiggle2<- NA
on_wiggle3 <- NA

off_wiggle1<- NA
off_wiggle2<- NA
off_wiggle3 <- NA

#message("please enter the number of experiments")
#Nexp <- scan(n=1,what= numeric())   #Number of flies
Nexp_total <- 40 # they wonÂ´t be more experiment than these #round(length(tested_flies$V1)/no_of_screens)*4
#effectsize <- as.numeric(rep(NA,3*Nexp_total)) # Saurabh variable defining effectsize score
effectsize_mat <- matrix(NA, 3*Nexp_total, no_of_screens)
just_rein_mat <- matrix(NA, 3*Nexp_total, no_of_screens)
used_traces <- as.numeric(rep(NA,no_of_screens))

#on_wiggle <- as.numeric(rep(NA,3*Nexp_total))
#off_wiggle <- as.numeric(rep(NA,3*Nexp_total))
on_wiggle_mat <- matrix(NA, 3*Nexp_total, no_of_screens)
off_wiggle_mat <- matrix(NA, 3*Nexp_total, no_of_screens)
all_PIs <- as.numeric(rep(NA,10))

for(j in 1:no_of_screens)
  
{
  count <- 0 #Defining a variable that keeps a count of the number of graphs taken 
  
  ## Groupname
  group_name <- all_screens[j]
  
  ## Number of experiments in the group. Each experiments runs three platforms at once
  Nexp <- length(tested_flies$V1[tested_flies$V2==group_name])
  
  ## Filenames from group
  data_files <- as.character(tested_flies$V1[tested_flies$V2==group_name])
  #effectsize <- as.numeric(rep(NA,3*Nexp))
  #on_wiggle <- as.numeric(rep(NA,3*Nexp))
  #off_wiggle <- as.numeric(rep(NA,3*Nexp))

  PI_platform <- matrix(NA, 3*Nexp, 10)    # Variable where PIs are saved for each of the ten segments
  combo_PI_platform <- matrix(NA, 3*Nexp, 5) # Variable where the mean of PI's are stored each double segments. There is always two equal consecutive segments
  
  ## Allocating vectors to save data
  on_wiggle <- as.numeric(rep(NA,3*Nexp))
  off_wiggle <- as.numeric(rep(NA,3*Nexp))
  diff_wiggle <- as.numeric(rep(NA,3*Nexp))
  
  effectsize <- as.numeric(rep(NA,3*Nexp)) # Saurabh variable defining effectsize score
  effectsize_reinf <- as.numeric(rep(NA,3*Nexp)) # Saurabh variable defining effectsize score

  # Start a for loop for the number of flies to analyze
  for(i in 1:Nexp){
    #group_name <- tested_flies$V2[i]   #Name of the group
    
    
    ############################## Import in a dataframe just the values. ################################
    
    #data <- read.table(file.choose(), header = FALSE, sep = "\t", quote = "\"" , dec = ".", fill = TRUE, skip = skip , comment.char = "", nrows = 24037-skip,col.names=c("n","t.s.","pos1","pos2","pos3"))
    
    
    #data <- read.table(as.character(tested_flies$V1[x+i]), header = FALSE, sep = "\t", quote = "\"" , dec = ".", fill = TRUE, skip = skip , comment.char = "", nrows = 24037-skip,col.names=c("n","t.s.","pos1","pos2","pos3"))
    #group_name <- tested_flies$V2[x+i]   #Name of the group
    data <- read.table(paste0(data_path,data_files[i]), header = FALSE, sep = "\t", quote = "\"" , dec = ".", fill = TRUE, skip = skip , comment.char = "", nrows = 24037-skip,col.names=c("n","t.s.","pos1","pos2","pos3"))
    
    
    
    ################################### Import in a dataframe the information of the experiments. ############################################
    
    
    info <-read.table(paste0(data_path,data_files[i]), header = FALSE, sep = "", 
                      col.names = paste0("V",seq_len(20)), fill = TRUE)
    info  <- info[1:20,]
    
    
    ######################################## Extracting some parameters from the meta data #########################################
    
    lengthExp <- length(data$t.s.)     # Number of data points
    # Side of the platform for turning the optogenetic light on, for each platform
    light_side1 <- c(as.character(info$V3[6]),as.character(info$V4[6]),as.character(info$V5[6]),as.character(info$V6[6]),as.character(info$V7[6]),as.character(info$V8[6]),as.character(info$V9[6]),as.character(info$V10[6]),as.character(info$V11[6]),as.character(info$V12[6])) 
    light_side2 <- c(as.character(info$V3[10]),as.character(info$V4[10]),as.character(info$V5[10]),as.character(info$V6[10]),as.character(info$V7[10]),as.character(info$V8[10]),as.character(info$V9[10]),as.character(info$V10[10]),as.character(info$V11[10]),as.character(info$V12[10])) 
    light_side3 <- c(as.character(info$V3[14]),as.character(info$V4[14]),as.character(info$V5[14]),as.character(info$V6[14]),as.character(info$V7[14]),as.character(info$V8[14]),as.character(info$V9[14]),as.character(info$V10[14]),as.character(info$V11[14]),as.character(info$V12[14])) 
    
    # If all true it means that the three experiments were ran with the lights switched on when platform moves right
    right_platform1 <- all(light_side1=="right")
    right_platform2 <- all(light_side2=="right")
    right_platform3 <- all(light_side3=="right")
    
    TimeExp <- data$t.s.[lengthExp]   # The total time it took for the experiment to complete
    data$Sampling<-c(0,diff(data$t.s., lag = 1)) # Calculating Inter Sample intervals (ISI)
    MaxSample<-max(data$Sampling)  # Checking what it the maximal ISI
    #plot(data[[2]]~data[[1]], type="l")  # This graph show how homogeneous was the sampling rate
    
    ######################################### Plot platform traces and decide which to keep for further analysis ############################################
    
    # Segment datapoint rechnen und plotten mit dem trace
    segment<- seq(from = 0,to = lengthExp, lengthExp/10)
    #sau_segment<- seq(from = 0,to = lengthExp, lengthExp/5) #saurabh  
    
    #plot(data$n,data$pos1, type = "l",xlab = "data points",ylab = "Position", main = "Platform 1")       ################# Amanda ####################
    #abline(v = segment, untf = FALSE, col="red",lwd=3)
    #message("please enter T for keeping the fly for analysis and F for deleting it and press enter")     ################# Amanda ####################
    #keep_fly1 <- scan(n=1,what= logical())
    
    #plot(data$n,data$pos2, type = "l",xlab = "data points",ylab = "Position", main = "Platform 2")       ################# Amanda ####################
    #abline(v = segment, untf = FALSE, col="red",lwd=3)                                                   ################# Amanda ####################
    #message("please enter T for keeping the fly for analysis and F for deleting it and press enter")          
    #keep_fly2 <- scan(n=1,what= logical())
    
    #plot(data$n,data$pos3, type = "l",xlab = "data points",ylab = "Position", main = "Platform 3")       ################# Amanda ####################
    #abline(v = segment, untf = FALSE, col="red",lwd=3)                                                   ################# Amanda ####################
    #message("please enter T for keeping the fly for analysis and F for deleting it and press enter")
    #keep_fly3 <- scan(n=1,what= logical())
    
    
    ############################################################### PI rechnen ##############################################
    
    
    
    
    
    # Save in state variable if the light is on or off taking care of hysteresis
    data$state1 <- ligth_state(data$pos1,Hysteresis)
    data$state2 <- ligth_state(data$pos2,Hysteresis)
    data$state3 <- ligth_state(data$pos3,Hysteresis)
    
    # Change the ON or OFF state if the platform were set to reinforce left
    if(right_platform1==FALSE){ data$state1 <- !data$state1}
    if(right_platform2==FALSE){ data$state2 <- !data$state2}
    if(right_platform3==FALSE){ data$state3 <- !data$state3}
    
    # Test Condition 1: if the platform signal stay too long flat, we want to discard results
    keep1<-all_false(flat_thres,data$pos1)
    keep2<-all_false(flat_thres,data$pos2)
    keep3<-all_false(flat_thres,data$pos3)
    
    ## Condition 2: Calculate PIs and sort the ones with good pretest
    
    # Allocate PI vectors
    PI_platform1 <- vector("numeric", length = 10)
    #combo_PI_platform1 <- vector("numeric", length = 5)
    PI_platform2 <- vector("numeric", length = 10)
    #combo_PI_platform2 <- vector("numeric", length = 5)
    PI_platform3 <- vector("numeric", length = 10)
    #combo_PI_platform3 <- vector("numeric", length = 5)
    
    # Allocate effect size vectors
    effectsize1 <- vector("numeric", 1)
    effectsize2 <- vector("numeric", 1)
    effectsize3 <- vector("numeric", 1)
    
    ## Calculate the wiggle as a measure of locomotory/activity state
    wiggle1<-abs(diff(data$pos1))
    wiggle2<-abs(diff(data$pos2))
    wiggle3<-abs(diff(data$pos3))
    
    ## Calculate the wiggle for each of the two treatments: lights on/off. 
    on_wiggle1 <- mean(wiggle1[data$state1[-24000]])
    off_wiggle1 <- mean(wiggle1[!data$state1[-24000]])
    on_wiggle2 <- mean(wiggle2[data$state2[-24000]])
    off_wiggle2 <- mean(wiggle2[!data$state2[-24000]])
    on_wiggle3 <- mean(wiggle3[data$state3[-24000]])
    off_wiggle3 <- mean(wiggle3[!data$state3[-24000]])
    
    # for(m in 2:length(data$n))
    # {
    #   if(data$state1[m]=="TRUE" && data$state1[m-1]== "TRUE")
    #   {
    #     on_wiggle1 <- on_wiggle1 + abs(data$pos1[m]- data$pos1[m-1])
    #   }
    # }
    # 
    # for(m in 2:length(data$n))
    # {
    #   if(data$state2[m]=="TRUE" && data$state2[m-1]== "TRUE")
    #   {
    #     on_wiggle2 <- on_wiggle2 + abs(data$pos2[m]- data$pos2[m-1])
    #   }
    # }
    # 
    # 
    # for(m in 2:length(data$n))
    # {
    #   if(data$state3[m]=="TRUE" && data$state3[m-1]== "TRUE")
    #   {
    #     on_wiggle3 <- on_wiggle3 + abs(data$pos3[m]- data$pos3[m-1])
    #   }
    # }
    # 
    # ###off_wiggle
    # for(m in 2:length(data$n))
    # {
    #   if(data$state1[m]=="FALSE" && data$state1[m-1]== "FALSE")
    #   {
    #     on_wiggle1 <- on_wiggle1 + abs(data$pos1[m]- data$pos1[m-1])
    #   }
    # }
    # 
    # for(m in 2:length(data$n))
    # {
    #   if(data$state2[m]=="FALSE" && data$state2[m-1]== "FALSE")
    #   {
    #     on_wiggle2 <- on_wiggle2 + abs(data$pos2[m]- data$pos2[m-1])
    #   }
    # }
    # 
    # 
    # for(m in 2:length(data$n))
    # {
    #   if(data$state3[m]=="FALSE" && data$state3[m-1]== "FALSE")
    #   {
    #     on_wiggle3 <- on_wiggle3 + abs(data$pos3[m]- data$pos3[m-1])
    #   }
    # }
    
    ## Apply the PI formula and save it
    for(oo in 1:10){
      PI_platform1[oo] <- calculate_PI(data$state1[segment[oo]:segment[oo+1]])
      PI_platform2[oo] <- calculate_PI(data$state2[segment[oo]:segment[oo+1]])
      PI_platform3[oo] <- calculate_PI(data$state3[segment[oo]:segment[oo+1]])
    }
    
    # for(oo in seq(1,11,2)){
    #   combo_PI_platform1[oo] <- calculate_PI(data$state1[segment[oo]:segment[oo+2]])
    #   combo_PI_platform2[oo] <- calculate_PI(data$state2[segment[oo]:segment[oo+2]])
    #   combo_PI_platform3[oo] <- calculate_PI(data$state3[segment[oo]:segment[oo+2]])
    # }
    
    ## Calculate the effect size by substracting the initial baseline (segments 1 and 2) to the training Periods (3,4,7,8)
    effectsize1 <- mean(PI_platform1[c(3,4,7,8)])-mean(PI_platform1[1:2])
    effectsize2 <- mean(PI_platform2[c(3,4,7,8)])-mean(PI_platform2[1:2])
    effectsize3 <- mean(PI_platform3[c(3,4,7,8)])-mean(PI_platform3[1:2])
    
    ## Get the PIs the training Periods (3,4,7,8) but without normalizing by initial baseline
    just_reinf1 <- mean(PI_platform1[c(3,4,7,8)])
    just_reinf2 <- mean(PI_platform2[c(3,4,7,8)])
    just_reinf3 <- mean(PI_platform3[c(3,4,7,8)])
    
    # Condition 3: we set a threshold of number of light switching to make sure the fly had the opportunity to explore both sides enough.
    # Tally light encounters in the first training period
    light_encounter1 <- sum(abs(diff(data$state1[(lengthExp/5):((lengthExp*2)/5)])))
    light_encounter2 <- sum(abs(diff(data$state2[(lengthExp/5):((lengthExp*2)/5)])))
    light_encounter3 <- sum(abs(diff(data$state3[(lengthExp/5):((lengthExp*2)/5)])))
    
    
    ## Here all conditions should be met: No flat signal, good baseline pretest PI and enough switching 
    if(keep1 && abs(mean(PI_platform1[1:2]))<PI_thres & light_encounter1>2){
      count<-count+1
      #pos <- pos + data$pos1
      PI_platform[((3*i)-2),] <- PI_platform1
      #combo_PI_platform[((3*i)-2),] <- combo_PI_platform1
      effectsize[((3*i)-2)] <- effectsize1
      effectsize_reinf[((3*i)-2)] <- just_reinf1
      on_wiggle[((3*i)-2)] <- on_wiggle1
      off_wiggle[((3*i)-2)] <- off_wiggle1
    }
    
    if(keep2 && abs(mean(PI_platform2[1:2]))<PI_thres & light_encounter2>2){
      count<-count+1
      #pos <- pos + data$pos2
      PI_platform[((3*i)-1),] <- PI_platform2
      #combo_PI_platform[((3*i)-1),] <- combo_PI_platform2 
      effectsize[((3*i)-1)] <- effectsize2
      effectsize_reinf[((3*i)-1)] <- just_reinf2
      on_wiggle[((3*i)-1)] <- on_wiggle2
      off_wiggle[((3*i)-1)] <- off_wiggle2
    }
    
    if(keep3 && abs(mean(PI_platform3[1:2]))<PI_thres & light_encounter3>2){
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
  
  ##  effectsize_mat[1:length(effectsize!=0.0),j] <- effectsize
  #  on_wiggle_mat[1:length(on_wiggle!=0.0),j]<-on_wiggle
  ##  off_wiggle_mat[1:length(off_wiggle!=0.0),j]<-off_wiggle
  effectsize_mat[1:length(effectsize),j] <- effectsize
  on_wiggle_mat[1:length(on_wiggle),j] <- on_wiggle
  off_wiggle_mat[1:length(off_wiggle),j] <- off_wiggle
  just_rein_mat[1:length(effectsize_reinf),j] <- effectsize_reinf
  all_PIs <- rbind(all_PIs,PI_platform)
  used_traces[j]<-count
  #effectsize <- vector("numeric", length = 3*Nexp_total) #resetting 
  #on_wiggle <- vector("numeric", length = 3*Nexp_total) #resetting
  #off_wiggle <- vector("numeric", length = 3*Nexp_total) #resetting
  
  
  #pos <- (pos)/count
  
  #plot(data$n, pos, type = "l",xlab = "data points",ylab = "mean_Position", main = "all_Platforms")
  #abline(v = segment, untf = FALSE, col="red",lwd=3)
  
  # saurabh effectsize and other parameters
  
  
  # mean_effectsize_mat <- apply(effectsize,2,function(effectsize){mean(effectsize, na.rm = TRUE)})
  #  std_dev_effectsize_mat <- apply(effectsize,2,function(effectsize){sd(effectsize, na.rm = TRUE)})
  
  
  
  # std_err_effectsize_mat <- std_dev_effectsize/count**0.5
  
  
  #barplot(mean_effectsize)
  #ggplot
  
  # Boxplot of the PIs
  boxplot(PI_platform, col="grey",xlab="",ylab="PI",main=group_name, ylim = c(-1, 1),names=c("Pretest","Pretest","Training","Training","Test","Test","Training","Training","Test","Test"), cex.lab=1.5, cex.axis = 1.2)
  abline(h = 0, untf = FALSE, col="black",lwd=3)
  group <- NULL
  for(s in 1:10){
    a <- rep(s,Nexp*3)
    group <- append(group,a)
  }
  stripchart(as.vector(PI_platform)~group,vertical = TRUE, method = "jitter",pch = 21, col = "maroon", bg = "bisque",add = TRUE) 
  
  accepted_flies <- which(!is.na(PI_platform[,1]))
  
  matplot(t(PI_platform[accepted_flies,]), type = c("b"),xlab=group_name,main=group_name, pch=1,col = 1:4)
  
}


#colnames(effectsize_mat) = c("blue1","blue3","blue4","darkgreen1","darkgreen2","darkgreen3","darkpink1","darkpink3","darkpink4","darkyellow1","darkyellow2","darkyellow3","darkyellow4","lightgreen1","lightgreen2","lightgreen4","lightyellow2","orange1","orange3","orange4","red1","red2","red3","white1","white2","positivecontrol")
colnames(effectsize_mat) <- as.character(all_screens)
#colnames(just_rein_mat) = c("blue1","blue3","blue4","darkgreen1","darkgreen2","darkgreen3","darkpink1","darkpink3","darkpink4","darkyellow1","darkyellow2","darkyellow3","darkyellow4","lightgreen1","lightgreen2","lightgreen4","lightyellow2","orange1","orange3","orange4","red1","red2","red3","white1","white2","positivecontrol")
colnames(just_rein_mat) <- as.character(all_screens)

#colnames(diff_wiggle_mat) = c("blue1","blue3","blue4","darkgreen1","darkgreen2","darkgreen3","darkpink1","darkpink3","darkpink4","darkyellow1","darkyellow2","darkyellow3","darkyellow4","lightgreen1","lightgreen2","lightgreen4","lightyellow2","orange1","orange3","orange4","red1","red2","red3","white1","white2","positivecontrol")
#colnames(diff_wiggle_mat) <- as.character(all_screens)

#boxplot(effectsize_mat,col="yelow",xlab="",ylab="reinforcement", names = c("blue1","blue3","blue4","darkgreen1","darkgreen2","darkgreen3","darkpink1","darkpink3","darkpink4","darkyellow1","darkyellow2","darkyellow3","darkyellow4","lightgreen1","lightgreen2","lightgreen4","lightyellow2","orange1","orange3","orange4","red1","red2","red3","white1","white2","positivecontrol"), cex.lab=1.0, cex.axis = 0.5, las=2)
boxplot(effectsize_mat,col="yellow",xlab="",ylab="reinforcement - pretest", names = as.character(all_screens), cex.lab=1.0, cex.axis = 1, las=2)
abline(h = 0, untf = FALSE, col="black",lwd=3)

#boxplot(effectsize_mat,col="blue",xlab="",ylab="reinforcement",names=c("red1","red2","darkyellow2","blue3","blue4","darkpink2","darkpink3","darkpink4","orange1","orange3","lightgreen1","positivecontrol"), cex.lab=1.0, cex.axis = 0.5, las=2)
#abline(h = 0, untf = FALSE, col="black",lwd=3)

#boxplot(just_rein_mat,col="yellow",xlab="",ylab="just_reinforcement",c("blue1","blue3","blue4","darkgreen1","darkgreen2","darkgreen3","darkpink1","darkpink3","darkpink4","darkyellow1","darkyellow2","darkyellow3","darkyellow4","lightgreen1","lightgreen2","lightgreen4","lightyellow2","orange1","orange3","orange4","red1","red2","red3","white1","white2","positivecontrol"), cex.lab=1.0, cex.axis = 0.5, las=2)
boxplot(just_rein_mat,col="yellow",xlab="",ylab="just_reinforcement", names = as.character(all_screens), cex.lab=1.0, cex.axis = 1, las=2)
abline(h = 0, untf = FALSE, col="black",lwd=3)

diff_wiggle_mat=on_wiggle_mat-off_wiggle_mat
#boxplot(diff_wiggle_mat,col="yellow",xlab="",ylab="wiggle difference(on - off)", c("blue1","blue3","blue4","darkgreen1","darkgreen2","darkgreen3","darkpink1","darkpink3","darkpink4","darkyellow1","darkyellow2","darkyellow3","darkyellow4","lightgreen1","lightgreen2","lightgreen4","lightyellow2","orange1","orange3","orange4","red1","red2","red3","white1","white2","positivecontrol"),cex.lab=1.0, cex.axis = 0.5, las = 2)
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

#barplot(mean_effectsize,col="yellow",xlab="",ylab="reinforcement",main = "Amanda", cex.lab=1.0, cex.axis = 0.5)
#abline(h = 0, untf = FALSE, col="black",lwd=3)

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


#barplot(mean_diff_wiggle,col="yellow",xlab="",ylab="wiggle difference(on - off)", cex.lab=1.0, cex.axis = 0.5)
#abline(h = 0, untf = FALSE, col="black",lwd=3)


## Total PI overall and number of experiments per line

boxplot(all_PIs,main="overall PI for the whole screen - mixed")

names(used_traces)<-all_screens
barplot(used_traces,las=2,main = "mixed experiments")


####################################
=======
rm(list=ls()) #clean memory
gc()          #collect garbage

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
           label=paste(samplesizes[boxes]))
}


################################################################ Set working directory where I saved the data   ##################################################

#setwd("C:/Users/LocalAdmin/Desktop/new_screens/all")
#setwd("C:/Users/LocalAdmin/Desktop/new_screens/AMANDA/amanda_all")
#setwd("C:/Users/LocalAdmin/Desktop/new_screens/SAURABH/saurabh_all")
setwd("D:/data/optogenetics/others")
codedir="B:/GitHub/screen-analysis-Joystick"
########################################## Initialize some parameters  ###########################################


#tested_flies <- read.table("C:/Users/LocalAdmin/Desktop/new_screens/all/all.txt", quote="\"", comment.char="")
#tested_flies <- read.table("C:/Users/LocalAdmin/Desktop/new_screens/AMANDA/amanda_all/amanda_all.txt", quote="\"", comment.char="")
#tested_flies <- read.table("C:/Users/LocalAdmin/Desktop/new_screens/SAURABH/saurabh_all/saurabh_all.txt", quote="\"", comment.char="")
tested_flies <- read.table("D:/data/optogenetics/others/all.txt", quote="\"", comment.char="")
all_screens <- unique(tested_flies$V2)
no_of_screens <- length(all_screens)
skip<-37
flat_thres <- 0.8
PI_thres <- 1 #flies rejected for pretest PI=1

groupedPIs <- list()

#j <- 0 #initializing the screen number variable
#x<-0
on_wiggle1<- NA
on_wiggle2<- NA
on_wiggle3 <- NA


off_wiggle1<- NA
off_wiggle2<- NA
off_wiggle3 <- NA


#message("please enter the number of experiments")
#Nexp <- scan(n=1,what= numeric())   #Number of flies

Nexp_total <- 40 # they won´t be more experiment than these #round(length(tested_flies$V1)/no_of_screens)*4


#effectsize <- as.numeric(rep(NA,3*Nexp_total)) # Saurabh variable defining effectsize score
effectsize_mat <- matrix(NA, 3*Nexp_total, no_of_screens)
just_rein_mat <- matrix(NA, 3*Nexp_total, no_of_screens)
used_traces <- as.numeric(rep(NA,no_of_screens))

#on_wiggle <- as.numeric(rep(NA,3*Nexp_total))
#off_wiggle <- as.numeric(rep(NA,3*Nexp_total))
on_wiggle_mat <- matrix(NA, 3*Nexp_total, no_of_screens)
off_wiggle_mat <- matrix(NA, 3*Nexp_total, no_of_screens)
all_PIs <- as.numeric(rep(NA,10))

for(j in 1:no_of_screens)
  
{
  #pos <- 0 #saurabhdefining variable for mean position
  count <- 0 #Saurabh defining a variable that keeps a count of the number of graphs taken 
  
  Nexp <- length(tested_flies$V1[tested_flies$V2==all_screens[j]])
  data_files <- as.character(tested_flies$V1[tested_flies$V2==all_screens[j]])
  #effectsize <- as.numeric(rep(NA,3*Nexp))
  #on_wiggle <- as.numeric(rep(NA,3*Nexp))
  #off_wiggle <- as.numeric(rep(NA,3*Nexp))
  #### This seems to me not well coded. By Christian Rohrsen
  
  # if(all_screens[j]=="blue1")
  # {
  #   Nexp <- length(tested_flies$V1[tested_flies$V2=="blue1"])
  #   x <- 0
  #   
  # }
  # 
  # if(all_screens[j]=="blue3")
  # {
  #   Nexp <- length(tested_flies$V1[tested_flies$V2=="blue3"])
  #   x <- length(tested_flies$V1[tested_flies$V2=="blue1"]) +x
  #   
  # }
  # 
  # if(all_screens[j]=="blue4")
  # {
  #   Nexp <- length(tested_flies$V1[tested_flies$V2=="blue4"])
  #   x <- length(tested_flies$V1[tested_flies$V2=="blue3"]) +x
  #   
  # }
  # 
  # if(all_screens[j]=="darkgreen1")
  # {
  #   
  #   Nexp <- length(tested_flies$V1[tested_flies$V2=="darkgreen1"])
  #   x <- length(tested_flies$V1[tested_flies$V2=="blue4"])+x
  #   
  # }
  # 
  # if(all_screens[j]=="darkgreen2")
  # {
  #   
  #   Nexp <- length(tested_flies$V1[tested_flies$V2=="darkgreen2"])
  #   x <- length(tested_flies$V1[tested_flies$V2=="darkgreen1"])+x
  #   
  # }
  # 
  # if(all_screens[j]=="darkgreen3")
  # {
  #   
  #   Nexp <- length(tested_flies$V1[tested_flies$V2=="darkgreen3"])
  #   x <- length(tested_flies$V1[tested_flies$V2=="darkgreen2"])+x
  #   
  # }
  # 
  # if(all_screens[j]== "darkpink1")
  # {
  #   
  #   Nexp <- length(tested_flies$V1[tested_flies$V2=="darkpink1"])
  #   x <- length(tested_flies$V1[tested_flies$V2=="darkgreen3"])+x
  #   
  # }
  # 
  # if(all_screens[j]=="darkpink2")
  # {
  #   
  #   Nexp <- length(tested_flies$V1[tested_flies$V2=="darkpink2"])
  #   x <- length(tested_flies$V1[tested_flies$V2== "darkpink2"])+x
  #   
  # }
  # 
  # 
  # if(all_screens[j]=="darkpink3")
  # {
  #   
  #   
  #   Nexp <- length(tested_flies$V1[tested_flies$V2=="darkpink3"])
  #   x <- length(tested_flies$V1[tested_flies$V2=="darkpink2"])+x
  #   
  # }
  # 
  # if(all_screens[j]=="darkpink4")
  # {
  #   
  #   
  #   Nexp <- length(tested_flies$V1[tested_flies$V2=="darkpink4"])
  #   x <- length(tested_flies$V1[tested_flies$V2=="darkpink3"])+x
  #   
  # }
  # 
  # if(all_screens[j]=="darkyellow1")
  # {
  #   
  #   Nexp <- length(tested_flies$V1[tested_flies$V2=="darkyellow1"])
  #   x <- length(tested_flies$V1[tested_flies$V2=="darkpink4"])+x
  #   
  # }
  # 
  # if(all_screens[j]=="darkyellow2")
  # {
  #   
  #   
  #   Nexp <- length(tested_flies$V1[tested_flies$V2=="darkyellow2"])
  #   x <- length(tested_flies$V1[tested_flies$V2=="darkyellow1"])+x
  #   
  # }
  # 
  # if(all_screens[j]=="darkyellow3")
  # {
  #   
  #   
  #   Nexp <- length(tested_flies$V1[tested_flies$V2=="darkyellow3"])
  #   x <- length(tested_flies$V1[tested_flies$V2=="darkyellow2"])+x
  #   
  # }
  # 
  # if(all_screens[j]=="darkyellow4")
  # {
  #   
  #   
  #   Nexp <- length(tested_flies$V1[tested_flies$V2=="darkyellow4"])
  #   x <- length(tested_flies$V1[tested_flies$V2=="darkyellow3"])+x
  #   
  # }
  # 
  # if(all_screens[j]=="lightgreen1")
  # {
  #   
  #   
  #   Nexp <- length(tested_flies$V1[tested_flies$V2=="lightgreen1"])
  #   x <- length(tested_flies$V1[tested_flies$V2=="darkyellow4"])+x
  #   
  # }
  # 
  # if(all_screens[j]=="lightgreen2")
  # {
  #   
  #   
  #   Nexp <- length(tested_flies$V1[tested_flies$V2=="lightgreen2"])
  #   x <- length(tested_flies$V1[tested_flies$V2=="lightgreen1"])+x
  #   
  # }
  # 
  # if(all_screens[j]=="lightgreen4")
  # {
  #   
  #   
  #   Nexp <- length(tested_flies$V1[tested_flies$V2=="lightgreen4"])
  #   x <- length(tested_flies$V1[tested_flies$V2=="lightgreen2"])+x
  #   
  # }
  # 
  # if(all_screens[j]=="lightyellow2")
  # {
  #   
  #   
  #   Nexp <- length(tested_flies$V1[tested_flies$V2=="lightyellow2"])
  #   x <- length(tested_flies$V1[tested_flies$V2=="lightgreen4"])+x
  #   
  # }
  # 
  # if(all_screens[j]=="orange1")
  # {
  #   
  #   
  #   Nexp <- length(tested_flies$V1[tested_flies$V2=="orange1"])
  #   x <- length(tested_flies$V1[tested_flies$V2=="lightyellow2"])+x
  #   
  # }
  # 
  # if(all_screens[j]=="orange2")
  # {
  #   
  #   
  #   Nexp <- length(tested_flies$V1[tested_flies$V2=="orange2"])
  #   x <- length(tested_flies$V1[tested_flies$V2=="orange1"])+x
  #   
  # }
  # 
  # if(all_screens[j]=="orange3")
  # {
  #   
  #   
  #   Nexp <- length(tested_flies$V1[tested_flies$V2=="orange3"])
  #   x <- length(tested_flies$V1[tested_flies$V2=="orange2"])+x
  #   
  # }
  # 
  # if(all_screens[j]=="orange4")
  # {
  #   
  #   
  #   Nexp <- length(tested_flies$V1[tested_flies$V2=="orange4"])
  #   x <- length(tested_flies$V1[tested_flies$V2=="orange3"])+x
  #   
  # }
  # 
  # if(all_screens[j]=="red1")
  # {
  #   
  #   
  #   Nexp <- length(tested_flies$V1[tested_flies$V2=="red1"])
  #   x <- length(tested_flies$V1[tested_flies$V2=="orange4"])+x
  #   
  # }
  # 
  # 
  # if(all_screens[j]=="red2")
  # {
  #   
  #   
  #   Nexp <- length(tested_flies$V1[tested_flies$V2=="red2"])
  #   x <- length(tested_flies$V1[tested_flies$V2=="red1"])+x
  #   
  # }
  # 
  # if(all_screens[j]=="red3")
  # {
  #   
  #   
  #   Nexp <- length(tested_flies$V1[tested_flies$V2=="red3"])
  #   x <- length(tested_flies$V1[tested_flies$V2=="red2"])+x
  #   
  # }
  # 
  # if(all_screens[j]=="white1")
  # {
  #   
  #   
  #   Nexp <- length(tested_flies$V1[tested_flies$V2=="white1"])
  #   x <- length(tested_flies$V1[tested_flies$V2=="red3"])+x
  #   
  # }
  # 
  # if(all_screens[j]=="white2")
  # {
  #   
  #   
  #   Nexp <- length(tested_flies$V1[tested_flies$V2=="white2"])
  #   x <- length(tested_flies$V1[tested_flies$V2=="white1"])+x
  #   
  # }
  # 
  # if(all_screens[j]=="positivecontrol")
  # {
  #   
  #   
  #   Nexp <- length(tested_flies$V1[tested_flies$V2=="positivecontrol1"])
  #   x <- length(tested_flies$V1[tested_flies$V2=="white2"])+x
  #   
  # }
  
  
  
  PI_platform <- matrix(NA, 3*Nexp, 10)    # Variable where PIs are saved
  combo_PI_platform <- matrix(NA, 3*Nexp, 5) # Saurabh variable where the mean of PI's are stored for 2-2 segments
  
  
  on_wiggle <- as.numeric(rep(NA,3*Nexp))
  off_wiggle <- as.numeric(rep(NA,3*Nexp))
  diff_wiggle <- as.numeric(rep(NA,3*Nexp))
  
  effectsize <- as.numeric(rep(NA,3*Nexp)) # Saurabh variable defining effectsize score
  effectsize_reinf <- as.numeric(rep(NA,3*Nexp)) # Saurabh variable defining effectsize score

  ######saurabh defining wiggle
  group_name <- all_screens[j]
  
  # Start a for loop for the number of files to analyze
  for(i in 1:Nexp){
    #group_name <- tested_flies$V2[i]   #Name of the group
    
    
    ############################## Import in a dataframe just the values. ################################
    
    #data <- read.table(file.choose(), header = FALSE, sep = "\t", quote = "\"" , dec = ".", fill = TRUE, skip = skip , comment.char = "", nrows = 24037-skip,col.names=c("n","t.s.","pos1","pos2","pos3"))
    
    
    #data <- read.table(as.character(tested_flies$V1[x+i]), header = FALSE, sep = "\t", quote = "\"" , dec = ".", fill = TRUE, skip = skip , comment.char = "", nrows = 24037-skip,col.names=c("n","t.s.","pos1","pos2","pos3"))
    #group_name <- tested_flies$V2[x+i]   #Name of the group
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
    #plot(data[[2]]~data[[1]], type="l")  # This graph show how homogeneous was the sampling rate
    
    ######################################### Plot platform traces and decide which to keep for further analysis ############################################
    
    # Segment rechnen und plotten mit dem trace
    
    segment<- seq(from = 0,to = lengthExp, lengthExp/10)
    #sau_segment<- seq(from = 0,to = lengthExp, lengthExp/5) #saurabh  
    
    #plot(data$n,data$pos1, type = "l",xlab = "data points",ylab = "Position", main = "Platform 1")       ################# Amanda ####################
    #abline(v = segment, untf = FALSE, col="red",lwd=3)
    #message("please enter T for keeping the fly for analysis and F for deleting it and press enter")     ################# Amanda ####################
    #keep_fly1 <- scan(n=1,what= logical())
    
    #plot(data$n,data$pos2, type = "l",xlab = "data points",ylab = "Position", main = "Platform 2")       ################# Amanda ####################
    #abline(v = segment, untf = FALSE, col="red",lwd=3)                                                   ################# Amanda ####################
    #message("please enter T for keeping the fly for analysis and F for deleting it and press enter")          
    #keep_fly2 <- scan(n=1,what= logical())
    
    #plot(data$n,data$pos3, type = "l",xlab = "data points",ylab = "Position", main = "Platform 3")       ################# Amanda ####################
    #abline(v = segment, untf = FALSE, col="red",lwd=3)                                                   ################# Amanda ####################
    #message("please enter T for keeping the fly for analysis and F for deleting it and press enter")
    #keep_fly3 <- scan(n=1,what= logical())
    
    
    ############################################################### PI rechnen ##############################################
    
    
    
    
    
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
    
    # for(m in 2:length(data$n))
    # {
    #   if(data$state1[m]=="TRUE" && data$state1[m-1]== "TRUE")
    #   {
    #     on_wiggle1 <- on_wiggle1 + abs(data$pos1[m]- data$pos1[m-1])
    #   }
    # }
    # 
    # for(m in 2:length(data$n))
    # {
    #   if(data$state2[m]=="TRUE" && data$state2[m-1]== "TRUE")
    #   {
    #     on_wiggle2 <- on_wiggle2 + abs(data$pos2[m]- data$pos2[m-1])
    #   }
    # }
    # 
    # 
    # for(m in 2:length(data$n))
    # {
    #   if(data$state3[m]=="TRUE" && data$state3[m-1]== "TRUE")
    #   {
    #     on_wiggle3 <- on_wiggle3 + abs(data$pos3[m]- data$pos3[m-1])
    #   }
    # }
    # 
    # ###off_wiggle
    # for(m in 2:length(data$n))
    # {
    #   if(data$state1[m]=="FALSE" && data$state1[m-1]== "FALSE")
    #   {
    #     on_wiggle1 <- on_wiggle1 + abs(data$pos1[m]- data$pos1[m-1])
    #   }
    # }
    # 
    # for(m in 2:length(data$n))
    # {
    #   if(data$state2[m]=="FALSE" && data$state2[m-1]== "FALSE")
    #   {
    #     on_wiggle2 <- on_wiggle2 + abs(data$pos2[m]- data$pos2[m-1])
    #   }
    # }
    # 
    # 
    # for(m in 2:length(data$n))
    # {
    #   if(data$state3[m]=="FALSE" && data$state3[m-1]== "FALSE")
    #   {
    #     on_wiggle3 <- on_wiggle3 + abs(data$pos3[m]- data$pos3[m-1])
    #   }
    # }
    
    
    for(oo in 1:10){
      PI_platform1[oo] <- calculate_PI(data$state1[segment[oo]:segment[oo+1]])
      PI_platform2[oo] <- calculate_PI(data$state2[segment[oo]:segment[oo+1]])
      PI_platform3[oo] <- calculate_PI(data$state3[segment[oo]:segment[oo+1]])
    }
    
    # for(oo in seq(1,11,2)){
    #   combo_PI_platform1[oo] <- calculate_PI(data$state1[segment[oo]:segment[oo+2]])
    #   combo_PI_platform2[oo] <- calculate_PI(data$state2[segment[oo]:segment[oo+2]])
    #   combo_PI_platform3[oo] <- calculate_PI(data$state3[segment[oo]:segment[oo+2]])
    # }
    
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
  
  ##  effectsize_mat[1:length(effectsize!=0.0),j] <- effectsize
  #  on_wiggle_mat[1:length(on_wiggle!=0.0),j]<-on_wiggle
  ##  off_wiggle_mat[1:length(off_wiggle!=0.0),j]<-off_wiggle
  effectsize_mat[1:length(effectsize),j] <- effectsize
  on_wiggle_mat[1:length(on_wiggle),j] <- on_wiggle
  off_wiggle_mat[1:length(off_wiggle),j] <- off_wiggle
  just_rein_mat[1:length(effectsize_reinf),j] <- effectsize_reinf
  all_PIs <- rbind(all_PIs,PI_platform)
  used_traces[j]<-count
  #effectsize <- vector("numeric", length = 3*Nexp_total) #resetting 
  #on_wiggle <- vector("numeric", length = 3*Nexp_total) #resetting
  #off_wiggle <- vector("numeric", length = 3*Nexp_total) #resetting
  
  
  #pos <- (pos)/count
  
  #plot(data$n, pos, type = "l",xlab = "data points",ylab = "mean_Position", main = "all_Platforms")
  #abline(v = segment, untf = FALSE, col="red",lwd=3)
  
  # saurabh effectsize and other parameters
  
  
  # mean_effectsize_mat <- apply(effectsize,2,function(effectsize){mean(effectsize, na.rm = TRUE)})
  #  std_dev_effectsize_mat <- apply(effectsize,2,function(effectsize){sd(effectsize, na.rm = TRUE)})
  
  
  
  # std_err_effectsize_mat <- std_dev_effectsize/count**0.5
  
  
  #barplot(mean_effectsize)
  #ggplot
  
  # Boxplot of the PIs
  boxplot(PI_platform, col="grey",xlab="",ylab="PI",main=group_name, ylim = c(-1, 1),names=c("Test","Test","Training","Training","Test","Test","Training","Training","Test","Test"), cex.lab=1.5, cex.axis = 1.2)
  abline(h = 0, untf = FALSE, col="black",lwd=3)
  group <- NULL
  for(s in 1:10){
    a <- rep(s,Nexp*3)
    group <- append(group,a)
  }
  stripchart(as.vector(PI_platform)~group,vertical = TRUE, method = "jitter",pch = 21, col = "maroon", bg = "bisque",add = TRUE) 
  
  accepted_flies <- which(!is.na(PI_platform[,1]))
  
  matplot(t(PI_platform[accepted_flies,]), type = c("b"),xlab=group_name,main=group_name, pch=1,col = 1:4)

colnames(PI_platform)<-c("pretest1","pretest2","training1","training2","test1","test2","training3","training4","test3","test4")
groupedPIs[[j]]=as.data.frame(na.omit(PI_platform)) #store the PIs without the excluded flies
    
}  #for nofScreens

names(groupedPIs) <- as.character(all_screens) #name the list elements in the list if PIs

#colnames(effectsize_mat) = c("blue1","blue3","blue4","darkgreen1","darkgreen2","darkgreen3","darkpink1","darkpink3","darkpink4","darkyellow1","darkyellow2","darkyellow3","darkyellow4","lightgreen1","lightgreen2","lightgreen4","lightyellow2","orange1","orange3","orange4","red1","red2","red3","white1","white2","positivecontrol")
colnames(effectsize_mat) <- as.character(all_screens)
#colnames(just_rein_mat) = c("blue1","blue3","blue4","darkgreen1","darkgreen2","darkgreen3","darkpink1","darkpink3","darkpink4","darkyellow1","darkyellow2","darkyellow3","darkyellow4","lightgreen1","lightgreen2","lightgreen4","lightyellow2","orange1","orange3","orange4","red1","red2","red3","white1","white2","positivecontrol")
colnames(just_rein_mat) <- as.character(all_screens)

#colnames(diff_wiggle_mat) = c("blue1","blue3","blue4","darkgreen1","darkgreen2","darkgreen3","darkpink1","darkpink3","darkpink4","darkyellow1","darkyellow2","darkyellow3","darkyellow4","lightgreen1","lightgreen2","lightgreen4","lightyellow2","orange1","orange3","orange4","red1","red2","red3","white1","white2","positivecontrol")
#colnames(diff_wiggle_mat) <- as.character(all_screens)

#boxplot(effectsize_mat,col="yelow",xlab="",ylab="reinforcement", names = c("blue1","blue3","blue4","darkgreen1","darkgreen2","darkgreen3","darkpink1","darkpink3","darkpink4","darkyellow1","darkyellow2","darkyellow3","darkyellow4","lightgreen1","lightgreen2","lightgreen4","lightyellow2","orange1","orange3","orange4","red1","red2","red3","white1","white2","positivecontrol"), cex.lab=1.0, cex.axis = 0.5, las=2)
boxplot(effectsize_mat,col="yellow",xlab="",ylab="reinforcement - pretest", names = as.character(all_screens), cex.lab=1.0, cex.axis = 1, las=2)
abline(h = 0, untf = FALSE, col="black",lwd=3)

#boxplot(effectsize_mat,col="blue",xlab="",ylab="reinforcement",names=c("red1","red2","darkyellow2","blue3","blue4","darkpink2","darkpink3","darkpink4","orange1","orange3","lightgreen1","positivecontrol"), cex.lab=1.0, cex.axis = 0.5, las=2)
#abline(h = 0, untf = FALSE, col="black",lwd=3)

#boxplot(just_rein_mat,col="yellow",xlab="",ylab="just_reinforcement",c("blue1","blue3","blue4","darkgreen1","darkgreen2","darkgreen3","darkpink1","darkpink3","darkpink4","darkyellow1","darkyellow2","darkyellow3","darkyellow4","lightgreen1","lightgreen2","lightgreen4","lightyellow2","orange1","orange3","orange4","red1","red2","red3","white1","white2","positivecontrol"), cex.lab=1.0, cex.axis = 0.5, las=2)
boxplot(just_rein_mat,col="yellow",xlab="",ylab="just_reinforcement", names = as.character(all_screens), cex.lab=1.0, cex.axis = 1, las=2)
abline(h = 0, untf = FALSE, col="black",lwd=3)

diff_wiggle_mat=on_wiggle_mat-off_wiggle_mat
#boxplot(diff_wiggle_mat,col="yellow",xlab="",ylab="wiggle difference(on - off)", c("blue1","blue3","blue4","darkgreen1","darkgreen2","darkgreen3","darkpink1","darkpink3","darkpink4","darkyellow1","darkyellow2","darkyellow3","darkyellow4","lightgreen1","lightgreen2","lightgreen4","lightyellow2","orange1","orange3","orange4","red1","red2","red3","white1","white2","positivecontrol"),cex.lab=1.0, cex.axis = 0.5, las = 2)
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

#barplot(mean_effectsize,col="yellow",xlab="",ylab="reinforcement",main = "Amanda", cex.lab=1.0, cex.axis = 0.5)
#abline(h = 0, untf = FALSE, col="black",lwd=3)

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


#barplot(mean_diff_wiggle,col="yellow",xlab="",ylab="wiggle difference(on - off)", cex.lab=1.0, cex.axis = 0.5)
#abline(h = 0, untf = FALSE, col="black",lwd=3)


## Total PI overall and number of experiments per line

boxplot(all_PIs,main="overall PI for the whole screen - mixed")

names(used_traces)<-all_screens
barplot(used_traces,las=2,main = "mixed experiments")


####################################

#### call RMarkdown to generate HTML file with evaluations #################################
rmarkdown::render(paste(codedir,"/firstscreen.Rmd", sep=""), 
                  output_file = "firstscreen.html", 
                  output_dir = "D:/data/optogenetics")
#### end RMarkdown for project evaluations #################################################
>>>>>>> Stashed changes:firstscreen.R
