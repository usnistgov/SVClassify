
### One Class Classification Script. Prepared for the svclassify Paper
### April 6 2015.  Hari Iyer, Statistical Engineering Division, ITL, NIST


### Step 1: Data Input and Data Preparation for L1 - One Class Classification

### Modify the path for home to suit your situation
rm(list = ls())	# Clear the work space. Desu
library(e1071)	# For SVM
home <- "/svclassify/"                          # Base directory for the project
tech <- c("ILL250","Moleculo","PacBio","PlatGen") #root names for technologies
### Make sure that this base directory has the following subdirectories

### Code:         /svclassify/code
### Results:      /svclassify/results
### Data:         /svclassify/data
###       Under the Data directory you need to have one subdirectory for each technology.
###       In our example we will have the following four subdirectories
###               /svclassify/data/ILL250
###               /svclassify/data/Moleculo
###               /svclassify/data/PacBio
###               /svclassify/data/PlatGen
### Inputs:       /svclassify/inputs   



####################################
### READ THE FOLLOWING explanations
####################################

###       Make sure you have installed the following packages:
###       (1) proxy (2) dendroextras (3) flashClust

###       In the input directory you should have a file called "filenames.txt" which has
###          a list of names of raw data files that should be included in the analyses. 
###          See the example provided in this directory.
###       In the data directory, under the subdirectory for each technology, store the raw data files
###          for each source of data (Randoms, PersonalisGold, etc). 
###          See example provided in these subdirectories
###       In the code directory the script files are provided. You do not need to change any of this
###          unless you want to get inside the code and make changes (at your own risk)
###       In the results directory you will find graph files and text files and csv files containing
###          various outputs, especially those discussed in the svclassify paper.


###########################################################################################################
########   DO NOT CHANGE ANYTHING ELSE  ###################################################################

print("Calculations Have Begun");flush.console()

codehome <- paste(home,"code/",sep="")            # Subdirectory where all code is stored
resultshome <- paste(home,"results/",sep="")      # Subdirectory where results files will be saved
datahome <- paste(home,"data/",sep="")            # Subdirectory containing additional subdirectories
                                                     # one subsubdirectory for each technology
                                                     # Data files are tab delimited text files
inputshome <- paste(home,"inputs/",sep="")        # Subdirectory containing INPUT information files

ntech <- length(tech)                             # number of technologies
nstrat <- 2*ntech+1                               # number of classification strategies

filenames <- read.table(paste(home,"inputs/fileNames.txt",sep=""),as.is=TRUE)
filenames <- c(as.matrix(filenames))              # converting into a character vector of file names
numfiles <- length(filenames)

source(paste(codehome,"myFunctions.R",sep=""))


### Create RData files from Hemant's text files
### Under the "Data" directory (datahome) there should be one subdirectory for each technology
### The subdirectory names should be exactly as the technology names in the vector tech above

### Create generic directory/filename for each technology, each set
fnames <- rep(NA,numfiles)
for (i in 1:numfiles){
fnames[i] <- paste(datahome,"FILE/",filenames[i],sep="")
}

fsave <- paste(datahome,"FILE.RData",sep="")

### SET names
SET = rep(NA,numfiles)
for (i in 1:numfiles){
temp <- unlist(strsplit(fnames[i],"FILE/"))[2]
SET[i] <- unlist(strsplit(temp,"_FILE"))[1]
}

strats = rep(NA,nstrat)
for (itech in 1:ntech){
strats[itech]=tech[itech]
}
strats[ntech+1]="Combined"
for (itech in 1:ntech){
strats[ntech+1+itech]= paste(itech," or More",sep="")
}



### READING the tab delimited text files and creating 
  # corresponding .RData files
print("Reading Data");flush.console()

for (i in 1:ntech){             # i loop start
#
for (j in 1:numfiles){          # j loop start
fname <- gsub("FILE",tech[i],fnames[j])
print(fname);flush.console()
temp <- read.table(fname,header=TRUE,as.is=TRUE,sep="\t")
if (j==1) cnames <- colnames(temp)
if (j >1) temp <- temp[cnames]
temp$GROUP <- j
if (j == 1) dat <- temp else dat <- rbind(dat,temp)
                      }         # j loop end
#
fname <- gsub("FILE",tech[i],fsave)
print(fname);flush.console()
save(dat,file=fname)
}                               # i loop end


###################################################################################

### Read all the input information for each technology
  # such as which annotations to use
  # (and, in future versions, which tail is relevant)
print("Reading Input Information");flush.console()

inputInfo <- NULL
infile <- "Annotations_FILE.csv"
for (i in 1:ntech){
readinput <- gsub("FILE",tech[i],infile)
temp <- read.table(paste(inputshome,readinput,sep=""),header=FALSE,skip=1,sep=",",as.is=TRUE)
colnames(temp) <- c("Annotation","OneClass","Unsupervised","Tail")
inputInfo <- c(inputInfo,list(temp))
}

###################################################################################

### Read and assemble the selected data for all technologies
  # Delete rows with missing values
  # Delete columns that are constant (standard deviation = 0)
print("Creating RData files");flush.console()

datpAll <- NULL

for (itech in 1:ntech){                  # itech loop
dname = gsub("FILE",tech[itech],fsave)
load(dname)
inInfo = inputInfo[[itech]]
idv = which(inInfo[,2]==1)
keepZ = inInfo[idv,1]
temp <- dat[,c("GROUP",keepZ)]          

#### Remove constant columns
 idsd <- which(apply(temp,2,sd,na.rm=TRUE)==0)
 if (length(idsd)>0) temp <- temp[,-idsd]

#### Remove rows with missing values in the selected variables
 idmiss <- unique(which(is.na(temp)==TRUE,arr.ind=TRUE)[,1])
 if (length(idmiss) > 0) datp <- temp[-idmiss,] else datp <- temp
 rm(temp)

datpAll <- c(datpAll,list(datp))
}                                        # end loop itech

save(datpAll,file=paste(datahome,"datpAllForOneClass.RData",sep=""))

###############################################################################

#### datpAll reloaded. For testing and debugging one may start here if datpAll is
   # already created and stored.
load(file=paste(datahome,"datpAllForOneClass.RData",sep=""))   

#### Form the combined data from all technologies
print("Combining Data from All Technologies");flush.console()

xmc <- datpAll[[1]][,-1]
if (ntech > 1){                           # start 'if'
 for (itech in 2:ntech){                  # start itech loop
 xmc <- cbind(xmc,datpAll[[itech]][,-1])
}                                         # end itech loop
}                                         # end 'if'
xmc_temp <- xmc
#### Precautionary processing in case same column occurs in two or more technologies
cxmc <- cor(xmc)
diag(cxmc) <- 0
idc <- unique(which(cxmc==1,arr.ind=TRUE)[,1])
if (length(idc)>0) {
   idcm <- idc[2:length(idc)]
   xmc <- xmc[,-idcm] }
####

## Now ROC for Individual Technologies
print("Starting ROC Calculations");flush.console()
GROUP = datpAll[[1]]$GROUP
x=seq(0,1,0.001)                             # quantile points 
nq=length(x)

nr <- nrow(datpAll[[1]])                     # all components have same number 
                                             # of rows but not cols
ROCMw <- array(NA,dim=c(ntech,nr,nq))
ROCw = array(NA,dim=c(ntech,numfiles,nq))
qAll = matrix(NA,nrow=ntech+1,ncol=nq)
dall = NULL
##############################################################################
## FOR SVM:
ROCMw_SVM <- array(NA,dim=c(ntech,nr,nq-2))
ROCw_SVM = array(NA,dim=c(ntech,numfiles,nq-2))
##############################################################################

for (itech in 1:ntech){                   # itech loop
   datp <- datpAll[[itech]]

## separating dependent and independent variables
 idg = which(colnames(datp)=="GROUP")
 xm = datp[,-idg]
 ym = datp[,idg]                             # single column dataframe


## transform and scale to group 1 (New Randoms)
xmIHS <- IHS(xm,1)
xmn <- scale(subset(xmIHS,ym==1))
mu <- matrix(attr(xmn,"scaled:center"),nrow=1)
sd <- attr(xmn,"scaled:scale")
xms <- scale(xmIHS,center=mu,scale=sd)

## distances from group 1 mean

zero = matrix(rep(0,ncol(xms)),nrow=1)
dmat = as.numeric(dist(xms,zero,method="manhattan"))
dall = c(dall,list(dmat))

## ROC calculations
 qAll[itech,]=quantile(dmat[GROUP==1],1-x)    # group 1 decreasing quantiles
 for (iq in 1:nq) { 
    ROCMw[itech,,iq]=(dmat>=qAll[itech,iq])*1  # q for quantile
}
##############################################################################
## FOR SVM: transform and scale to group 1 (New Randoms) with tail definitions
xm_SVM <- xm
xmN_SVM <- xm_SVM[ym==1,]
inInfo = inputInfo[[itech]]
idv = which(inInfo$OneClass==1)
taildef_SVM <- inInfo$Tail[idv]		# Tail column is the 4th column

for (j in 1:length(taildef_SVM)){
if (taildef_SVM[j]==0){			# Scaling according to tail starts here
med_SVM <- median(xmN_SVM[,j])
xm_SVM[,j] <- xm_SVM[,j]-med_SVM
xmN_SVM[,j] <- xmN_SVM[,j]-med_SVM
leftdatidxN_SVM <- which(xmN_SVM[,j]<0)
rightdatidxN_SVM <- which(xmN_SVM[,j]>0)
leftdatidx_SVM <- which(xm_SVM[,j]<0)
rightdatidx_SVM <- which(xm_SVM[,j]>0)
xm_SVM[leftdatidx_SVM,j] <- -xm_SVM[leftdatidx_SVM,j]/sqrt(sum(xmN_SVM[leftdatidxN_SVM,j]^2)/length(xmN_SVM[leftdatidxN_SVM,j]))
xm_SVM[rightdatidx_SVM,j] <- xm_SVM[rightdatidx_SVM,j]/sqrt(sum(xmN_SVM[rightdatidxN_SVM,j]^2)/length(xmN_SVM[rightdatidxN_SVM,j]))
} else if (taildef_SVM[j]==1) {
max_SVM <- max(xmN_SVM[,j])
xm_SVM[,j] <- (max_SVM-xm_SVM[,j])/sqrt(sum((max_SVM-xmN_SVM[,j])^2)/length(xmN_SVM[,j]))
xm_SVM[which(xm_SVM[,j]<0),j] <- 0
} else if (taildef_SVM[j]==2) {
min_SVM <- min(xmN_SVM[,j])
xm_SVM[,j] <- (xm_SVM[,j]-min_SVM)/sqrt(sum((xmN_SVM[,j]-min_SVM)^2)/length(xmN_SVM[,j]))
xm_SVM[which(xm_SVM[,j]<0),j] <- 0
}
}

xm_SVM <- 2/(1+cosh(xm_SVM))		# SVM scaling finishes
## For SVM: One-class SVM
nu_SVM <- x[2:(nq-1)]
predict_SVM <- matrix(nrow = nrow(xm_SVM), ncol = length(nu_SVM))
#weight_list <- matrix(nrow = length(nu_SVM), ncol = ncol(xm_SVM))
#colnames(weight_list) <- colnames(xm_SVM)
#rho_list <- rep(0,length(nu_SVM))
myTable <- data.frame(xm_SVM[GROUP==1,],GROUP[GROUP==1])
colnames(myTable) <- c(colnames(xm_SVM),"GROUP")				# Use the quantiles as nu parameters
for (j in 1:length(nu_SVM)){				# Train SVM for each nu, j loop
myModel <- svm( GROUP ~., data = myTable, scale = FALSE, type = "one-classification", kernel = "linear", tolerance = 1e-9, nu = nu_SVM[j], cross = 5, probability = FALSE, cachesize = 1000)
#weight <- rep(0,each=ncol(xm_SVM))
#for (k in 1:length(myModel$coefs)){			# calculate the orientation of the classifier hyperplane
#weight <- weight + myModel$SV[k,]*myModel$coefs[k,]
#}								# k loop ends
#weight_list[j,] <- weight
#rho_list[j] <- myModel$rho				
predict_SVM[,j] <- as.numeric(predict(myModel,xm_SVM)==FALSE)
}								# j loop ends

## For SVM: Transfer the SVM predictions with nu series to the rate comparable to L1 output
for (iq in 1:(nq-2)) {
ROCMw_SVM[itech,,iq] <- predict_SVM[,iq]
}
##############################################################################

}                                           # end itech

for (itech in 1:ntech){                     # itech loop
for (igp in 1:numfiles){                    # igp loop
       ROCw[itech,igp,]=apply(ROCMw[itech,GROUP==igp,],2,mean)
}                                           # end igp loop
}                                           # end itech
##############################################################################
## For SVM: 
for (itech in 1:ntech){                     # itech loop
for (igp in 1:numfiles){                    # igp loop
       ROCw_SVM[itech,igp,]=apply(ROCMw_SVM[itech,GROUP==igp,],2,mean)
}                                           # end igp loop
}
##############################################################################    

#### Creating ROCs for the various "1 or more" etc strategies

temp = ROCMw[1,,]
if (ntech > 1) {                            # begin 'if'
  for (itech in 2:ntech){                   # itech loop
  temp = temp+ROCMw[itech,,]
                        }                   # end itech loop
               }                            # end 'if'

ROCMw4pt = array(NA,dim=c(ntech,nr,nq))
for (itech in 1:ntech){                     # itech loop
ROCMw4pt[itech,,] = (temp >= itech)*1
                      }                     # end itech loop

ROCw4pt = array(NA,dim=c(ntech,numfiles,nq))
for (istrat in 1:ntech){
for (igp in 1:numfiles){
       ROCw4pt[istrat,igp,]=apply(ROCMw4pt[istrat,GROUP==igp,],2,mean)
}}

#############################################################################
## For SVM:
temp = ROCMw_SVM[1,,]
if (ntech > 1) {                            # begin 'if'
  for (itech in 2:ntech){                   # itech loop
  temp = temp+ROCMw_SVM[itech,,]
                        }                   # end itech loop
               }                            # end 'if'

ROCMw4pt_SVM = array(NA,dim=c(ntech,nr,nq-2))
for (itech in 1:ntech){                     # itech loop
ROCMw4pt_SVM[itech,,] = (temp >= itech)*1
                      }                     # end itech loop

ROCw4pt_SVM = array(NA,dim=c(ntech,numfiles,nq-2))
for (istrat in 1:ntech){
for (igp in 1:numfiles){
       ROCw4pt_SVM[istrat,igp,]=apply(ROCMw4pt_SVM[istrat,GROUP==igp,],2,mean)
}}

#############################################################################

#### Creating ROCs for combined data from all technologies
xmcIHS <- IHS(xmc,1)
xmcn <- scale(subset(xmcIHS,ym==1))
mu <- matrix(attr(xmcn,"scaled:center"),nrow=1)
sd <- attr(xmcn,"scaled:scale")
xmcs <- scale(xmcIHS,center=mu,scale=sd)

zero = matrix(rep(0,ncol(xmcs)),nrow=1)
dmat = as.numeric(dist(xmcs,zero,method="Manhattan"))
dall = c(dall,list(dmat))

qAll[ntech+1,]=quantile(dmat[GROUP==1],1-x)
ROCMc = array(NA,dim=c(1,nr,nq))
for (i in 1:length(x)) ROCMc[1,,i]=(dmat>=qAll[ntech+1,i])*1
ROCc = array(NA,dim=c(1,numfiles,nq))
for (igp in 1:numfiles){
       ROCc[1,igp,]=apply(ROCMc[1,GROUP==igp,],2,mean)
}
#############################################################################
## For SVM: transform and scale to group 1 (New Randoms) with tail definitions

xm_SVM <- xmc
xmN_SVM <- xm_SVM[ym==1,]
inInfo = inputInfo[[1]]
idv = which(inInfo$OneClass==1)
taildef_SVM <- inInfo$Tail[idv]
for (itech in 2:ntech){ 
inInfo = inputInfo[[itech]]
idv = which(inInfo$OneClass==1)
taildef_SVM <- c(taildef_SVM,inInfo$Tail[idv])
}						# Tail column is the 4th column
cxmc <- cor(xmc_temp)
diag(cxmc) <- 0
idc <- unique(which(cxmc==1,arr.ind=TRUE)[,1])
if (length(idc)>0) {
   idcm <- idc[2:length(idc)]
   taildef_SVM <- taildef_SVM[-idcm] }

for (j in 1:length(taildef_SVM)){
if (taildef_SVM[j]==0){			# Scaling according to tail starts here
med_SVM <- median(xmN_SVM[,j])
xm_SVM[,j] <- xm_SVM[,j]-med_SVM
xmN_SVM[,j] <- xmN_SVM[,j]-med_SVM
leftdatidxN_SVM <- which(xmN_SVM[,j]<0)
rightdatidxN_SVM <- which(xmN_SVM[,j]>0)
leftdatidx_SVM <- which(xm_SVM[,j]<0)
rightdatidx_SVM <- which(xm_SVM[,j]>0)
xm_SVM[leftdatidx_SVM,j] <- -xm_SVM[leftdatidx_SVM,j]/sqrt(sum(xmN_SVM[leftdatidxN_SVM,j]^2)/length(xmN_SVM[leftdatidxN_SVM,j]))
xm_SVM[rightdatidx_SVM,j] <- xm_SVM[rightdatidx_SVM,j]/sqrt(sum(xmN_SVM[rightdatidxN_SVM,j]^2)/length(xmN_SVM[rightdatidxN_SVM,j]))
} else if (taildef_SVM[j]==1) {
max_SVM <- max(xmN_SVM[,j])
xm_SVM[,j] <- (max_SVM-xm_SVM[,j])/sqrt(sum((max_SVM-xmN_SVM[,j])^2)/length(xmN_SVM[,j]))
xm_SVM[which(xm_SVM[,j]<0),j] <- 0
} else if (taildef_SVM[j]==2) {
min_SVM <- min(xmN_SVM[,j])
xm_SVM[,j] <- (xm_SVM[,j]-min_SVM)/sqrt(sum((xmN_SVM[,j]-min_SVM)^2)/length(xmN_SVM[,j]))
xm_SVM[which(xm_SVM[,j]<0),j] <- 0
}
}

xm_SVM <- 2/(1+cosh(xm_SVM))		# SVM scaling finishes
## For SVM: One-class SVM
nu_SVM <- x[2:(nq-1)]
predict_SVM <- matrix(nrow = nrow(xm_SVM), ncol = length(nu_SVM))
myTable <- data.frame(xm_SVM[GROUP==1,],GROUP[GROUP==1])
colnames(myTable) <- c(colnames(xm_SVM),"GROUP")				# Use the quantiles as nu parameters
for (j in 1:length(nu_SVM)){				# Train SVM for each nu, j loop
myModel <- svm( GROUP ~., data = myTable, scale = FALSE, type = "one-classification", kernel = "linear", tolerance = 1e-9, nu = nu_SVM[j], cross = 5, probability = FALSE, cachesize = 1000)				
predict_SVM[,j] <- as.numeric(predict(myModel,xm_SVM)==FALSE)
}								# j loop ends

## For SVM: Transfer the SVM predictions with nu series to the rate comparable to L1 output
ROCMc_SVM = array(NA,dim=c(1,nr,nq-2))
for (iq in 1:(nq-2)) {
ROCMc_SVM[1,,iq] <- predict_SVM[,iq]
}
ROCc_SVM = array(NA,dim=c(1,numfiles,nq-2))
for (igp in 1:numfiles){
       ROCc_SVM[1,igp,]=apply(ROCMc_SVM[1,GROUP==igp,],2,mean)
}
##############################################################################

#### Putting all ROC numbers in an array
ROC = array(NA,dim=c(nstrat,numfiles,nq))
ROC[1:ntech,,]=ROCw
ROC[(ntech+1),,]=ROCc
ROC[(ntech+2):nstrat,,]=ROCw4pt
##############################################################################
## For SVM:
ROC_SVM = array(0,dim=c(nstrat,numfiles,nq-1))
ROC_SVM[1:ntech,,2:(nq-1)]=ROCw_SVM
ROC_SVM[(ntech+1),,2:(nq-1)]=ROCc_SVM
ROC_SVM[(ntech+2):nstrat,,2:(nq-1)]=ROCw4pt_SVM
##############################################################################
## For SVM: calculating percentage to OLD RANDOMS
PStanding_SVM <- matrix(1,nrow=nr,ncol=nstrat)
for (itech in 1:ntech){
for (i in 1:nr){
for (iq in 1:(nq-2))
if (ROCMw_SVM[itech,i,iq]==0) PStanding_SVM[i,itech] <- ROCw_SVM[itech,2,iq]
}
} 

for (i in 1:nr){
for (iq in 1:(nq-2)) {
if (ROCMc_SVM[1,i,iq]==0) PStanding_SVM[i,ntech+1] <- ROCc_SVM[1,2,iq]
}
}

for (istrat in 1:ntech){
for (i in 1:nr){
for (iq in 1:(nq-2))
if (ROCMw4pt_SVM[istrat,i,iq]==0) PStanding_SVM[i,ntech+1+istrat] <- ROCw4pt_SVM[istrat,2,iq]
}
}
PStanding_SVM <- cbind(PStanding_SVM,GROUP)
colnames(PStanding_SVM)<-c("ILL250","Moleculo","PacBio","PlatGen","Combined",
                       "1orMore","2orMore","3orMore","4orMore","SET") 
save(PStanding_SVM,file=paste(resultshome,"PStanding_SVM.RData",sep=""))
 
##############################################################################

#################################################################################
print("Plotting ROC curves. PDF file in results directory");flush.console()

### ROC plots

pdf(file=paste(resultshome,"ROCplots.pdf",sep=""),height=8.5,width=11)
clr = c("gray75","lightpink","lightgreen","blue","red","cyan",
        "brown","magenta","orange")

plot(c(0,100),c(0,100),type="n",xlab="",ylab="",cex.main=1,main="Legend Colors")
legend(7.5,97.5,legend=strats,fill=clr,cex=2)

for (igp in 1:numfiles){
for (istrat in 1:nstrat){
uuu = 100*ROC[istrat,2,]
vvv = 100*ROC[istrat,igp,]
if (istrat == 1) plot(uuu,vvv,type="l",col=clr[istrat],ylim=c(0,100),
xlab="False Positive (%) Using Old Randoms",xlim=c(0,20),cex.main=1,
ylab="True Positive Percent",main=paste("ROC Curves for ",SET[igp])) else
lines(uuu,vvv,col=clr[istrat]) 
abline(0,1,col="gray40",lty=2)
abline(h=c(20,40,60,80,90,95),col="gray40",lty=3)
axis(side=2,at=90,labels=90,cex.axis=0.8)
abline(v=c(5,10,20),lty=3,col="gray40")
}
}
dev.off(dev.cur())
#################################################################################

#### Plots of distances

pdf(file=paste(resultshome,"PlotOfAllDistances.pdf",sep=""),width=11,height=8.5)
par(mfrow=c(2,2))
for (igp in 1:numfiles){
for (itech in 1:ntech){
idx = which(GROUP %in% c(1,igp))
xc = dall[[itech]][idx]
gp = GROUP[idx]
colr=c(rep(8,4000),rep(2,length(gp)-4000))
cexr=c(rep(0.5,4000),rep(1,length(gp)-4000))
plot(xc,col=colr,cex=cexr,xlab="",ylab="Distance",
    main=paste("Distances for ",tech[itech],"\n ",SET[1]," (Gray) and\n",
    SET[igp]," (Red)"),cex.main=0.75)
#readline()
}}
dev.off(dev.cur())


#################### Find the Percentile Standings of Each Sequence ##########
print("Calculating Percentile Scores. This takes a while");flush.console()

PStanding <- matrix(NA,nrow=nr,ncol=nstrat)
gpsize = table(GROUP)
for (itech in 1:(ntech+1)){
print(itech);flush.console()
dmat = dall[[itech]]
for (ir in 1:nr){
PStanding[ir,itech]= sum(dmat[ir] >= dmat[GROUP==2])/gpsize[2]
}}

### now percentiles for 1orMore, etc.

Dall <- cbind(dall[[1]],dall[[2]],dall[[3]],dall[[4]])
Dall0 <- Dall[GROUP==2,]
for (ir in 1:nr){
if (ir%%500 == 0) print("Working Very Hard");flush.console()
temp1 <- apply(t(Dall[ir,]>=t(Dall0)),1,sum)
for (itech in 1:4){
istrat = ntech+1+itech
PStanding[ir,istrat] <- mean(temp1>=itech)
}}

##Calibrate the 1ormore, 2ormore, columns relative to OLD RANDOMS                        }
temp <- PStanding
for (icol in (ntech+2):nstrat){
print("Working");flush.console()
for (ir in 1:nr){
temp[ir,icol]=mean(PStanding[ir,icol]>=PStanding[GROUP==2,icol])
}}

temp <- cbind(temp,SET=GROUP)
PStanding <- as.data.frame(round(temp,4))
colnames(PStanding)<-c("ILL250","Moleculo","PacBio","PlatGen","Combined",
                       "1orMore","2orMore","3orMore","4orMore","SET")

save(PStanding,file=paste(resultshome,"PStanding.RData",sep=""))

### Write the info in excel files (csv files)

for (itech in 1:ntech){                  # itech loop
dname = gsub("FILE",tech[itech],fsave)
load(dname)
datx <- cbind(dat,PStanding)

write.table(datx,
    file=paste(resultshome,tech[itech],"_ROC.csv",sep=""),
    col.names=TRUE,row.names=FALSE,sep=",")
## For SVM:
datx_SVM <- cbind(dat,PStanding_SVM)

write.table(datx_SVM,
    file=paste(resultshome,tech[itech],"_ROC_SVM.csv",sep=""),
    col.names=TRUE,row.names=FALSE,sep=",")
}

# Write Pstanding to excel file
write.table(PStanding,
    file=paste(resultshome,"PercentileScores.csv",sep=""),
    col.names=TRUE,row.names=FALSE,sep=",")

write.table(PStanding_SVM,
    file=paste(resultshome,"PercentileScores_SVM.csv",sep=""),
    col.names=TRUE,row.names=FALSE,sep=",")

##################################################################
## Compare L1 and SVM results
load(file=paste(resultshome,"PStanding.RData",sep=""))
load(file=paste(resultshome,"PStanding_SVM.RData",sep="")) 	# load for debug

compare_L1_SVM <- array(NA,dim = c(length(x)-2,nr,nstrat)) 	# compare L1 and SVM for each SVM
compare_summary <- array(NA,dim = c(nstrat,length(x)-2,5)) 	# summarize to get the confusion matrix

for (i in 2:(length(x)-1)){
decision_L1 <- 1*(PStanding[,1:nstrat]>=1-x[i])
decision_SVM <- 1*(PStanding_SVM[,1:nstrat]>=1-x[i])
compare_L1_SVM[i-1,,] <- 1.5*(decision_L1^2+decision_SVM^2)+0.5*(decision_SVM-decision_L1)
}
for (istrat in 1:nstrat){
compare_summary[istrat,,1] <- 1-x[2:(length(x)-1)]
for (i in 2:(length(x)-1)){
compare_summary[istrat,i-1,2] <- sum(1*compare_L1_SVM[i-1,,istrat]==0) 	# 0: L1- SVM-
compare_summary[istrat,i-1,3] <- sum(1*compare_L1_SVM[i-1,,istrat]==1)	# 1: L1+ SVM-
compare_summary[istrat,i-1,4] <- sum(1*compare_L1_SVM[i-1,,istrat]==2)	# 2: L1- SVM+
compare_summary[istrat,i-1,5] <- sum(1*compare_L1_SVM[i-1,,istrat]==3)	# 3: L1+ SVM+
}
}

save(compare_L1_SVM,compare_summary,file=paste(resultshome,"Compare_L1_SVM.RData",sep="")) 

print("DONE")
############################## THE END ####################################