
### Unsupervised Clustering and MDS Analysis script. Prepared for the svclassify Paper
### April 6 2015.  Hari Iyer, Statistical Engineering Division, ITL, NIST


### Step 1: Data Input and Data Preparation for L1 - One Class Classification

### Modify the path for home to suit your situation

home <- "svclassify/"                          # Base directory for the project
tech <- c("ILL250","Moleculo","PacBio","PlatGen") #root names for technologies

### Make sure that this base directory has the following subdirectories

### Code:         /svclassify/code
### Results:      /svclassify/results
### Data:         /svclassify/data
###       Under the data directory you need to have one subdirectory for each technology.
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
idv = which(inInfo[,3]==1)
keepZ = inInfo[idv,1]
temp <- dat[,c("GROUP",keepZ)]          

#### Remove constant columns
 print(itech);flush.console()
 idsd <- which(apply(temp,2,sd,na.rm=TRUE)==0)
 if (length(idsd)>0) temp <- temp[,-idsd]

#### Remove rows with missing values in the selected variables
 idmiss <- unique(which(is.na(temp)==TRUE,arr.ind=TRUE)[,1])
 if (length(idmiss) > 0) datp <- temp[-idmiss,] else datp <- temp
 rm(temp)

datpAll <- c(datpAll,list(datp))
}                                        # end loop itech

save(datpAll,file=paste(datahome,"datpAllForUnsupervised.RData",sep=""))

###############################################################################

#### datpAll reloaded. For testing and debugging one may start here if datpAll is
   # already created and stored.
load(file=paste(datahome,"datpAllForUnsupervised.RData",sep=""))   

#### Form the combined data from all technologies
print("Combining Data from All Technologies");flush.console()
  
GROUP <- datpAll[[1]]$GROUP
datp1 <- datpAll[[1]][GROUP<=numfiles,-1]
datp2 <- datpAll[[2]][GROUP<=numfiles,-1]
datp3 <- datpAll[[3]][GROUP<=numfiles,-1]
datp4 <- datpAll[[4]][GROUP<=numfiles,-1]
colnames(datp1) <- paste("1_",colnames(datp1),sep="")
colnames(datp2) <- paste("2_",colnames(datp2),sep="")
colnames(datp3) <- paste("3_",colnames(datp3),sep="")
colnames(datp4) <- paste("4_",colnames(datp4),sep="")
GROUP <- GROUP[GROUP<=numfiles]

datp <- as.data.frame(cbind(GROUP,datp1,datp2,datp3,datp4))
corm <- cor(datp)
idcor <- which(abs(corm-diag(diag(corm))) > 0.9999,arr.ind=TRUE)
#idcor <- which(abs(corm-diag(diag(corm))) == 1,arr.ind=TRUE)
idx <- NULL
for (jj in 1:nrow(idcor)){
if (idcor[jj,1]>idcor[jj,2]) idx <- c(idx,idcor[jj,1])
                         }
datp <- datp[,-idx]

### datp now has the right metrics from all the technologies with no duplicate columns

## separating dependent and independent variables
idg = which(colnames(datp)=="GROUP")
xmall = datp[-idg]
ymall = datp[,idg]     # single column dataframe

########################## MDS ANALYSIS ############################################

xmallIHS <- IHS(xmall,1)
#xmIHS <- xmallIHS[ymall!=numfiles,]

xms <- scale(xmallIHS)   # scaled to the entire data set (not just randoms)
mu <- matrix(attr(xms,"scaled:center"),nrow=1)
sd <- attr(xms,"scaled:scale")
#xmsnew <- scale(xmnewIHS,center=mu,scale=sd)

print("Calculating Distance Matrix - Takes a while");flush.console()
dmat = dist(xms,method="Manhattan")
#save(dmat,file="dmatForUnsupervised.RData")

#load("dmatForUnsupervised.RData")

aa = flashClust(dmat,method="ward")
hcd = as.dendrogram(aa)

postscript(file=paste(resultshome,"Dendrogram.eps",sep=""),width=11,height=8.5)
plot(hcd,leaflab="none",ylab="Cluster Dissimilarity Index",
     main = "Hierarchical Clustering (Ward Method)
L1 Distance - IHS Transformed and Standardized")

h=10000
abline(h=h,lty=2,col=2)
mtext(side=1,at=c(10,10),line=2,adj=0,
  text="Horizontal Dotted Red Line shows the 8 cluster split",col=4)
dev.off(dev.cur())

#bb=cut(hcd,h=h)
#plot(bb$upper,center=TRUE,main="Clusters Based on L1 Distances (All Tech)")
#mtext(side=1,at=c(20,10),line=-10,adj=0,
#  text="Only the 8 major branches are shown.\n Further subclusters are not shown.",col=4)

groups1 <- slice(aa, h = h)[order(aa$order)]
(lg = max(groups1))
(tab = table(groups1,GROUP))
#save(groups1,file="h10000_7files_8Clusters_Feb082015.RData")
A1 = names(groups1)[groups1 == 1]
A2 = names(groups1)[groups1 == 2]
A3 = names(groups1)[groups1 == 3]
A4 = names(groups1)[groups1 == 4]
A5 = names(groups1)[groups1 == 5]
A6 = names(groups1)[groups1 == 6]
A7 = names(groups1)[groups1 == 7]
A8 = names(groups1)[groups1 == 8]

################################ Distances for the new 39 GOld from clusters #

df1 = data.frame(matrix(tab,nrow=lg))
colnames(df1)=SET
rownames(df1)=paste("Cluster",1:lg)
#colnames(df1)=c("NEW RANDOMS","OLD RANDOMS","RM_LINE","RM_LTR","RM_SINE",
#                "PERSONALIS_GOLD","1KG_PILOT")
print(df1)
write.table(df1,file=paste(resultshome,"ClusterTable.csv",sep=""),sep=",")


#### MDS analysis

print("MDS Calculations. Takes a LONG time");
print("Take a break and Come back");flush.console()

mds = cmdscale(dmat,10)
yl = min(mds[,2])
yu = max(mds[,2])
xl = min(mds[,1])
xu = max(mds[,1])
zl = min(mds[,3])
zu = max(mds[,3])

clrs = c("black","red","green","blue","cyan","magenta","yellow","gray")
         #"aquamarine","bisque","deeppink1","yellow","brown",
         #"cyan","chocolate","blue")
         #"red",,"darkgoldenrod1","firebrick1","chartreuse",
         #"blueviolet","darkgreen","darkseagreen1")
         #,"indianred1","khaki"

postscript(file=paste(resultshome,"MDS1vsMDS2.ps",sep=""),width=11,height=8.5)
plot(mds,col=clrs[groups1],#main="MDS-1 vs MDS-2 plot for 8 Clusters",
     xlab="MDS-1",ylab="MDS-2",ylim=c(yl,yu),xlim=c(xl,xu))
legend(-150,yu,legend=c("Cluster 1","Cluster 2","Cluster 3","Cluster 4",
        "Cluster 5","Cluster 6","Cluster 7","Cluster 8"),
        fill=clrs,cex=0.7)
dev.off(dev.cur())

postscript(file=paste(resultshome,"MDS1vsMDS3.ps",sep=""),width=11,height=8.5)
plot(mds[,c(1,3)],col=clrs[groups1],xlab="MDS-1",ylab="MDS-3",
     ylim=c(zl,zu),xlim=c(xl,xu))
legend(-150,zu,legend=c("Cluster 1","Cluster 2","Cluster 3","Cluster 4",
        "Cluster 5","Cluster 6","Cluster 7","Cluster 8"),
        fill=clrs,cex=0.7)
dev.off(dev.cur())





