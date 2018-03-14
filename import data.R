#Script to import data and clinical response 


library(data.table)

#Import raw file

blind.code<-read.csv("/data/metaderm/blind decoded.csv")

easi.area<-read.csv("/data/metaderm/easi area.csv")
#####
#easi.severity
data.tmp<-read.csv("/data/metaderm/easi severity.csv")
storage<-data.tmp[,1:6]
storage<-cbind(3,storage)
visit<-c(3,4,5,6,7,8)
s4<-cbind(4,data.tmp[,c(1,2,7:10)])
s5<-cbind(5,data.tmp[,c(1,2,11:14)])
s6<-cbind(6,data.tmp[,c(1,2,15:18)]) 
s7<-cbind(7,data.tmp[,c(1,2,19:22)])
s8<-cbind(8,data.tmp[,c(1,2,23:26)])
c.names<-c("visit","PID","Site","easi.Erythema","easi.Papulation","easi.Excoriation","easi.Lichenification")
colnames(storage)<-colnames(s4)<-colnames(s5)<-colnames(s6)<-colnames(s7)<-colnames(s8)<-c.names
easi.severity<-rbind(storage,s4,s5,s6,s7,s8)
#####
#scorad.area
scorad.area<-read.csv("/data/metaderm/scorad area.csv")
#####
#scorad.partC
data.tmp<-read.csv("/data/metaderm/scorad part c.csv")
s3<-cbind(3,data.tmp[,1:3])
s4<-cbind(4,data.tmp[,c(1,4,5)])
s5<-cbind(5,data.tmp[,c(1,6,7)])
s6<-cbind(6,data.tmp[,c(1,8,9)])
s7<-cbind(7,data.tmp[,c(1,10,11)])
s8<-cbind(8,data.tmp[,c(1,12,13)])

col.names<-c("visit","PID","P","S")
colnames(s3)<-colnames(s4)<-colnames(s5)<-colnames(s6)<-colnames(s7)<-colnames(s8)<-col.names
scorad.partC<-rbind(s3,s4,s5,s6,s7,s8)
#####
#scorad severity
data.tmp<-read.csv("/data/metaderm/scorad severity.csv")
s3<-cbind(3,data.tmp[,c(1:2,3:8)])
s4<-cbind(4,data.tmp[,c(1:2,9:14)])
s5<-cbind(5,data.tmp[,c(1:2,15:20)])
s6<-cbind(6,data.tmp[,c(1:2,21:26)])
s7<-cbind(7,data.tmp[,c(1:2,27:32)])
s8<-cbind(8,data.tmp[,c(1:2,33:38)])

col.names<-c("visit","PID","area", "Erythema","Oedema","Crust","Excoriation","Lichenification","Dryness")
colnames(s3)<-colnames(s4)<-colnames(s5)<-colnames(s6)<-colnames(s7)<-colnames(s8)<-col.names
scorad.severity<-rbind(s3,s4,s5,s6,s7,s8)

tape.strip<-read.csv("/data/metaderm/tape data.csv")

mastr.clin<-read.csv("/data/metaderm/master clinical.csv")

#BL - v3
#W12 - v7 last treatment TP
#W6 - v8 Follow-up TP

#Scorad score is A/5 + 7B/2 + C

#A - scorad.area
#B - scorad.severity
#C- Scorad.partC



# Computing Scorad Scores.  Also store data in long mode.

pids<-unique(scorad.area[,"PID"])
visits<-c(3,4,5,6,7,8)
time<-c(0,2,4,8,12,16)

storage.mat<-matrix(data=NA,nrow=length(pids)*length(visits),ncol = 9)
colnames(storage.mat)<-c("PID","visit", "time", "A","B","C","scorad","easi","pga")
#PID Visit Time Scorad A B C Final
counter<-1

for(i in pids){
  
  for(j in visits){
    data.tmp<-scorad.severity[scorad.severity$PID==i & scorad.severity$visit==j,-c(1:3)]
    storage.mat[counter,1]<-i
    storage.mat[counter,2]<-j
    storage.mat[counter,5]<-sum(data.tmp)
    
    storage.mat[counter,6]<-sum(scorad.partC[scorad.partC$PID==i & scorad.partC$visit==j,-c(1,2)])
    
    counter<-counter+1
  }
  
  
}
tmp<-NULL
for(i in pids){
  data.tmp<-scorad.area[scorad.area$PID==i,-c(1:2)]
  values<-colSums(data.tmp)
  out<-cbind(i,cbind(c(3,4,5,6,7,8),values))
  tmp<-rbind(tmp,out)
}

storage.mat[,4]<-as.numeric(tmp[,3])

tmp<-NULL
for(i in pids){
  data.tmp<-mastr.clin[mastr.clin$PID==i,-c(1)]
  scorad<-cbind(i,cbind(c(3,4,5,6,7,8),as.numeric(data.tmp[1:6])))
  easi<-cbind(i,cbind(c(3,4,5,6,7,8),as.numeric(data.tmp[7:12])))
  pga<-cbind(i,cbind(c(3,4,5,6,7,8),as.numeric(data.tmp[13:18])))
  
  tmp<-rbind(tmp,cbind(scorad,easi,pga))
}

storage.mat[,7]<-as.numeric(tmp[,3])
storage.mat[,8]<-as.numeric(tmp[,6])
storage.mat[,9]<-as.numeric(tmp[,9])

smat<-as.data.frame(storage.mat)
setDT(smat)

save(smat,tape.strip,file="/data/metaderm/metaderm data.Rdata",compress = "xz")
