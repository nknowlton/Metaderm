 #plotting and inital analysis


load(file="/data/metaderm/metaderm data.Rdata")

library(ggplot2)
library(data.table)

setDT(tape.strip)
setDT(blind.code)
colnames(tape.strip)<-c("PID","Protein","trt","site","time","site.trt","data")

ts.dat<-merge(tape.strip,blind.code)
visit<-as.character(ts.dat[,time])
visit[visit=="BL"]<-3
visit[visit=="W2"]<-4
visit[visit=="W4"]<-5
visit[visit=="W8"]<-6
visit[visit=="W12"]<-7
visit[visit=="W16"]<-8
visit<-as.numeric(visit)
ts.dat[,visit:=visit]
ts.dat[,data2:=log2(2^data+1)]

g<-ggplot(smat,aes(visit,scorad))+geom_point()

ggplot(ts.dat,aes(visit,data))+geom_point()

mean.plot.data<-ts.dat[,mean(data2,na.rm = TRUE),by=list(visit,Protein,trt,site)]

ggplot(mean.plot.data,aes(visit,V1,colour=trt,fill=site))+geom_point()+facet_wrap(~Protein)

nolesionremoved<-ts.dat[ts.dat$site=="L"]
mean.plot.data.notrt<-nolesionremoved[,mean(data2,na.rm = TRUE),by=list(visit,Protein,trt)]
ggplot(mean.plot.data.notrt,aes(visit,V1,colour=trt))+geom_point()+facet_wrap(~Protein)
#Run order check
dim(ts.dat)
plot(1:39168,ts.dat[,data2])
cor.test(1:39168,ts.dat[,data2])



#####
# Scorad Plots

score.dat<-merge(smat,blind.code)

score.mean.data<-score.dat[,mean(scorad,na.rm = TRUE),by=list(visit,Class)]
ggplot(score.mean.data,aes(visit,V1,colour=Class))+geom_point()

# Create a ratio of clean to infected skin

who<-unique(ts.dat$PID)
what<-unique(ts.dat$Protein)
time<-unique(ts.dat$time)
junk.mat<-NULL
counter=1

for (i in who) {
  print(i)
  for (j in what) {
    #print(j)
    for (k in time){
      tmp.dat<-ts.dat[ts.dat$PID==i & ts.dat$Protein==j & ts.dat$time==k,]
     # ratio.mat[counter,]<-as.data.frame(tmp.dat[1,])
      tmp.dat[1,7]<-tmp.dat[tmp.dat$site=='L',data2]/tmp.dat[tmp.dat$site=='N',data2]
      tmp.dat<-tmp.dat[1,]
      tmp.dat[1,4]="ratio"
      tmp.dat[1,6]="ratio"
       
      junk.mat<-rbind(junk.mat,tmp.dat)
      
      counter=counter+1
    } #end of k
  } # end of j
} # end of i

colnames(junk.mat)<-c("PID","Protein","trt", "site", "time", "site.trt","ratio.data","Class","visit", "data2")
ratio.mat<-junk.mat
rm(junk.mat)


save(ratio.mat,file="/data/metaderm/metaderm ratio data.Rdata",compress = "xz")

mean.plot.data2<-ratio.mat[,mean(ratio.data,na.rm = TRUE),by=list(visit,Protein,trt)]
ggplot(mean.plot.data2,aes(visit,V1,colour=trt))+geom_point()+facet_wrap(~Protein)
######
# Inital analysis

# Normalised ratio.data across proteins at week 8 only

w12.dat<-ratio.mat[ratio.mat$time=="W12"]

prots<-unique(w12.dat$Protein)
storage.mat1<-data.frame(matrix(data=NA,nrow=length(prots),ncol=6))
colnames(storage.mat1)<-c("Protein","week","Drug Mean","Veh Mean","Raw p","BH p")

storage.mat1[,"week"]<-8

counter=1
for(i in prots) {
  print(i)
  storage.mat1[counter,1]<-i
  tmp.dat<-w12.dat[w12.dat$Protein==i,]
  out<-t.test(ratio.data~Class,data=tmp.dat)
  storage.mat1[counter,5]<-out$p.value
  storage.mat1[counter,3]<-out$estimate[1]
  storage.mat1[counter,4]<-out$estimate[2]
  
  counter=counter+1
  
}

#BH values

storage.mat1[,6]<-p.adjust(storage.mat1$`Raw p`,method="BH")
write.csv(storage.mat1,file="/data/metaderm/W12 Normalised Protein Ttest.csv")


######
# Repeat for week 1
# Rationale is to compare baseline to week 8. Proteins should be non-sig at baseline and sig at w8.
# Normalised ratio.data across proteins at week 1 only

w1.dat<-ratio.mat[ratio.mat$time=="BL"]

prots<-unique(w1.dat$Protein)
storage.mat2<-data.frame(matrix(data=NA,nrow=length(prots),ncol=6))
colnames(storage.mat2)<-c("Protein","week","Drug Mean","Veh Mean","Raw p","BH p")

storage.mat2[,"week"]<-0

counter=1
for(i in prots) {
  print(i)
  storage.mat2[counter,1]<-i
  tmp.dat<-w1.dat[w1.dat$Protein==i,]
  out<-t.test(ratio.data~Class,data=tmp.dat)
  storage.mat2[counter,5]<-out$p.value
  storage.mat2[counter,3]<-out$estimate[1]
  storage.mat2[counter,4]<-out$estimate[2]
  
  counter=counter+1
  
}

#BH values

storage.mat2[,6]<-p.adjust(storage.mat2$`Raw p`,method="BH")
write.csv(storage.mat2,file="/data/metaderm/BL Normalised Protein Ttest.csv")

######
# Time 0 minus T8 by trt

diff.dat<-w1.dat
diff.dat[,7]<-w1.dat[,7]-w12.dat[,7]

w1.dat<-ratio.mat[ratio.mat$time=="BL"]

prots<-unique(w1.dat$Protein)
storage.mat3<-data.frame(matrix(data=NA,nrow=length(prots),ncol=6))
colnames(storage.mat3)<-c("Protein","week","Drug Mean","Veh Mean","Raw p","BH p")

storage.mat3[,"week"]<-"BL-W8"

counter=1
for(i in prots) {
  print(i)
  storage.mat3[counter,1]<-i
  tmp.dat<-diff.dat[diff.dat$Protein==i,]
  out<-t.test(ratio.data~Class,data=tmp.dat)
  storage.mat3[counter,5]<-out$p.value
  storage.mat3[counter,3]<-out$estimate[1]
  storage.mat3[counter,4]<-out$estimate[2]
  
  counter=counter+1
  
}

#BH values

storage.mat3[,6]<-p.adjust(storage.mat3$`Raw p`,method="BH")
write.csv(storage.mat3,file="/data/metaderm/BLminusW12 Normalised Protein Ttest.csv")

#######
# Scorad by Protein associations

prots<-unique(ratio.mat$Protein)
storage.mat<-data.frame(matrix(data=NA,nrow=length(prots),ncol=4))
colnames(storage.mat)<-c("Protein","Cor","Raw p","BH p")

counter=1
for (i in prots){
  tmp.mat<-cbind(score.dat,ratio.mat[ratio.mat$Protein==i,])
  
  a<-cor.test(tmp.mat$scorad,tmp.mat$ratio.data)
  est<-a$estimate
  p<-a$p.value
 
  #if(p<0.05) print(cbind(i,est))
  storage.mat[counter,1]<-i
  storage.mat[counter,2]<-est
  storage.mat[counter,3]<-p
  
   counter=counter+1
}

storage.mat[,4]<-p.adjust(storage.mat$`Raw p`,method="BH")
write.csv(storage.mat,file="/data/metaderm/Cor Normalised Protein SCORAD.csv")
