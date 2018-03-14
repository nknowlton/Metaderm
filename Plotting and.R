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



#Run order check
dim(ts.dat)
plot(1:39168,ts.dat[,data2])
cor.test(1:39168,ts.dat[,data2])



#####
# Scorad Plots

score.dat<-merge(smat,blind.code)

score.mean.data<-score.dat[,mean(scorad,na.rm = TRUE),by=list(visit,Class)]
ggplot(score.mean.data,aes(visit,V1,colour=Class))+geom_point()
