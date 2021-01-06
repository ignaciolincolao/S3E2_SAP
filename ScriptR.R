library(rayshader)
library(ggplot2)
library(tidyverse)
library(reshape2)
library(MBA)
library(dplyr)

###########################
##### Cargar datos
###########################

sequential_old.data_111 = read.csv("~/Proyectos/S3E2_SAP/sequential_recalentamiento_V2/cmake-build-debug/info_save_test_1_1_1.txt", header = FALSE, sep=",")
sequential_old.data_111 = as.data.frame(sequential_old.data_111)
colnames(sequential_old.data_111)[1]<-"time"
colnames(sequential_old.data_111)[2]<-"Z"
colnames(sequential_old.data_111)[3]<-"meandist"
colnames(sequential_old.data_111)[4]<-"segregation"
colnames(sequential_old.data_111)[5]<-"costquota"
colnames(sequential_old.data_111)[6]<-"cycle"
colnames(sequential_old.data_111)[13]<-"seed"
as.double(unlist(sequential_old.data_111["time"]))
as.double(unlist(sequential_old.data_111["Z"]))
as.double(unlist(sequential_old.data_111["meandist"]))
as.double(unlist(sequential_old.data_111["segregation"]))
as.double(unlist(sequential_old.data_111["costquota"]))
as.numeric(unlist(sequential_old.data_111["cycle"]))
as.numeric(unlist(sequential_old.data_111["seed"]))



sequential_old.data_153025 = read.csv("~/Proyectos/S3E2_SAP/sequential_recalentamiento_V2/cmake-build-debug/info_save_test_15_30_25.txt", header = FALSE, sep=",")
sequential_old.data_153025 = as.data.frame(sequential_old.data_153025)
colnames(sequential_old.data_153025)[1]<-"time"
colnames(sequential_old.data_153025)[2]<-"Z"
colnames(sequential_old.data_153025)[3]<-"meandist"
colnames(sequential_old.data_153025)[4]<-"segregation"
colnames(sequential_old.data_153025)[5]<-"costquota"
colnames(sequential_old.data_153025)[6]<-"cycle"
colnames(sequential_old.data_153025)[13]<-"seed"
as.double(unlist(sequential_old.data_153025["time"]))
as.double(unlist(sequential_old.data_153025["Z"]))
as.double(unlist(sequential_old.data_153025["meandist"]))
as.double(unlist(sequential_old.data_153025["segregation"]))
as.double(unlist(sequential_old.data_153025["costquota"]))
as.numeric(unlist(sequential_old.data_153025["cycle"]))
as.numeric(unlist(sequential_old.data_153025["seed"]))


sequential.data_111 = read.csv("~/Proyectos/S3E2_SAP/sequential_recalentamiento/cmake-build-debug/info_save_test_1_1_1.txt", header = FALSE, sep=",")
sequential.data_111 = as.data.frame(sequential.data_111)
colnames(sequential.data_111)[1]<-"time"
colnames(sequential.data_111)[2]<-"Z"
colnames(sequential.data_111)[3]<-"meandist"
colnames(sequential.data_111)[4]<-"segregation"
colnames(sequential.data_111)[5]<-"costquota"
colnames(sequential.data_111)[6]<-"cycle"
colnames(sequential.data_111)[13]<-"seed"
as.double(unlist(sequential.data_111["time"]))
as.double(unlist(sequential.data_111["Z"]))
as.double(unlist(sequential.data_111["meandist"]))
as.double(unlist(sequential.data_111["segregation"]))
as.double(unlist(sequential.data_111["costquota"]))
as.numeric(unlist(sequential.data_111["cycle"]))
as.numeric(unlist(sequential.data_111["seed"]))

sequential.data_153025 = read.csv("~/Proyectos/S3E2_SAP/sequential_recalentamiento/cmake-build-debug/info_save_test_15_30_25.txt", header = FALSE, sep=",")
sequential.data_153025 = as.data.frame(sequential.data_153025)
colnames(sequential.data_153025)[1]<-"time"
colnames(sequential.data_153025)[2]<-"Z"
colnames(sequential.data_153025)[3]<-"meandist"
colnames(sequential.data_153025)[4]<-"segregation"
colnames(sequential.data_153025)[5]<-"costquota"
colnames(sequential.data_153025)[6]<-"cycle"
colnames(sequential.data_153025)[13]<-"seed"
as.double(unlist(sequential.data_153025["time"]))
as.double(unlist(sequential.data_153025["Z"]))
as.double(unlist(sequential.data_153025["meandist"]))
as.double(unlist(sequential.data_153025["segregation"]))
as.double(unlist(sequential.data_153025["costquota"]))
as.numeric(unlist(sequential.data_153025["cycle"]))
as.numeric(unlist(sequential.data_153025["seed"]))

sequential.data_153025_v3 = read.csv("~/Proyectos/S3E2_SAP/sequential_recalentamiento/cmake-build-debug/info_save_test_v3_15_30_25.txt", header = FALSE, sep=",")
sequential.data_153025_v3 = as.data.frame(sequential.data_153025_v3)
colnames(sequential.data_153025_v3)[1]<-"time"
colnames(sequential.data_153025_v3)[2]<-"Z"
colnames(sequential.data_153025_v3)[3]<-"meandist"
colnames(sequential.data_153025_v3)[4]<-"segregation"
colnames(sequential.data_153025_v3)[5]<-"costquota"
colnames(sequential.data_153025_v3)[6]<-"cycle"
colnames(sequential.data_153025_v3)[13]<-"seed"
as.double(unlist(sequential.data_153025_v3["time"]))
as.double(unlist(sequential.data_153025_v3["Z"]))
as.double(unlist(sequential.data_153025_v3["meandist"]))
as.double(unlist(sequential.data_153025_v3["segregation"]))
as.double(unlist(sequential.data_153025_v3["costquota"]))
as.numeric(unlist(sequential.data_153025_v3["cycle"]))
as.numeric(unlist(sequential.data_153025_v3["seed"]))

#########################
##### rutas
#########################

sequential.rute_153025 = "/home/ignacio/Proyectos/S3E2_SAP/sequential_recalentamiento/cmake-build-debug/save_example_15_30_25/"
sequential.rute_153025_v3 = "/home/ignacio/Proyectos/S3E2_SAP/sequential_recalentamiento/cmake-build-debug/save_example_15_30_25_v3/"
sequential.rute_111 = "/home/ignacio/Proyectos/S3E2_SAP/sequential_recalentamiento/cmake-build-debug/save_example_1_1_1/"
sequential_old.rute_153025 = "/home/ignacio/Proyectos/S3E2_SAP/sequential_recalentamiento_V2/cmake-build-debug/save_example_15_30_25/"
sequential_old.rute_111 = "/home/ignacio/Proyectos/S3E2_SAP/sequential_recalentamiento_V2/cmake-build-debug/save_example_1_1_1/"


###########################
##### Agrupa datos
###########################

## Desviación standar
sequential_old.sd_111 <- sapply(sequential_old.data_111[,1:8],sd)
sequential_old.sd_153025 <- sapply(sequential_old.data_153025[,1:8],sd)
sequential.sd_111 <- sapply(sequential.data_111[,1:8],sd)
sequential.sd_153025 <- sapply(sequential.data_153025[,1:8],sd)
sequential.sd_153025_v3 <- sapply(sequential.data_153025_v3[,1:8],sd)
## Promedio
sequential_old.mean_111 <- sapply(sequential_old.data_111[,1:8], mean)
sequential_old.mean_153025 <- sapply(sequential_old.data_153025[,1:8], mean)
sequential.mean_111 <- sapply(sequential.data_111[,1:8], mean)
sequential.mean_153025 <- sapply(sequential.data_153025[,1:8], mean)
sequential.mean_153025_v3 <- sapply(sequential.data_153025_v3[,1:8],mean)
## mediana
sequential_old.median_111 <- sapply(sequential_old.data_111[,1:8], median)
sequential_old.median_153025 <- sapply(sequential_old.data_153025[,1:8], median)
sequential.median_111 <- sapply(sequential.data_111[,1:8], median)
sequential.median_153025 <- sapply(sequential.data_153025[,1:8], median)
sequential.median_153025_v3 <- sapply(sequential.data_153025_v3[,1:8],median)

########################
###### Graficar ejecución
#######################

number = 29
ruta = paste(sequential.rute_153025,number,"-info-graficos.txt", sep="")

mydata = as.matrix(read.csv(ruta), header = FALSE, sep=",")
colnames(mydata)[1] <- "count"
colnames(mydata)[2] <- "Distancia"
colnames(mydata)[3] <- "Segregacion"
colnames(mydata)[4] <- "Costo_Cupo"
colnames(mydata)[5] <- "Costo_Solucion"
par(mfrow=c(2,2))
#par(mfrow=c(1,1))
plot(mydata[,"count"],mydata[,"Distancia"], ylab=bquote(bar(d)), xlab="iterations", type="l",col="blue")
plot(mydata[,"count"],mydata[,"Segregacion"], ylab="D", xlab="iterations" , type="l",col="blue")
plot(mydata[,"count"],mydata[,"Costo_Cupo"], ylab="C", xlab="iterations", type="l",col="blue")
plot(mydata[,"count"],mydata[,"Costo_Solucion"], ylab="Z", xlab="iterations", type="l",col="blue")



#######
### APARTE
######

mydata = as.matrix(read.csv("/home/ignacio/Proyectos/S3E2_SAP/sequential_recalentamiento/cmake-build-debug/save/2020-09-29 T:17-33-info-graficos.txt"), header = FALSE, sep=",")
colnames(mydata)[1] <- "count"
colnames(mydata)[2] <- "Distancia"
colnames(mydata)[3] <- "Segregacion"
colnames(mydata)[4] <- "Costo_Cupo"
colnames(mydata)[5] <- "Costo_Solucion"
par(mfrow=c(2,2))
#par(mfrow=c(1,1))
plot(mydata[,"count"],mydata[,"Distancia"], ylab=bquote(bar(d)), xlab="iterations", type="l",col="blue")
plot(mydata[,"count"],mydata[,"Segregacion"], ylab="D", xlab="iterations" , type="l",col="blue")
plot(mydata[,"count"],mydata[,"Costo_Cupo"], ylab="C", xlab="iterations", type="l",col="blue")
plot(mydata[,"count"],mydata[,"Costo_Solucion"], ylab="Z", xlab="iterations", type="l",col="blue")


##########
### Geom
##########

mydata2 = data.frame(mydata) 
ggplot(data = mydata2, aes(x = count, y = Costo_Solucion))+
  labs(title="")+
  ylab("Z")+
  xlab("Iterations")+
  geom_smooth(color="blue")+
  scale_y_log10()

  ggplot(data = mydata2, aes(x = count, y = Costo_Cupo))+
  labs(title="")+
  ylab("C")+
  xlab("Iterations")+
  geom_smooth(color="blue")+
  scale_y_log10()


ggplot(data = mydata2, aes(x = count, y = Distancia))+
  labs(title="")+
  ylab(bquote(bar(d)))+
  xlab("Iterations")+
  geom_smooth(color="blue")

ggplot(data = mydata2, aes(x = count, y = Segregacion))+
  labs(title="")+
  ylab("D")+
  xlab("Iterations")+
  geom_smooth(color="blue")



#######################
####### Grafica en save
######################
mydata = as.matrix(read.csv("/home/ignacio/Proyectos/S3E2_SAP/sequential_recalentamiento/cmake-build-debug/save/2020-09-29 T:16-11-info-graficos.txt"), header = FALSE, sep=",")
colnames(mydata)[1] <- "count"
colnames(mydata)[2] <- "Distancia"
colnames(mydata)[3] <- "Segregacion"
colnames(mydata)[4] <- "Costo_Cupo"
colnames(mydata)[5] <- "Costo_Solucion"
par(mfrow=c(2,2))
#par(mfrow=c(1,1))
plot(mydata[,"count"],mydata[,"Distancia"], ylab=bquote(bar(d)), xlab="iterations", type="l",col="blue")
plot(mydata[,"count"],mydata[,"Segregacion"], ylab="D", xlab="iterations" , type="l",col="blue")
plot(mydata[,"count"],mydata[,"Costo_Cupo"], ylab="C", xlab="iterations", type="l",col="blue")
plot(mydata[,"count"],mydata[,"Costo_Solucion"], ylab="Z", xlab="iterations", type="l",col="blue")


mydata = as.matrix(read.csv("/home/ignacio/Proyectos/S3E2_SAP/sequential_recalentamiento/cmake-build-debug/save/2020-09-29 T:16-48-info-graficos.txt"), header = FALSE, sep=",")
colnames(mydata)[1] <- "count"
colnames(mydata)[2] <- "Distancia"
colnames(mydata)[3] <- "Segregacion"
colnames(mydata)[4] <- "Costo_Cupo"
colnames(mydata)[5] <- "Costo_Solucion"
par(mfrow=c(2,2))
#par(mfrow=c(1,1))
plot(mydata[,"count"],mydata[,"Distancia"], ylab=bquote(bar(d)), xlab="iterations", type="l",col="blue")
plot(mydata[,"count"],mydata[,"Segregacion"], ylab="D", xlab="iterations" , type="l",col="blue")
plot(mydata[,"count"],mydata[,"Costo_Cupo"], ylab="C", xlab="iterations", type="l",col="blue")
plot(mydata[,"count"],mydata[,"Costo_Solucion"], ylab="Z", xlab="iterations", type="l",col="blue")


#####################################################
############### Version inicial con baja temperatura
#####################################################


mydata = as.matrix(read.csv("/home/ignacio/Proyectos/S3E2_SAP/sequential_recalentamiento/cmake-build-debug/save_example_15_30_25_v3/1-info-graficos.txt"), header = FALSE, sep=",")
colnames(mydata)[1] <- "count"
colnames(mydata)[2] <- "Distancia"
colnames(mydata)[3] <- "Segregacion"
colnames(mydata)[4] <- "Costo_Cupo"
colnames(mydata)[5] <- "Costo_Solucion"
as.numeric(unlist(mydata["count"]))
as.double(unlist(mydata["Distancia"]))
as.double(unlist(mydata["Segregacion"]))
as.double(unlist(mydata["Costo_Cupo"]))
as.double(unlist(mydata["Costo_Solucion"]))
par(mfrow=c(2,2))
#par(mfrow=c(1,1))
plot(mydata[,"count"],mydata[,"Distancia"], ylab=bquote(bar(d)), xlab="iterations", type="l",col="blue")
plot(mydata[,"count"],mydata[,"Segregacion"], ylab="D", xlab="iterations" , type="l",col="blue")
plot(mydata[,"count"],mydata[,"Costo_Cupo"], ylab="C", xlab="iterations", type="l",col="blue")
plot(mydata[,"count"],mydata[,"Costo_Solucion"], ylab="Z", xlab="iterations", type="l",col="blue")



mydata2 = data.frame(mydata) 
ggplot(data = mydata2, aes(x = count, y = Costo_Solucion))+
  labs(title="")+
  ylab("Z")+
  xlab("Iterations")+
  geom_smooth(model="lm") 

ggplot(data = mydata2, aes(x = count, y = Distancia))+
  labs(title="")+
  ylab(bquote(bar(d)))+
  xlab("Iterations")+
  geom_smooth(model="lm") 

ggplot(data = mydata2, aes(x = count, y = Segregacion))+
  labs(title="")+
  ylab("D")+
  xlab("Iterations")+
  geom_smooth(model="lm") 

ggplot(data = mydata2, aes(x = count, y = Costo_Cupo))+
  labs(title="")+
  ylab("C")+
  xlab("Iterations")+
  geom_smooth(model="lm") 


########################################
#### Final graficos 
######################################

mydata = as.matrix(read.csv("/home/ignacio/Proyectos/S3E2_SAP/sequential_recalentamiento/cmake-build-debug/save_example_15_30_25_v3/2-info-graficos.txt"), header = FALSE, sep=",")
colnames(mydata)[1] <- "count"
colnames(mydata)[2] <- "Distancia"
colnames(mydata)[3] <- "Segregacion"
colnames(mydata)[4] <- "Costo_Cupo"
colnames(mydata)[5] <- "Costo_Solucion"
par(mfrow=c(2,2))
#par(mfrow=c(1,1))
plot(mydata[,"count"],mydata[,"Distancia"], ylab=bquote(bar(d)), xlab="iterations", type="l",col="blue")
plot(mydata[,"count"],mydata[,"Segregacion"], ylab="D", xlab="iterations" , type="l",col="blue")
plot(mydata[,"count"],mydata[,"Costo_Cupo"], ylab="C", xlab="iterations", type="l",col="blue")
plot(mydata[,"count"],mydata[,"Costo_Solucion"], ylab="Z", xlab="iterations", type="l",col="blue")




mydata2 = data.frame(mydata) 
ggplot(data = mydata2, aes(x = count, y = Costo_Solucion))+
  labs(title="")+
  ylab("Z")+
  xlab("iterations")+
  geom_smooth(model="lm") 

ggplot(data = mydata2, aes(x = count, y = Distancia))+
  labs(title="")+
  ylab(bquote(bar(d)))+
  xlab("iterations")+
  geom_smooth(model="lm") 

ggplot(data = mydata2, aes(x = count, y = Segregacion))+
  labs(title="")+
  ylab("D")+
  xlab("iterations")+
  geom_smooth(model="lm") 

ggplot(data = mydata2, aes(x = count, y = Costo_Cupo))+
  labs(title="")+
  ylab("C")+
  xlab("iterations")+
  geom_smooth(model="lm") 


ggplot(data=mydata2, aes(x=count))+
  labs(title="")+
  geom_line(aes(y = Costo_Solucion), color="red",linetype="dotdash")+
  ylab("values")+
  geom_line(aes(y = Segregacion), color="blue",linetype="longdash")+
  geom_line(aes(y = Costo_Cupo), color="orange",linetype="twodash")

  
  

########################
########### Comparar la cantidad de alumnos movidos
########################

bestSolution_csv = as.matrix(read.csv("/home/ignacio/Proyectos/S3E2_SAP/sequential_recalentamiento/cmake-build-debug/save_example_15_30_25_v3/2-info-graficos_bestSolution.txt", header = FALSE, sep=","))
bestSolutionV2 = select(as.data.frame(bestSolution_csv),-V29854)


#porrcentajedealumnoscambiados=(length(which(bestSolutionV2[1,]!=bestSolutionV2[2,]))+563)/ncol(bestSolutionV2)
porrcentajedemovimientos=(length(which(bestSolutionV2[1,]!=bestSolutionV2[2,]))+563)/ncol(bestSolutionV2)
porrcentajedemovimientos
segregacion_reducida=0.555358-0.2419633
segregacion_reducida
0.555358-0.294476


school_activates = rep(FALSE,85)

for (x in bestSolutionV2[2,]) {
  school_activates[x+1]=TRUE
}
school_activates


options(scipen = 999)   
sequential.sd_153025_v3
sequential.mean_153025_v3
sequential.median_153025_v3


########################
####### Calcula Duncan
########################
alumnos=as.matrix(read.csv("/home/ignacio/Proyectos/S3E2_SAP/sequential_recalentamiento/cmake-build-debug/alumnos_utm.txt", header = FALSE, sep=","))
alumnos = alumnos[-1,]

totalVuln=length(which(alumnos[,4]==1))
totalNoVuln=length(which(alumnos[,4]!=1))
aluVulCol=0
aluNoVulCol=0
totalsesc=0
for (n in c(0:84)) {
  aluVulCol=0
  aluNoVulCol=0
  #aluNoVulCol=length(which(alumnos[which(bestSolutionV2[1,]==n),4]==0))
  #=length(which(alumnos[which(bestSolutionV2[1,]==n),4]==1))
  aluNoVulCol=length(which(alumnos[which(bestSolutionV2[2,]==n),4]==0))
  aluVulCol=length(which(alumnos[which(bestSolutionV2[2,]==n),4]==1))
  totalsesc = totalsesc+((1/2)*abs((aluVulCol/totalVuln)-(aluNoVulCol/totalNoVuln))  )
  #totalsesc = totalsesc+abs((aluVulCol/totalVuln)-(aluNoVulCol/totalNoVuln))
}
totalsesc
#(1/2)*totalsesc

