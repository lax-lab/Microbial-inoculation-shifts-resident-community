setwd("E:/onedrive/小论文/光合菌Rre")
#phylum求和
data=read.table("clipboard",row.names = 1,header=T,sep='\t')
sum_phylum=aggregate( .~ phylum, data = data, sum)
sum_phylum=sum_phylum[,1:13]
write.csv(sum_phylum,"sum_phylum.csv")

##relative abundance
data=read.table("clipboard",header=T,sep='\t')
data= melt(data, id=1)
data$variable_o= factor(data$variable, levels = c('Control','R.palustris','B.subtilis', 'Mix'))
data$phylum_o= factor(data$phylum, levels = c('Others','Gemmatimonadetes','Ignavibacteriae','Cyanobacteria','Planctomycetes','Bacteroidetes','Actinobacteria','Firmicutes','Acidobacteria','Chloroflexi', 'Proteobacteria'))

phylum_colors <- c("#E280B6","#808080","#2F528F","#9A0000","#BA74B7","#8E7900","#E1B1CC","#E4932D","#7F83D0", "#4DB14C","#A09FDB","#BD38A3","#A5ADB8","#BD38A3")
p<-ggplot(data,aes(variable_o,value,fill=phylum_o))+
  geom_bar(stat="identity")+
  scale_fill_manual(values = phylum_colors)+
  theme_bw()+ylab("Relative ahundance")+
  theme(axis.text.x = element_text(size=13,angle = 30, hjust= 1, vjust = 1),axis.text.y = element_text(size=13))+
  theme(axis.title.y = element_text(size = 15))
graph2ppt(file="figs by xiao.pptx", width=4, height=3,append =TRUE)
##微生物酶活
data=read.table("clipboard",header=T,sep='\t')

data$Treatment_o = factor(data$Treatment, levels=c('Control','R.palustris','B.subtilis', 'Mix'))
data1=subset(data,Enzyme== "Shannon diversity")
p1=ggplot(data1,aes(x=Treatment_o,y=Mean))+
  geom_bar(stat="identity",width =0.7,fill= "DimGray" )+
  geom_errorbar(aes(ymin=Mean-Errorbar, ymax=Mean+Errorbar), width=.2,
                position=position_dodge(.9))+
  theme_bw()+
  theme(axis.text.x = element_text(size=13,angle = 30, hjust= 1, vjust = 1),axis.text.y = element_text(size=13))+
  theme(axis.title.y = element_text(size = 15))+
  theme(plot.title = element_text(size = 16))+
  xlab("")+ylab("Shannon diversity")+
  coord_cartesian(ylim=c(6,7.5))
p1
graph2ppt(file="figs by xiao.pptx", width=2.5, height=3,append =TRUE)

##NMDS
data=read.table("clipboard",header=T,sep='\t')
data$Group1_o= factor(data$Group1, levels = c('Control','R.palustris','B.subtilis','Mix'))
ggplot(data, aes(x=NMDS1, y=NMDS2, color=Group1_o)) +
  geom_point(size=4.5)+theme_bw()+
  theme(panel.grid.major = element_blank (), panel.grid.minor = element_blank (), axis.line = element_line (colour = "black"))+
  theme(axis.text.x = element_text(size = 15),axis.text.y = element_text(size = 15))+
  theme(axis.title.x = element_text(size =15),axis.title.y = element_text(size = 15))+
  theme(legend.text=element_text(size=11))+
  scale_shape_manual(values=c(1, 16))+
  scale_color_manual(values=c('#440154','#3B528B','#5DC863','#FDE725'))
graph2ppt(file="figs by xiao.pptx", width=4.5, height=3, append=TRUE)


data=data.matrix(data)
dune.env=read.table("clipboard",header = T,sep='\t',row.names=1)
dune.dist <- vegdist(otu)
attach(dune.env)
dune.ano <- anosim(dune.dist, dune.env$Group2)
summary(dune.ano)
plot(dune.ano)
(sim <- with(dune.env, simper(otu, dune.env$Group2)))
a=summary(sim,ordered = FALSE )
write.csv(a$Rhizophere_Rhizoplane,"3.csv")

##SEM
data=read.table("clipboard",header=T,sep='\t')
data$ID_o= factor(data$ID, levels = c('SOM','Soil nitrogen','Alpha diversity','Beta diversity','Microbial activity'))
p1=ggplot(data,aes(x=ID_o,y=s4))+
  geom_bar(stat="identity",width =0.7,fill="#B60000" )+
  theme_bw()+
  theme(axis.text.x = element_text(size=34,angle = 34, hjust= 1, vjust = 1),axis.text.y = element_text(size=13))+
  theme(axis.title.y = element_text(size = 36))+
  theme(plot.title = element_text(size = 34))+
  xlab("")+ylab("Increase in MSE (%)")
p1
graph2ppt(file="SEM.pptx", width=9, height=8, append=TRUE)

##OD##
data=read.table("clipboard",header=T,sep='\t')
data$ID_o = factor(data$ID, levels=c('R.palustris','B.subtilis', 'Mix'))
p1=ggplot(data,aes(x=ID_o,y=mean))+
  geom_bar(stat="identity",width =0.7,fill= "DimGray" )+
  geom_errorbar(aes(ymin=mean-error, ymax=mean+error), width=.2,
                position=position_dodge(.9))+
  theme_bw()+
  theme(axis.text.x = element_text(size=13,angle = 30, hjust= 1, vjust = 1),axis.text.y = element_text(size=13))+
  theme(axis.title.y = element_text(size = 15))+
  theme(plot.title = element_text(size = 16))+
  xlab("")+ylab("OD (660 nm)")+
  coord_cartesian(ylim=c(0,3))
p1
graph2ppt(file="figs by xiao.pptx", width=2.5, height=3,append =TRUE)

data=read.table("clipboard",header=T,sep='\t')
data$Treatment_o = factor(data$Treatment, levels=c('R.palustris','B.subtilis'))
p1=ggplot(data,aes(x=Treatment_o,y=Mean))+
  geom_bar(stat="identity",width =0.7,fill= "DimGray" )+
  geom_errorbar(aes(ymin=Mean-Error, ymax=Mean+Error), width=.2,
                position=position_dodge(.9))+
  facet_grid(variable~.,scales = "free")+
  theme_bw()+
  theme(axis.text.x = element_text(size=13),axis.text.y = element_text(size=13))+
  theme(axis.title.y = element_text(size = 15))+
  theme(plot.title = element_text(size = 16))+
  xlab("")+ylab("")
graph2ppt(file="figs by xiao.pptx", width=2, height=3,append =TRUE)

##linear regression###
data=read.table("clipboard",header=T,sep='\t')
ggplot(data, aes(x = Mix,y=Root.dry.weight)) + 
  geom_point(size=3,alpha=0.6) + 
  theme_bw()+
  theme(axis.text.x = element_text(size=18),axis.text.y = element_text(size=18))+
  theme(axis.title.x = element_text(size = 20),axis.title.y = element_text(size =20))+
  xlab("Abudance of indicator species")+ylab("Yield (t/ha)")+
  geom_smooth(method = lm, formula = y ~ x)

  coord_cartesian(ylim=c(0.12,0.26))
 
graph2ppt(file="linear.pptx", width=3, height=2.5,append =TRUE)
