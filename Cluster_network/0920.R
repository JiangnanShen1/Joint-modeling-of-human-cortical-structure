####ggnetwork
library(ggplot2)
library(geomnet)
library(ggrepel)
library(sna)
library(dplyr)


volume<-read.csv(".csv")
#area
volume$group <- c(1, 2, 2, 3, 1, 1, 1, 1, 3, 3, 2, 3, 2, 1, 1, 2 ,2, 2, 2, 3 ,1 ,2, 2, 1, 2, 2, 2, 1, 2, 1, 2, 2 ,2 ,1 ,2, 2 ,3, 1, 1 ,1 ,1, 3, 3, 2, 3, 2, 1 ,1, 2 ,2, 2, 2, 3, 1, 2, 2, 1, 2, 2, 2, 1, 2, 1, 2, 2, 2)
#volume
#volume$group<-c(2,3,3,4,5,2,2,2,4,4,3,4,3,2,5,3,3,3,3,4,1,3,3,1,3,3,3,1,3,1,3,3,3,2,3,3,4,5,2,2,2,4,4,3,4,3,2,5,3,3,3,3,4,1,3,3,1,3,3,3,1,3,1,3,3,3)
net<-read.csv("network.csv")
net<-net%>%select(-X)%>%filter(net==1)%>%filter(Trait1!=Trait2)%>%left_join(volume,by=c("Trait1"="Node"))

#volume$group2 <- c(1,2,3,4,1,1,5,1,2,4,3,4,3,1,1,3,3,3,3,4,5,2,3,5,2,3,3,5,4,5,3,1,2)


net <- net[1:(nrow(net)/2),]
aaddd<-matrix(0,66,66)
for (i in 1:length(seq(344,409))) {
  n1<-net%>%filter(Trait1==seq(344,409)[i])%>%select(Trait2)
  aaddd[n1$Trait2-343,i]=1
}

set.seed(10312016)
guse<-gplot(aaddd)
x<-guse[,1]
y<-guse[,2]


x1<-(x-min(x))/(max(x)-min(x))
y1<-(y-min(y))/(max(y)-min(y))



net$group<-as.character(net$group)
net$group2<-as.character(net$group2)
net$Trait1<-as.character(net$Trait1)
net$Trait2<-as.character(net$Trait2)
label=unique(net$Trait1)


net$color_net <- as.numeric(net$group)
net$color_net[which(net$color_net==1)] <- "PaleGreen3"
net$color_net[which(net$color_net==2)] <- "LightBlue2"
net$color_net[which(net$color_net==3)] <- "PaleVioletRed1"

volume$color_net <- volume$group
volume$color_net[which(volume$color_net==1)] <- "PaleGreen3"
volume$color_net[which(volume$color_net==2)] <- "LightBlue2"
volume$color_net[which(volume$color_net==3)] <- "PaleVioletRed1"

set.seed(10312016)
ggplot() + 
  geom_net(data=net,aes(from_id =Trait1, to_id = Trait2,group=group,colour =group,shape=as.character(group2)),layout.alg = 'fruchtermanreingold',
           linewidth = 0, ealpha = 0.25,  curvature = 0.05,
           size = 3, vjust = -0.8) +
  geom_text_repel(aes(x1,y1),color=as.character(volume$color_net),label=label,size=3)+ 
  theme_net() +
  theme(legend.position = "bottom") +
  scale_colour_manual('Trait group',  values=c("PaleGreen3","LightBlue2","PaleVioletRed1")) +
  scale_shape_manual(values=c(15,16,17,18,19))
  xlim(c(-0.1,1.1))
  
ss <-    ggplot() + 
    geom_net(data=net,aes(from_id =Trait1, to_id = Trait2,group=group,colour =group),layout.alg = 'fruchtermanreingold',
             linewidth = 0.3, ealpha = 0.25,  curvature = 0.05,
             size = 4, vjust = -0.8) +
    geom_text_repel(aes(x1,y1),color=as.character(volume$color_net),label=label,size=3)+ 
    theme_net() +
    theme(legend.position = "bottom") +
    scale_colour_manual('Trait group',  values=c("PaleGreen3","LightBlue2","PaleVioletRed1")) +
  xlim(c(-0.1,1.1))

ss2 <-  ggplot() + 
  geom_net(data=net,aes(from_id =Trait1, to_id = Trait2,group=group,colour =as.character(group2)),layout.alg = 'fruchtermanreingold',
           linewidth = 0.3, ealpha = 0,  curvature = 0.05,
           size = 2, vjust = -0.8,ecolour="white") +
  geom_text_repel(aes(x1,y1),color="red",alpha=0,label=label,size=3)+ 
  theme_net() +
  theme(legend.position = "bottom") +
  scale_colour_manual(values=c("PaleGreen3","LightBlue2","PaleVioletRed1","red","blue")) +
  xlim(c(-0.1,1.1))


 
 set.seed(10312016)
 vp <- viewport()
 print(ss)
 set.seed(10312016)
 print(ss2,vp=vp) 
 
 
 #################################
 ggsave("ss.jpg",ss)
 sss <- load.image("ss.jpg")
 
 theme_transp_overlay <- theme(
   panel.background = element_rect(fill="transparent",color=NA),
   plot.background = element_rect(fill="transparent",color=NA)
 )
 
 set.seed(10312016)
 plot(ss)+annotation_custom(ggplotGrob(ss2))