##################################
  #install.packages("igraph")
library(igraph)

measure<-"volume"
node_info<- read.csv(paste0("/gpfs/gibbs/pi/zhao/js4377/UKB/Cluster/Data/",measure,"_nodes_information.csv"))
gnet<- read.csv(paste0("/gpfs/gibbs/pi/zhao/js4377/UKB/Cluster/Data/",measure,"_network.csv"))
pnet<-  read.csv(paste0("/gpfs/gibbs/pi/zhao/js4377/UKB/Cluster/Data/intercept_",measure,"_network.csv"))
#  randomn<-function(vec num,edge num,clus, rep){
#    modu<- c()
# for(i in 1:rep){
#   random net<-erdos.renyi.game(vec num,edge num, type="gnm",directed = TRUE)
#   modu<-c(modu,modularity(random net,clus))}
#   return(modu)
# }
# cal_p<-function(observed value, data){zscore<-(observed value-mean(data))/sd(data)return(2*(1-pnorm(abs(z score))))

####Given the genetic correlation, return the top per percent connection
data_frame_percentage<- function(gg_df, per){ #col gg_df: Trait1, Trait2, weight
  gg_df<-gg_df[order(-gg_df$weight),]
  return(gg_df[1:round(per*nrow(gg_df)),])
} 

# ##Given the network, return the modularity and random network modularity
# modu<- function(gg_df){ #col gg_df: Trait1, Trait2, weight
#   gg<- data_frame_percentage(gg_df, per)
#   ggraph<- graph_from_data_frame(gg,directed =F)
#   gen_modu<-modularity(ggraph, membership = node_info$cluster_group, weights = gg$weight)
#   #generate random network
# }


gg_df<-gnet[,c("Trait1", "Trait2", "correlation_rg")]
colnames(gg_df)[3]<- "weight"


res_modu<- data.frame(per=numeric(), gm=numeric(), rm=numeric())
per_list<- seq(from = 0.001, to = 0.2, by = 0.001)
for(per in per_list){
  gg<- data_frame_percentage(gg_df, per)
  ggraph<- graph_from_data_frame(gg,directed =F)
  V <- names(as.list(V(ggraph)))
  membership<- node_info[node_info$Node%in%V,]$cluster_group
  gen_modu<-modularity(ggraph, membership = membership, weights = gg$weight)
  ran_modu<-c()
  for(i in 1:100){
    ran_modu<-c(ran_modu,modularity(ggraph, membership = sample(membership), weights = gg$weight))
  }
  
  res_modu<- rbind(res_modu, data.frame(per, gm=gen_modu, rm=mean(ran_modu)))
}

colnames(res_modu)<- c("per", "Genetic Correlations", "Random Networks")
df_long <- tidyr::pivot_longer(res_modu, cols = c("Genetic Correlations", "Random Networks"), 
                               names_to = "type", values_to = "modularity")


p<- ggplot(data = df_long, aes(x = per, y = modularity, color = type)) +
  geom_smooth(se = FALSE, size = 1.5) + # 设置线条粗细
  scale_color_manual(values = c("Genetic Correlations" = "blue", "Random Networks" = "grey"),
                     #labels = c("Genetic Correlations"= "Genetic Correlations", "Random Networks" = "Random Networks")
                     ) +
  labs(
       x = "Fraction of Strongest Connection",
       y = "Modularity",
       color = "Legend") + # 设置图例标题
  theme_minimal() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "top",
        legend.title=element_blank()) # 设置图例位置


#per<- 0.1
# p<-ggplot(data = res_modu, aes(x = per)) +
#   geom_line(aes(y = gm), color = "blue") +
#   geom_line(aes(y = rm), color = "black") +
#   labs(title = "Genetic Correlations Random Networks",
#        x = "Percentage",
#        y = "Modularity") +
#   theme_minimal()
# 
# print(p)
# p<-ggplot(data = res_modu, aes(x = per)) +
#   geom_smooth(aes(y = gm), color = "blue", se = FALSE) +
#   geom_smooth(aes(y = rm), color = "black", se = FALSE) +
#   labs(title = "Genetic Correlations Random Networks",
#        x = "Fraction of Stongest connection",
#        y = "Modularity") +
#   theme_minimal() +
#   theme(panel.grid.major = element_blank(), 
#         panel.grid.minor = element_blank(),
#         panel.background = element_blank(),
#         axis.line = element_line(colour = "black"))

pdf(paste0("/gpfs/gibbs/pi/zhao/js4377/UKB/Cluster/Figure/",measure, "_network.pdf"),width=4, height=4)
print(p)
dev.off()


























set.seed(123) 
random_g <- degree.sequence.game(degree_sequence, method = "vl")
while(!is_connected(random_g)) {
  random_g <- degree.sequence.game(degree_sequence, method = "vl")
}
modularity(random_g, membership = node_info$cluster_group)


pg<- subset(pnet, net!=0)[,c("Trait1","Trait2")]
pgraph<- graph_from_data_frame(pg,directed =F)
pen_modu<- modularity(pgraph, membership = node_info$cluster_group)


##############################
measure <- "area"
node_info<- read.csv(paste0("/gpfs/gibbs/pi/zhao/js4377/UKB/Cluster/Data/",measure,"_nodes_information.csv"))
gnet<- read.csv(paste0("/gpfs/gibbs/pi/zhao/js4377/UKB/Cluster/Data/",measure,"_network.csv"))
pnet<-  read.csv(paste0("/gpfs/gibbs/pi/zhao/js4377/UKB/Cluster/Data/intercept_",measure,"_network.csv"))

data<- data.frame(genetic_cor=gnet$correlation_rg, phenotypic_cor= pnet$Intercept_rg)
correlation <- cor(data$genetic_cor, data$phenotypic_cor)

# # Create the scatter plot
# ggplot(data, aes(x = genetic_cor, y = phenotypic_cor)) +
#   geom_point(alpha = 0.5) + # Add transparency to points
#  # geom_smooth(method = "lm", color = "blue") + # Add linear regression line
#   geom_text(aes(label = paste("r =", round(correlation, 2), ", p < 2.2e-16")), 
#             x = max(data$genetic_cor), y = max(data$phenotypic_cor), hjust = 1, vjust = 1, color = "red") +
#   labs(x = "Interregional genetic correlations", 
#        y = "Interregional phenotypic correlations") +
#   theme_minimal()
# 
# # Print the plot
# print(p)


#Given fraction of the strongest connection, calculate the modularity between the genetic correlation and cluster.

