
#load data
load("sel.RData")
sel_w=sel
load("sel_stand.RData")
sel_uni=sel

load("cov.RData")
cov_w=cov_sel
load("cov_stand.RData")
cov_uni=cov_sel
names(cov_w)=names(cov_uni)=c("LASSO.min","LASSO.1se","LASSO.BIC","AdapL.min",
                              "AdapL.1se","SCAD","Dant","RelaxL", "SqrtL", 
                              "ScalL", "DC.VS")






pdf("barplots_riboflavin_ggplot.pdf", width=10, height=2.5)

library(ggplot2) 
library(gridExtra)
sel=data.frame( names=c("LASSO.min","LASSO.1se","LASSO.BIC","AdapL.min",
                        "AdapL.1se","SCAD","Dant","RelaxL", "SqrtL", 
                        "ScalL", "DC.VS"),
                sel=sel_w )
sel$names=factor( sel$names, levels=sel$names )
colores=c()
for(kk in 1:11){
  if( sel_w[kk] <=5 ){
    colores[kk]="lightsalmon"
  } else if( 6<=sel_w[kk] & sel_w[kk]<14 ){
    colores[kk]="mediumpurple1"
  } else {
    colores[kk]="steelblue1"
  }
}
col_met=c("black","black","black","black","black","black","black","black","black","black","black")
p<-ggplot(data=sel, aes(x=sel, y=names)) +  geom_bar(stat="identity", fill=colores) +
  geom_text(aes(label=sel), hjust=1.3, color="white", size=4) + theme_minimal() +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
        axis.text.x=element_blank(), axis.text.y=element_text(size=12,colour=col_met), 
        plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
#p

sel_stand=data.frame( names=c("LASSO.min","LASSO.1se","LASSO.BIC","AdapL.min",
                              "AdapL.1se","SCAD","Dant","RelaxL", "SqrtL", 
                              "ScalL", "DC.VS"),
                      sel=sel_uni )
sel_stand$names=factor( sel_stand$names, levels=sel_stand$names )
colores=c()
for(kk in 1:11){
  if( sel_uni[kk] <=5 ){
    colores[kk]="lightsalmon"
  } else if( 6<=sel_uni[kk] & sel_uni[kk]<=14 ){
    colores[kk]="mediumpurple1"
  } else {
    colores[kk]="steelblue1"
  }
}
# colores=c( rep("steelblue1",3), "mediumpurple1", "lightsalmon", "steelblue1", "mediumpurple1", rep("steelblue1",4))
p_stand<-ggplot(data=sel_stand, aes(x=sel, y=names)) +  geom_bar(stat="identity", fill=colores) +
  geom_text(aes(label=sel), hjust=1.3, color="white", size=4) + theme_minimal() +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
        axis.text.x=element_blank(), axis.text.y=element_text(size=12, hjust=1.1, colour=col_met), 
        plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank()) 
#p_stand

grid.arrange(p, p_stand, ncol=2)
dev.off()

#-----------------------------------------------------------------------------------------

