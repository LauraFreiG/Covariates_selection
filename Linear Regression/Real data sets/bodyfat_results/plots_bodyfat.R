
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



pdf("barplots_bodyfat_ggplot.pdf", width=10*2/3, height=2)
library(ggplot2) 
library(gridExtra)
dat=data.frame( names=c("LASSO.min","LASSO.1se","LASSO.BIC","AdapL.min",
                        "AdapL.1se","SCAD","Dant","RelaxL", "SqrtL", 
                        "ScalL", "DC.VS"),
                sel=sel_w )
dat$names=factor( dat$names, levels=dat$names )
colores=c()
for(kk in 1:11){
  if( sel_w[kk] <=3 ){
    colores[kk]="lightsalmon"
  } else if( 4<=sel_w[kk] & sel_w[kk]<=8 ){
    colores[kk]="mediumpurple1"
  } else {
    colores[kk]="steelblue1"
  }
}
p<-ggplot(data=dat, aes(x=sel, y=names)) +  geom_bar(stat="identity", fill=colores) +
  geom_text(aes(label=sel), hjust=1.8, color="white", size=3.2) + theme_minimal() +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
        axis.text.x=element_blank(), axis.text.y=element_text(size=9), 
        plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
        coord_cartesian(xlim = c(0, 13.35))
#p

dat_uni=data.frame( names=c("LASSO.min","LASSO.1se","LASSO.BIC","AdapL.min",
                        "AdapL.1se","SCAD","Dant","RelaxL", "SqrtL", 
                        "ScalL", "DC.VS"),
                sel=sel_uni )
dat_uni$names=factor( dat_uni$names, levels=dat_uni$names )
colores=c()
for(kk in 1:11){
  if( sel_uni[kk] <=3 ){
    colores[kk]="lightsalmon"
  } else if( 4<=sel_uni[kk] & sel_uni[kk]<=8 ){
    colores[kk]="mediumpurple1"
  } else {
    colores[kk]="steelblue1"
  }
}
p_uni<-ggplot(data=dat_uni, aes(x=sel, y=names)) +  geom_bar(stat="identity", fill=colores) +
  geom_text(aes(label=sel), hjust=1.8, color="white", size=3.2) + theme_minimal() +
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(),
        axis.text.x=element_blank(), axis.text.y=element_text(size=9), 
        plot.title = element_text(hjust = 0.5),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
        coord_cartesian(xlim = c(0, 13.35))

grid.arrange(p, p_uni, ncol=2)
dev.off()
