##Function for plotting QQplots of gene expression profiles
qqplot_geneProfile<-function(val1,val2,method1,method2){
    a<-cor.test(val1,val2)
    rho<-round(a$estimate,digits=2)
    dat=as.data.frame(cbind(val1,val2))
    
    ggplot(dat, aes(x=val1, y=val2)) +    
        geom_point(aes(color="grey60"), show.legend=F, size=2)+
        geom_smooth(method = lm, color="black", se=F)+
        theme_classic () +
        scale_color_manual(values="grey60")+
        ggtitle("Pearson correlation of Gene-Expression profiles")+
        labs(title="", x=paste0(method1," "), y = paste0(method2, " "))+
        annotate("text", x=15, y=10, label=paste0("rho = ",rho), size=4)+
        theme(plot.title = element_text(face="bold", colour = "black", size=9))+
        theme(axis.title.y = element_text(margin = margin(r=0.2, unit = "cm"),size=12, vjust = -0.5))+
        theme(axis.title.x = element_text(margin = margin(t=0.2, unit = "cm"),size=12))+
        theme(axis.text.x = element_text(colour = "black", size = 10))+
        theme(axis.text.y = element_text(colour = "black", size = 10))+
        theme(plot.margin = margin(t = 0.2, r = 0.1, b = 0.2, l = 0.1, unit = "cm"))
}  