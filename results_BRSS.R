library(tidyverse)
library(gridExtra)

# Create summary tables and plots for simulation studies using BRSS

# load simulation results
# normal
load("./brss.normal.result.Rdata")
# lognormal
#load("./brss.lognormal.result.Rdata")
# uniform
#load("./brss.uniform.result.Rdata")

params=expand.grid(nx=c(20,40,80),m=c(2,4,5),
                   AUC=c(0.6,0.8,0.9,0.95),corxy=c(1,0.9,0.7))

for (ii in 1:length(sim.res)){
  sim.res[[ii]]$nx = params$nx[ii]
  sim.res[[ii]]$m = params$m[ii]
  sim.res[[ii]]$AUC = params$AUC[ii]
  sim.res[[ii]]$corxy = params$corxy[ii]
}

tb.res = NULL
for (ii in 1:length(sim.res)){
  tb.res=rbind(tb.res,c(params$corxy[ii],params$nx[ii],
                        params$AUC[ii],params$m[ii],
                        sim.res[[ii]]$cp[3],
                        sim.res[[ii]]$cp[2],
                        sim.res[[ii]]$cp[1],
                        sim.res[[ii]]$len[3],
                        sim.res[[ii]]$len[2],
                        sim.res[[ii]]$len[1]))
}

tb.res=data.frame(tb.res)
colnames(tb.res)=c("cor","nx","AUC","m","srs-cp","ker-cp","el-cp","srs-len","ker-len","el-len")

##############################################################
# Create Summary Table
##############################################################
retable = matrix(NA,36,14)
res1 = tb.res %>% 
  filter(cor==1) %>% 
  arrange(nx)
retable[1,1]=res1[1,4]

corval = c(1,0.9,0.7)
nxval = c(20,40,80)
aucval = c(0.6,0.8,0.9,0.95)

for (ii in 1:length(corval)){
  for (jj in 1:length(nxval)){
    for (kk in 1:length(aucval)){
      res.temp = tb.res %>% 
        filter(cor==corval[ii],nx==nxval[jj],AUC==aucval[kk])
      # cp
      retable[(12*(ii-1))+(4*(jj-1))+(kk-1)+1,1]=res.temp[1,5] # SRS-EL
      retable[(12*(ii-1))+(4*(jj-1))+(kk-1)+1,2:4]=res.temp[,6] # RSS-KER
      retable[(12*(ii-1))+(4*(jj-1))+(kk-1)+1,5:7]=res.temp[,7] # RSS-KER
      # cp
      retable[(12*(ii-1))+(4*(jj-1))+(kk-1)+1,7+1]=res.temp[1,8] # SRS-EL
      retable[(12*(ii-1))+(4*(jj-1))+(kk-1)+1,7+2:4]=res.temp[,9] # RSS-KER
      retable[(12*(ii-1))+(4*(jj-1))+(kk-1)+1,7+5:7]=res.temp[,10] # RSS-KER
    }
  }
}

# save as csv files
# normal
#write.csv(retable,"./brss.res.normal.csv")
# lognormal
#write.csv(retable,"./brss.res.lognormal.csv")
# uniform
#write.csv(retable,"./brss.res.uniform.csv")


##############################################################
# Create Plots
##############################################################
## coverage probability
for (nn in c(20, 40, 80)){
  for (rr in c(1, 0.9, 0.7)) {
    tb.res.temp = tb.res%>%
      filter(nx==nn,cor==rr)%>%
      select(AUC,m,"srs-cp","ker-cp","el-cp") 
    
    tb.res.temp.m2 = tb.res.temp %>%
      filter(m==2) %>%
      select(-m) %>%
      rename("SRS-EL"="srs-cp","RSS-KER (m=2)" = "ker-cp", "RSS-EL (m=2)"="el-cp") %>%
      pivot_longer(-AUC,names_to = "type",values_to = "cp")
    
    tb.res.temp.m4 = tb.res.temp %>%
      filter(m==4) %>%
      select(AUC,"ker-cp","el-cp") %>%
      rename( "RSS-KER (m=4)" = "ker-cp", "RSS-EL (m=4)"="el-cp") %>%
      pivot_longer(-AUC,names_to = "type",values_to = "cp")
    
    
    tb.final = data.frame(rbind(tb.res.temp.m2,tb.res.temp.m4))
    tb.final$AUC=ifelse(tb.final$AUC==0.6,1,ifelse(tb.final$AUC==0.8,2,ifelse(tb.final$AUC==0.9,3,4)))
    
    tb.final=tb.final %>%
      mutate(
        pshape=case_when(
          type == "SRS-EL" ~ 1,
          type == "RSS-EL (m=2)" ~ 2,
          type == "RSS-EL (m=4)" ~ 3,
          type == "RSS-KER (m=2)" ~ 4,
          type == "RSS-KER (m=4)" ~ 5,
          TRUE ~ 6
        )
      )
    tb.final$pshape=as.factor(tb.final$pshape)
    
    if (rr==1) {
      fig = tb.final %>%
        ggplot(aes(x=AUC,y=cp,col=pshape)) +
        geom_line(aes(linetype=pshape))+
        scale_linetype_manual(values=c("solid","longdash",
                                       "dotted","dashed","dotdash"))+
        geom_point(aes(shape=pshape),size=2)+
        scale_shape_manual(values=0:4)+
        labs(tx ="AUC", y="Coverage Probability", title=bquote(n[x] == .(nn) ~ rho == .(rr)) )+
        scale_y_continuous(limits=c(0, 1))+
        scale_x_continuous(breaks = 1:4, labels = c("0.6","0.8","0.9","0.95")) +
        theme_classic()+
        theme(legend.position="none",
              plot.title = element_text(hjust = 0.5,vjust=-2),
              plot.margin = margin(1.5,1.5,1.5,1.5)
        )
    } else{
      fig = tb.final %>%
        ggplot(aes(x=AUC,y=cp,col=pshape)) +
        geom_line(aes(linetype=pshape))+
        scale_linetype_manual(values=c("solid","longdash",
                                       "dotted","dashed","dotdash"))+
        geom_point(aes(shape=pshape),size=2)+
        scale_shape_manual(values=0:4)+
        labs(tx ="AUC", y="", title=bquote(n[x] == .(nn) ~ rho == .(rr)))+
        scale_y_continuous(limits=c(0, 1))+
        scale_x_continuous(breaks = 1:4, labels = c("0.6","0.8","0.9","0.95")) +
        theme_classic()+
        theme(legend.position="none",
              plot.title = element_text(hjust = 0.5,vjust=-2),
              plot.margin = margin(1.5,1.5,1.5,1.5)
        )
    }
    assign(paste0("c",nn,rr),fig)
  }
}

## length of confidence interval
for (nn in c(20, 40, 80)){
  for (rr in c(1, 0.9, 0.7)) {
    tb.res.temp = tb.res%>%
      filter(nx==nn,cor==rr)%>%
      select(AUC,m,"srs-len","ker-len","el-len") 
    
    tb.res.temp.m2 = tb.res.temp %>%
      filter(m==2) %>%
      select(-m) %>%
      rename("SRS-EL"="srs-len","RSS-KER (m=2)" = "ker-len", "RSS-EL (m=2)"="el-len") %>%
      pivot_longer(-AUC,names_to = "type",values_to = "len")
    
    tb.res.temp.m4 = tb.res.temp %>%
      filter(m==4) %>%
      select(AUC,"ker-len","el-len") %>%
      rename( "RSS-KER (m=4)" = "ker-len", "RSS-EL (m=4)"="el-len") %>%
      pivot_longer(-AUC,names_to = "type",values_to = "len")
    
    
    tb.final = data.frame(rbind(tb.res.temp.m2,tb.res.temp.m4))
    tb.final$AUC=ifelse(tb.final$AUC==0.6,1,ifelse(tb.final$AUC==0.8,2,ifelse(tb.final$AUC==0.9,3,4)))
    
    tb.final=tb.final %>%
      mutate(
        pshape=case_when(
          type == "SRS-EL" ~ 1,
          type == "RSS-EL (m=2)" ~ 2,
          type == "RSS-EL (m=4)" ~ 3,
          type == "RSS-KER (m=2)" ~ 4,
          type == "RSS-KER (m=4)" ~ 5,
          TRUE ~ 6
        )
      )
    tb.final$pshape=as.factor(tb.final$pshape)
    
    if (rr==1) {
      fig = tb.final %>%
        ggplot(aes(x=AUC,y=len,col=pshape)) +
        geom_line(aes(linetype=pshape))+
        scale_linetype_manual(values=c("solid","longdash",
                                       "dotted","dashed","dotdash"))+
        geom_point(aes(shape=pshape),size=2)+
        scale_shape_manual(values=0:4)+
        labs(tx ="AUC", y="Length", title=bquote(n[x] == .(nn) ~ rho == .(rr)))+
        scale_y_continuous(limits=c(0, 0.4))+
        scale_x_continuous(breaks = 1:4, labels = c("0.6","0.8","0.9","0.95")) +
        theme_classic()+
        theme(legend.position="none",
              plot.title = element_text(hjust = 0.5,vjust=-2),
              plot.margin = margin(1.5,1.5,1.5,1.5))
    } else{
      fig = tb.final %>%
        ggplot(aes(x=AUC,y=len,col=pshape)) +
        geom_line(aes(linetype=pshape))+
        scale_linetype_manual(values=c("solid","longdash",
                                       "dotted","dashed","dotdash"))+
        geom_point(aes(shape=pshape),size=2)+
        scale_shape_manual(values=0:4)+
        labs(tx ="AUC", y="", title=bquote(n[x] == .(nn) ~ rho == .(rr)))+
        scale_y_continuous(limits=c(0, 0.4))+
        scale_x_continuous(breaks = 1:4, labels = c("0.6","0.8","0.9","0.95")) +
        theme_classic()+
        theme(legend.position="none",
              plot.title = element_text(hjust = 0.5,vjust=-2),
              plot.margin = margin(1.5,1.5,1.5,1.5))
    }
    assign(paste0("l",nn,rr),fig)
  }
}
fig.total=grid.arrange(c201,c200.9,c200.7,c401,c400.9,c400.7,c801,c800.9,c800.7,l201,l200.9,l200.7,l401,l400.9,l400.7,l801,l800.9,l800.7,nrow=6)
