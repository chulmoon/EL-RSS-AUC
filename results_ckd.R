library(tidyverse)
library(gridExtra)

# create figures
load("./sim.ckd.res.Rdata")
params=expand.grid(nx=c(20,40,60,80),m=c(2,4))

for (ii in 1:length(sim.res)){
  sim.res[[ii]]$nx = params$nx[ii]
  sim.res[[ii]]$m = params$m[ii]
}

tb.res = NULL
for (ii in 1:length(sim.res)){
  tb.res=rbind(tb.res,c(params$nx[ii],
                        params$m[ii],
                        sim.res[[ii]]$cp[3],
                        sim.res[[ii]]$cp[2],
                        sim.res[[ii]]$cp[1],
                        sim.res[[ii]]$len[3],
                        sim.res[[ii]]$len[2],
                        sim.res[[ii]]$len[1]))
}
tb.res=data.frame(tb.res)
colnames(tb.res)=c("nx","m","srs-cp","ker-cp","el-cp","srs-len","ker-len","el-len")

# plot
## coverage probability
tb.final=NULL
for (nn in c(20, 40, 60, 80)){
  tb.res.temp = tb.res%>%
    filter(nx==nn)%>%
    select(m,"srs-cp","el-cp") 
  
  tb.res.temp.m2 = tb.res.temp %>%
    filter(m==2) %>%
    select(-m) %>%
    rename("SRS-EL"="srs-cp", "RSS-EL (m=2)"="el-cp") %>%
    pivot_longer(contains("-"),names_to = "type",values_to = "cp")
  
  tb.res.temp.m4 = tb.res.temp %>%
    filter(m==4) %>%
    select("el-cp") %>%
    rename("RSS-EL (m=4)"="el-cp") %>%
    pivot_longer(starts_with("RSS"),names_to = "type",values_to = "cp")
  
  tb.final = rbind(tb.final,data.frame(rbind(tb.res.temp.m2,tb.res.temp.m4),n=nn))
}

tb.final=tb.final %>%
  mutate(
    pshape=case_when(
      type == "SRS-EL" ~ 1,
      type == "RSS-EL (m=2)" ~ 2,
      type == "RSS-EL (m=4)" ~ 3,
      TRUE ~ 4
    )
  )
tb.final$pshape=as.factor(tb.final$pshape)

fig1 = tb.final %>%
  ggplot(aes(x=n,y=cp,col=pshape)) +
  geom_line(aes(linetype=pshape))+
  scale_linetype_manual(values=c("solid","longdash",
                                 "dotted","dashed","dotdash"))+
  geom_point(aes(shape=pshape),size=2)+
  scale_shape_manual(values=0:2)+
  labs(x = bquote(n[x]), y="Coverage Probability")+
  scale_y_continuous(limits=c(0.925, 0.975))+
  theme_classic()+
  theme(legend.position="none",
        plot.title = element_text(hjust = 0.5,vjust=-2),
        plot.margin = margin(-1, 0.1, 0.1, 0.1))

## coverage probability
tb.final=NULL
for (nn in c(20, 40, 60, 80)){
  tb.res.temp = tb.res%>%
    filter(nx==nn)%>%
    select(m,"srs-len","el-len") 
  
  tb.res.temp.m2 = tb.res.temp %>%
    filter(m==2) %>%
    select(-m) %>%
    rename("SRS-EL"="srs-len", "RSS-EL (m=2)"="el-len") %>%
    pivot_longer(contains("-"),names_to = "type",values_to = "len")
  
  tb.res.temp.m4 = tb.res.temp %>%
    filter(m==4) %>%
    select("el-len") %>%
    rename("RSS-EL (m=4)"="el-len") %>%
    pivot_longer(starts_with("RSS"),names_to = "type",values_to = "len")
  
  tb.final = rbind(tb.final,data.frame(rbind(tb.res.temp.m2,tb.res.temp.m4),n=nn))
}


tb.final=tb.final %>%
  mutate(
    pshape=case_when(
      type == "SRS-EL" ~ 1,
      type == "RSS-EL (m=2)" ~ 2,
      type == "RSS-EL (m=4)" ~ 3,
      TRUE ~ 4
    )
  )
tb.final$pshape=as.factor(tb.final$pshape)

fig2 = tb.final %>%
  ggplot(aes(x=n,y=len,col=pshape)) +
  geom_line(aes(linetype=pshape))+
  scale_linetype_manual(values=c("solid","longdash",
                                 "dotted","dashed","dotdash"))+
  geom_point(aes(shape=pshape),size=2)+
  scale_shape_manual(values=0:2)+
  labs(x = bquote(n[x]) , y="Length")+
  scale_y_continuous(limits=c(0.12, 0.32))+
  theme_classic()+
  theme(legend.position="none",
        plot.title = element_text(hjust = 0.5,vjust=-2),
        plot.margin = margin(-1, 0.1, 0.1, 0.1))

f.tot=grid.arrange(fig1,fig2,nrow=1)
ggsave("C:/Users/Chul/Box/r_Papers/EL_RSS_for_AUC/result/ckd_result_three.pdf",f.tot,width=6,height=3)

