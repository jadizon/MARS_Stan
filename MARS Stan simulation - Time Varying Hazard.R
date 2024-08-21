rm(list=ls())
library(tidyverse)
library(cmdstanr)
library(data.table)
library(survival)
library(broom)
library(here)
Sys.setenv(Home="C:/Users/...")
stan_file_path = file.path(here("MARS_Stan.stan"))##NOTE: Put Stan file somewhere in the C:/ drive to avoid errors
options(error=NULL)




simulation<-function(repeats){
 time_start<<-timestamp()
 sim_summary<<-list()
 for(repz in 1:repeats){
 clusters<-20
 r.tvh <- matrix(c(.01000, 0.0000, 
 0.0000, .01000), 
 nrow = 2, ncol = 2, byrow = TRUE)*
 matrix(c(1.0000, 0.3000, 
 0.3000, 1.0000), 
 nrow = 2, ncol = 2, byrow = TRUE) *
 matrix(c(.01000, 0.0000, 
 0.0000, .01000), 
 nrow = 2, ncol = 2, byrow = TRUE)
 
 rand_eff <- MASS::mvrnorm(clusters, mu = c(0, 0), Sigma = r.tvh)
 sim.vars.tvh <- list()
 sim.vars.tvh$study_id=c(1:clusters)
 sim.vars.tvh$sample_size=sample.int(100, clusters, replace=TRUE)+50
 sim.vars.tvh$followup_time=5*sample.int(8, clusters, replace=TRUE)+80
 sim.vars.tvh$alpha=rand_eff[, 1]
 sim.vars.tvh$beta=rand_eff[, 2]
 sim.vars.tvh$b=-0.60
 sim.vars.tvh$total_sample=sum(sim.vars.tvh$sample_size)
 
 
 
 
 sim.data.tvh. <- data.table(rownum=1:sum(sim.vars.tvh$sample_size), 
 study_id=rep(sim.vars.tvh$study_id, sim.vars.tvh$sample_size), 
 end=rep(sim.vars.tvh$followup_time, sim.vars.tvh$sample_size), 
 alpha=rep(sim.vars.tvh$alpha, sim.vars.tvh$sample_size), 
 beta=rep(sim.vars.tvh$beta, sim.vars.tvh$sample_size), 
 x=rbinom(sim.vars.tvh$total_sample, 1, 0.55), 
 censor_time = -log(runif(sim.vars.tvh$total_sample))/0.01
 )%>% 
 data.frame()
 
 
 for (i in 1:dim(sim.data.tvh.)[1]){
 if (i==1) sim.data.tvh <- sim.data.tvh.
 time<-0
 proceed<-1
 x<-sim.data.tvh[["x"]][i]
 alpha<-sim.data.tvh[["alpha"]][i]
 beta<-sim.data.tvh[["beta"]][i]
 censor_time<-sim.data.tvh[["censor_time"]][i]
 end_time<-sim.data.tvh[["end"]][i]
 for (t in 1:10){
 if (proceed == 0) break
 u<-runif(1)
 b<-1-0.175^(t*12/45-1)
 survtime<- -log(u)/(0.01*exp(alpha + (b+beta)*x) );
 if (survtime > 12){
 if(t==10) time = time+survtime
 else time = time+12
 proceed <-1
 }
 else {
 time = survtime + time
 proceed <-0
 }
 sim.data.tvh[["interval"]][i] <- t 
 }
 sim.data.tvh[["death"]][i] <- ifelse(time<min(censor_time, end_time), 1, 0)
 sim.data.tvh[["time"]][i] <- ifelse(sim.data.tvh[["death"]][i]==1, time, min(censor_time, end_time))
 sim.data.tvh[["censor"]][i] <- ifelse(sim.data.tvh[["death"]][i]==1, 0, 1)
 rm(u, b, survtime, proceed, x, t, time, alpha, beta, i, censor_time, end_time)
 } 
 k<<-sim.data.tvh
 summary(survfit(Surv(time, death) ~ factor(x), data=sim.data.tvh), times = c(12, 36, 60, 120))
 # time n.risk n.event survival std.err lower 95% CI upper 95% CI
 # 12 35540 4841 0.887 0.00154 0.884 0.890
 # 36 22011 6757 0.698 0.00238 0.693 0.702
 # 60 13659 4144 0.550 0.00278 0.545 0.556
 # 120 513 3857 0.301 0.00482 0.292 0.310
 # 
 # factor(x)=1 
 # time n.risk n.event survival std.err lower 95% CI upper 95% CI
 # 12 48125 501 0.990 0.000431 0.989 0.991
 # 36 33744 4442 0.885 0.001548 0.882 0.888
 # 60 19256 8206 0.642 0.002565 0.637 0.647
 # 120 344 9520 0.160 0.004370 0.152 0.169
 ########################################################################################
 ########################################################################################
 ########################################################################################
 
 
 sim.tvh.t1 <- sim.data.tvh[sim.data.tvh$study_id<=9, ]
 
 
 ########################################################################################
 ########################################################################################
 ########################################################################################
 
 
 sim.tvh.hr <- sim.data.tvh[sim.data.tvh$study_id>9 & sim.data.tvh$study_id<=13, ]
 
 t2_studies.tvh <- min(sim.tvh.hr$study_id):max(sim.tvh.hr$study_id)
 
 for (i in t2_studies.tvh){
 if (i==min(t2_studies.tvh)) {
 sim.tvh.t2 <- data.table()
 }
 cox <- data.table(coxph(Surv(time, death) ~ factor(x), data=sim.tvh.hr[sim.tvh.hr$study_id==i, ])|> 
 tidy(conf.int = TRUE, exponentiate = F) |> 
 select(term, estimate, starts_with("conf")))
 
 sim.tvh.t2 <- rbind(sim.tvh.t2, cox)
 
 
 
 if(i==max(sim.tvh.hr$study_id)){
 sim.tvh.t2 <- data.frame(study_id=sim.vars.tvh$study_id[t2_studies.tvh], followup_time=sim.vars.tvh$followup_time[t2_studies.tvh], sim.tvh.t2)
 colnames(sim.tvh.t2) <- c("study_id", "time", "", "log_hr", "log_hr_l", "log_hr_u")
 }
 
 };rm(i)
 
 
 
 ########################################################################################
 ########################################################################################
 ########################################################################################
 
 
 
 
 t3.df <- sim.data.tvh[sim.data.tvh$study_id>13, ]
 t3_studies.tvh <- min(t3.df$study_id):max(t3.df$study_id)
 t3.dfa <- sim.data.tvh[sim.data.tvh$study_id>13, ]
 
 
 for (k in t3_studies.tvh){
 if(k==min(t3_studies.tvh)){
 t3.rates <- list()
 t3.studies <-list()
 t3.followup_time <-list()
 t3.group <-list()
 t3.se<-list()
 
 
 }
 
 survfit_x0 <- survfit(Surv(time, death) ~ factor(x), data=t3.dfa[t3.dfa$study_id==k & t3.dfa$x==0, ])
 survfit_x1 <- survfit(Surv(time, death) ~ factor(x), data=t3.dfa[t3.dfa$study_id==k & t3.dfa$x==1, ])
 surv_x0<- summary(survfit_x0, times = c(12, 36, 60, min(120, max(survfit_x0$time))))
 surv_x1<- summary(survfit_x1, times = c(12, 36, 60, min(120, max(survfit_x1$time))))
 t3.rates<-append(t3.rates, c(unlist(surv_x0$surv), unlist(surv_x1$surv)))
 t3.followup_time<-append(t3.followup_time, c(unlist(surv_x0$time), unlist(surv_x1$time)))
 t3.se<-append(t3.se, c(unlist(surv_x0$std.err), unlist(surv_x1$std.err)))
 t3.studies <-append(t3.studies, rep(k, 8))
 t3.group <-append(t3.group, c(rep(0, length(surv_x0$time)), rep(1, length(surv_x1$time))))
 
 }
 sim.tvh.t3<-data.frame(study_id=unlist(t3.studies), 
 x=unlist(t3.group), 
 event_prob=unlist(t3.rates), 
 se=unlist(t3.se), 
 followup_time=unlist(t3.followup_time)
 )%>%
 mutate(logit=log(event_prob/(1-event_prob)))%>%
 mutate(logit_se=(se^2)/(((event_prob)^2)*((1-event_prob)^2)))%>%
 group_by(study_id, x)%>%
 filter(followup_time== max(followup_time))
 sim.tvh.t3<-sim.tvh.t3[sim.tvh.t3$event_prob !=0, ]
 rm(t3.df, t3.dfa, t3.followup_time, t3.rates, t3.group, t3.se, t3.studies, survfit_x0, survfit_x1, surv_x0, surv_x1)
 
 
 ########################################################################################
 ########################################################################################
 ########################################################################################
 
 J=10
 J_step=12
 
 T1.os <- sim.tvh.t1
 T2.os <- sim.tvh.t2
 T3.os <- sim.tvh.t3
 
 idlist <- sort(unique(c(T1.os$study_id, T2.os$study_id, T3.os$study_id)))
 T1.os$studyid <- sapply(T1.os$study_id, function(x) which(idlist==x))
 T2.os$studyid <- sapply(T2.os$study_id, function(x) which(idlist==x))
 T3.os$studyid <- sapply(T3.os$study_id, function(x) which(idlist==x))
 n_Study <- length(idlist)
 T1.os<-data.frame(T1.os)
 T2.os<-data.frame(T2.os)
 T3.os<-data.frame(T3.os)
 
 idlist2 <- data.frame(studyid=c(sort(unique(c(T1.os$studyid, T3.os$studyid))), sort(unique(c(T2.os$studyid)))), new_ids=1:n_Study)
 
 T1.os$s_id1 <- merge(T1.os, idlist2, by="studyid")[, "new_ids"]
 T1.os$t1 <- ceiling((T1.os$time/12))
 T1.os.long <- data.frame(T1.os)%>%
 rowwise()%>%
 slice(rep(1:n(), each = t1))%>%
 data.frame()%>%
 mutate(t1_true=t1)%>%
 arrange(rownum)%>%
 group_by(rownum)%>%
 mutate(t1=row_number())%>%
 mutate(censor=ifelse(t1==t1_true, censor, 0))%>%
 mutate(death=ifelse(t1==t1_true, death, 0))%>%
 rowwise()%>%
 mutate(jt1=t1*12)%>%
 mutate(jt2=(t1-1)*12)%>%
 mutate(time2=min(t1*12, time))%>%
 mutate(time_in_interval=min(time2-jt2, jt1-jt2))%>%
 mutate(time_in_interval=ifelse(time_in_interval<0, 0, time_in_interval))%>%
 mutate(time_in_interval=time_in_interval/12)
 
 
 
 T2.os$s_id2 <- merge(T2.os, idlist2, by="studyid")[, "new_ids"]
 T2.os$log_hr_se <- (T2.os$log_hr_u - T2.os$log_hr_l)/3.92
 
 T3.os$s_id3 <- merge(T3.os, idlist2, by="studyid")[, "new_ids"]
 N1.long <- dim(T1.os.long)[1]
 study1.long <- T1.os.long$s_id1
 x1.long <- T1.os.long$x
 censor.long <- T1.os.long$censor
 death.long <- T1.os.long$death
 t.long <- T1.os.long$time
 t1.long <- ceiling(T1.os.long$time)
 
 study_n1 <-length(unique(T1.os.long$s_id1))
 
 N2 <- dim(T2.os)[1]
 study2 <- T2.os$s_id2
 
 study_n2 <-length(unique(T2.os$s_id2))
 
 N3 <- dim(T3.os)[1]
 study3 <- T3.os$s_id3
 
 study_n3 <-length(unique(T3.os$s_id3))
 
 data.long <<- list(n_Study=n_Study, J=J, J_step=12, 
 study_n1=study_n1, N1=N1.long, study1=study1.long, x1=x1.long, time=T1.os.long$time, t1=t1.long, death=death.long, int1=ceiling(T1.os.long$t1), censor=censor.long, time_in_interval=T1.os.long$time_in_interval, 
 study_n2=study_n2, N2=N2, study2=study2, y2=T2.os$log_hr, s2=T2.os$log_hr_se, t2=T2.os$time, int2=ceiling(T2.os$time/12), 
 study_n3=study_n3, N3=N3, study3=study3, x3=T3.os$x, y3=T3.os$logit, s3=T3.os$logit_se, t3=T3.os$followup_time, int3=ceiling(T3.os$followup_time/12), event_prob=T3.os$event_prob, event_se=T3.os$se
 )
 T1.os.long<<-T1.os.long
 
 T1.os.reduced <- data.frame(T1.os.long)%>%
 group_by(t1, death, s_id1, x, time_in_interval)%>%
 summarise(count=n())
 
 N1.long <- dim(T1.os.reduced)[1]
 study1.long <- T1.os.reduced$s_id1
 x1.long <- T1.os.reduced$x
 death.long <- T1.os.reduced$death

 study_n1 <-length(unique(T1.os.reduced$s_id1))
 
 N2 <- dim(T2.os)[1]
 study2 <- T2.os$s_id2
 
 study_n2 <-length(unique(T2.os$s_id2))
 
 N3 <- dim(T3.os)[1]
 study3 <- T3.os$s_id3
 
 study_n3 <-length(unique(T3.os$s_id3))
 
 data.reduced <<- list(n_Study=n_Study, J=J, J_step=12, 
 study_n1=study_n1, N1=N1.long, study1=study1.long, x1=x1.long, t1=t1.long, death=death.long, int1=ceiling(T1.os.reduced$t1), time_in_interval=T1.os.reduced$time_in_interval, count=T1.os.reduced$count, 
 study_n2=study_n2, N2=N2, study2=study2, y2=T2.os$log_hr, s2=T2.os$log_hr_se, t2=T2.os$time, int2=ceiling(T2.os$time/12), 
 study_n3=study_n3, N3=N3, study3=study3, x3=T3.os$x, y3=T3.os$logit, s3=T3.os$logit_se, t3=T3.os$followup_time, int3=ceiling(T3.os$followup_time/12), event_prob=T3.os$event_prob, event_se=T3.os$se)
 
 Sys.setenv(Home="C:/Users/...")

 options(error=NULL)
 set_cmdstan_path("C:/Users/.../.cmdstan/cmdstan-2.32.2")
 mod_reduced<-cmdstan_model(stan_file_path, pedantic = F, compile=T);
 fit <<- mod_reduced$sample(
 data = data.reduced, 
 iter_warmup = 500, 
 iter_sampling = 250, 
 seed = 696999, 
 chains = 4, 
 parallel_chains = 4, 
 adapt_delta = 0.9, 
 refresh = 50,  
 max_treedepth=10
 )
 if(repz==1) orig_reduced<-data.frame(t(unlist(fit$summary(variables="HR_x1")[c(2)])))
 else orig_reduced<-rbind(orig_reduced, data.frame(t(unlist(fit$summary(variables="HR_x1")[c(2)]))))
 sim_b<-list(fit$summary(variables=c("HR_x1",   "log_lambda", "surv_rate_x0", "surv_rate_x1", "RMST_delta", "RMST_ratio", "RMST_x1", "RMST_x0")))
 sim_b$sim_num <- repz
 sim_b$cox_estimate_frequentist <-unlist(coxph(Surv(sim.tvh.t1$time, sim.tvh.t1$death) ~ sim.tvh.t1$x)[1])
 sim_summary<<-append(sim_summary, sim_b)
 saveRDS(sim_summary, here("simulation_output_tvh_20.rds"))


 for (i in 1:repz){
 if (i==1){
 mars_estimates<-list()
 cox_estimates<-list()
 mean<-list()
 surv_est_x0<-list()
 surv_est_x1<-list()
 
 }
 surv_est_x0<-rbind(surv_est_x0, unlist(sim_summary[[(i*3)-2]][[2]][57:66]))
 surv_est_x1<-rbind(surv_est_x1, unlist(sim_summary[[(i*3)-2]][[2]][67:76]))
 mars_estimates<-rbind(mars_estimates, unlist(sim_summary[[(i-1)*3+1]][[2]][1:10]))
 
 cox_estimates<-rbind(cox_estimates, unlist(sim_summary[[i*3]]))
 mean[i]<-mean(unlist(sim_summary[[(i*3)-2]][[2]][1:10]))[1]
 if (i==repz) {
 mars_estimates<-data.frame(matrix(unlist(mars_estimates), ncol=10), mean=unlist(mean), cox=unlist(cox_estimates), matrix(unlist(surv_est_x0), ncol=10), matrix(unlist(surv_est_x1), ncol=10))
 
 }
 }
 
 
 names(mars_estimates)[1:10] <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
 names(mars_estimates)[13:22] <- c("s0.01", "s0.02", "s0.03", "s0.04", "s0.05", "s0.06", "s0.07", "s0.08", "s0.09", "s0.10")
 names(mars_estimates)[23:32] <- c("s1.01", "s1.02", "s1.03", "s1.04", "s1.05", "s1.06", "s1.07", "s1.08", "s1.09", "s1.10")
 names(mars_estimates)
 ref<-list
 mars_est_long <- data.frame(name=unlist(pivot_longer(mars_estimates[, 1.:10], everything())["name"] ))
 
 mars_est_long$hr <- unlist(pivot_longer(mars_estimates[, 1.:10], everything())["value"] )
 mars_est_long$b <- log(mars_est_long$hr)
 mars_est_long$surv_est_x0 <- unlist(pivot_longer(mars_estimates[, 13.:22], everything())["value"] )
 mars_est_long$surv_est_x1 <- unlist(pivot_longer(mars_estimates[, 23.:32], everything())["value"] )
 
 mars_est_long$name<-as.numeric(mars_est_long$name)
 mars_est_long$x.start<-mars_est_long$name-.1
 mars_est_long$x.end<-mars_est_long$name+.1
 mars_est_long$ref<-exp(1-0.175^(mars_est_long$name*12/45-1))
 mars_est_long$y.start<-mars_est_long$ref-.01
 mars_est_long$y.end<-mars_est_long$ref+.01
 
 mars_est_long <- arrange(mars_est_long, name)
 mars_est_long2 <<-mars_est_long
 #boxplot( value ~ name, data = mars_est_long, ylim = c(0, 1), yaxs = "i", yaxt="n");axis(2, at = c(0, .25, .50, 0.75, 1)); 
 if(length(dev.list())!=0) dev.off()
 
 p<-ggplot(data=mars_est_long2, aes(x=name, y=hr, group=name))+
 geom_boxplot()+
 geom_segment(aes(x=x.start, xend=x.end, y=ref, yend=ref, color="red"))+
 scale_y_continuous(trans="log", breaks=c(.1, .2, .4, 1, 2))+
 theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
 panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
 
 
 print(p)
 }
 time_end<<-timestamp()
 return(sim_summary) 
}



sims<-simulation(repeats=1000)


sim_b<-list(fit$summary(variables=c("b",   "log_lambda", "surv_rate_x0", "surv_rate_x1", "RMST_delta", "RMST_ratio", "RMST_x1", "RMST_x0")))
sims<-sim_summary


for (i in 1:(length(sims)/3)){
 if (i==1){
 mars_estimates<-list()
 cox_estimates<-list()
 mean<-list()
 surv_est_x0<-list()
 surv_est_x1<-list()
 mars_estimates_randa<-list()
 mars_estimates_randb<-list()
 }
 surv_est_x0<-rbind(surv_est_x0, unlist(sims[[(i*3)-2]][[2]][21:30]))
 surv_est_x1<-rbind(surv_est_x1, unlist(sims[[(i*3)-2]][[2]][31:40]))
 mars_estimates<-rbind(mars_estimates, unlist(sims[[(i-1)*3+1]][[2]][1:10]))
 cox_estimates<-rbind(cox_estimates, unlist(sims[[i*3]]))
 mean[i]<-mean(unlist(sims[[(i*3)-2]][[2]][1:10]))[1]
 if (i == (length(sims)/3)) {
 mars_estimates<-data.frame(matrix(unlist(mars_estimates), ncol=10), 
                            mean=unlist(mean), 
                            cox=unlist(cox_estimates), 
                            matrix(unlist(surv_est_x0), ncol=10), 
                            matrix(unlist(surv_est_x1), ncol=10))
 

 }
}


names(mars_estimates)[1:10] <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
names(mars_estimates)[13:22] <- c("s0.01", "s0.02", "s0.03", "s0.04", "s0.05", "s0.06", "s0.07", "s0.08", "s0.09", "s0.10")
names(mars_estimates)[23:32] <- c("s1.01", "s1.02", "s1.03", "s1.04", "s1.05", "s1.06", "s1.07", "s1.08", "s1.09", "s1.10")

names(mars_estimates)
ref<-list
mars_est_long <- data.frame(name=unlist(pivot_longer(mars_estimates[, 1.:10], everything())["name"] ))

mars_est_long$hr <- unlist(pivot_longer(mars_estimates[, 1.:10], everything())["value"] )
mars_est_long$b <- log(mars_est_long$hr)
mars_est_long$surv_est_x0 <- unlist(pivot_longer(mars_estimates[, 13.:22], everything())["value"] )
mars_est_long$surv_est_x1 <- unlist(pivot_longer(mars_estimates[, 23.:32], everything())["value"] )

mars_est_long$name<-as.numeric(mars_est_long$name)
mars_est_long$x.start<-mars_est_long$name-.1
mars_est_long$x.end<-mars_est_long$name+.1
mars_est_long$ref<-exp(1-0.175^(mars_est_long$name*12/45-1))
mars_est_long$y.start<-mars_est_long$ref-.01
mars_est_long$y.end<-mars_est_long$ref+.01

mars_est_long <- arrange(mars_est_long, name)
#boxplot( value ~ name, data = mars_est_long, ylim = c(0, 1), yaxs = "i", yaxt="n");axis(2, at = c(0, .25, .50, 0.75, 1)); 
dev.off()
ggplot(data=mars_est_long, aes(x=name, y=hr, group=name))+
 geom_boxplot()+
 geom_segment(aes(x=x.start, xend=x.end, y=ref, yend=ref, color="red"))+
 scale_y_continuous(trans="log", breaks=c(.1, .2, .4, 1, 2))+
 theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
 panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

mean(mars_est_long$hr)
survival_probabilities<-data.frame(mars_est_long)%>%
 group_by(name)%>%
 summarise( s_x0=mean(surv_est_x0), sd_x0=sd(surv_est_x0), s_x1=mean(surv_est_x1), sd_x1=sd(surv_est_x1), .groups="keep")


