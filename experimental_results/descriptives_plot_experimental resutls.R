library(readxl)
library(dplyr)
library(ggplot2)
library(reshape2)
library(stringr)
library(writexl)

ser <- as.data.frame(read_excel("experimental_results/simulated_experiments_results.xlsx"))

ser = ser[ser$d<30,]

ser = ser[!duplicated(ser), ]

cases = ser[,c('obs','d','dof')]

cases_unique = cases[!duplicated(cases), ]
nrow(cases_unique)

obs_list <- c(5,10,30,50,100)
d_list <- c(2,5,10,30,50)
dof_list <- c(3,6,11,31,51,101)
df_all = data.frame()
for (obs in obs_list){
  for (d in d_list[d_list<=obs]){
    for (dof in dof_list[dof_list>d]){
      df_all <- rbind(df_all,data.frame(obs = obs, d = d, dof = dof))
      
    }
  }
}

dplyr::anti_join(df_all , cases_unique)



#ser = ser[ser$obs!=ser$d,]

#ser <- ser[ser$dof!=101,]

#ser$d = paste0('d = ',ser$d)
ser$dof = paste0('dof = ',ser$dof)

counts_rep = as.data.frame(ser %>% group_by(obs, d,dof, method, prior) %>%summarise(n=n()))
mean(counts_rep$n ==10)

ser[ser$time==60,]
ser = ser[complete.cases(ser$PE_dof), ]


ser_df_grouped <- as.data.frame(ser %>% group_by(obs, d,dof, method, prior) %>%summarise(max(time)))
ser_df_grouped[ser_df_grouped$method == 'stan_all_model',][250:350,]

ser['PE_ALL'] = ser$PE_dof + ser$PE_V

##################################################################
#####Tables in the paper####

ser_df_method_prior_grouped <- as.data.frame(ser[,c(4:ncol(ser))] %>% 
                                               group_by(method, prior) %>% 
                                               summarise(across(.cols = is.numeric, .fns = list(Mean = mean, SD = sd), 
                                                                  na.rm = TRUE, .names = "{col}_{fn}")))

ser_df_method_prior_grouped$method <- factor(ser_df_method_prior_grouped$method, 
                                             levels = c("Newton_Raphson","MwG_RwM","MwG_Slice","MwG_HMC","stan_all_model"))


ser_df_method_prior_grouped$prior <- factor(ser_df_method_prior_grouped$prior, 
                                             levels = c("uniform","exponential","gamma","inverse_gamma","log_normal"))

sorted_df_report <- ser_df_method_prior_grouped[
  order( ser_df_method_prior_grouped[,1], ser_df_method_prior_grouped[,2] ),
c("method","prior","time_Mean","time_SD","PE_dof_Mean","PE_dof_SD","PE_V_Mean","PE_V_SD", 
  "PE_ALL_Mean", "PE_ALL_SD","ESS_Mean","ESS_SD","Geweke_Mean","Geweke_SD")]

sorted_df_report

sorted_df_report = data.frame(lapply(sorted_df_report, function(y) if(is.numeric(y)) round(y, 3) else y)) 
write_xlsx(sorted_df_report,"experimental_results/sorted_df_report_low_dim.xlsx")

sorted_df_report[order(sorted_df_report$Geweke_SD),]





##################################################################
#####plot####
plot_ser_df <- ser[,c(1:5,7)]
#plot_ser_df <-as.data.frame(ser %>% group_by(obs,d,dof,method, prior) %>%summarise_all(mean))[,c(1:5,7)]

plot_ser_df['method_prior'] = paste0(plot_ser_df$method,'_',plot_ser_df$prior)

plot_ser_df = plot_ser_df[plot_ser_df$method_prior %in% c('MwG_HMC_log_normal',
                                                          'MwG_RwM_inverse_gamma',
                                                          'MwG_Slice_log_normal',
                                                          'Newton_Raphson_NA',
                                                          'stan_all_model_log_normal'),]


plot_ser_df$obs <- factor(plot_ser_df$obs, 
                          levels = c(5,10,30,50,100))

plot_ser_df$d <- factor(plot_ser_df$d, 
                        levels = c(2,5,10,30,50))

plot_ser_df$dof <- factor(plot_ser_df$dof, 
                          levels = c('dof = 3','dof = 6','dof = 11','dof = 31','dof = 51','dof = 101'))

plot_ser_df$method_prior <- factor(plot_ser_df$method_prior, 
                             levels = c("Newton_Raphson_NA",
                                        "MwG_RwM_inverse_gamma",
                                        "MwG_Slice_log_normal",
                                        "MwG_HMC_log_normal",
                                        "stan_all_model_log_normal"))

levels(plot_ser_df$method_prior) <- c('ML', 'RWM_Gibbs_inverse_gamma',
                                      'Slice_Gibbs_log_normal',
                                      'HMC_Gibbs_log_normal',
                                      "Stan_log_normal")


  

plot_title = paste0("PE of estimated degrees of freedom by \n method,  and dimension")

plot_dof = ggplot(plot_ser_df[plot_ser_df$obs==100,],
                  aes(x = d, y = PE_dof, fill = method_prior)) +
  geom_boxplot() +
  facet_wrap(~ dof, scales = "free_y", ncol = 6)+
  scale_y_continuous(name = "PE") +
  scale_x_discrete(labels = abbreviate, name = "dimensions")+
  labs(title=plot_title) + 
  theme(legend.spacing.y = unit(0.5, 'cm'),
        plot.title = element_text(face = "bold",size = 20,hjust = 0.5),
        strip.text = element_text(size = 15),
        axis.title=element_text(size=15,face="bold"),
        axis.text = element_text(size = 10),
        legend.text = element_text(size=15),
        legend.title = element_text(size=15,face="bold"))  +
  ## important additional element
  guides(fill = guide_legend(byrow = TRUE,title="Methods"))
  #scale_fill_brewer(palette="Set2",labels=c('Bisction', 'Newton-Raphson', 'Root \nSimulated Annealing',
  #                                          'Likelihood \nSimulated Annealing','Metropolis-within-Gibbs'))

plot_dof
#plot_name = paste0("plots/dof_boxplot.jpg")

jpeg(plot_name,width=3920, height=2000, res=300)
print(plot_dof)
dev.off()
  


for (obs_i in unique(plot_ser_df$obs)){
  
  plot_ser_df_obs = plot_ser_df[plot_ser_df$obs==obs_i,]
  
  plot_title = paste0("PE of estimated degrees of freedom by method and dimension\n for ",obs_i," observations")
  
  plot_dof = ggplot(plot_ser_df_obs,
                    aes(x = d, y = PE_dof, fill = method_prior)) +
    geom_boxplot() +
    facet_wrap(~ dof, scales = "free_y", ncol = 6)+
    scale_y_continuous(name = "PE") +
    scale_x_discrete(labels = abbreviate, name = "dimensions")+
    labs(title=plot_title) + 
    theme(legend.spacing.y = unit(0.5, 'cm'),
          plot.title = element_text(face = "bold",size = 20,hjust = 0.5),
          strip.text = element_text(size = 15),
          axis.title=element_text(size=15,face="bold"),
          axis.text = element_text(size = 10),
          legend.text = element_text(size=15),
          legend.title = element_text(size=15,face="bold"))  +
    guides(fill = guide_legend(byrow = TRUE,title="Methods"))+
    scale_fill_brewer(palette="Set2",labels=c('MLE', 'RWM within Gibbs \ninverse_gamma', 'SlS within Gibbs \nlog_normal',
                                              'HMC within Gibbs \nlog_normal','NUTS (Stan)\nlog_normal'))
  
  
  
  plot_name = paste0("plots/dof_boxplot_",obs_i,".jpg")
  
  jpeg(plot_name,width=3920, height=2000, res=300)
  print(plot_dof)
  dev.off()
  
}






