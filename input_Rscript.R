#INPUT:
##argument 1: row number of simulated file.
##argument 2: simulated rates file. 
###The code is sensitive to the column names of simulated rates file, please use the exactly same column #names: 
##num_tree: A integer used for numbering trees. It is used to call the phylogenetic tree, thus, it needs to ##be modified according to the tree name of the input trees.
##model_for_simulation: the true model that the input phylogenetic tree simulated with.
##model_for_estimation: the model that will be applied in the analysis.
##The following are the columns that include the initial parameters estimated by the true rates. It can ##also be modified to fixed initial parameters or optimized parameters.
#####true_lamb_par1, true_lamb_linear_par1, true_lamb_linear_par2, true_lamb_expo_par1, #true_lamb_expo_par2, true_mu_par1, true_mu_linear_par1, true_mu_linear_par2, #true_mu_expo_par1, true_mu_expo_par2
##The following are the columns that include simulated rates.
#####simulated_initial_speciation_rate, simulated_initial_extinction_rate, #simulated_initial_diversification_rate, simulated_final_speciation_rate, #simulated_final_extinction_rate, simulated_final_diversification_rate

#OUTPUT:
##The first five columns include the output of fit_bd function.
##The other columns include estimated initial (past) and final (present) speciation, extinction and diversification rates as well as the deviation and relative deviation from the simulated rates. fit_bd function.


##################***************RSCRIPT**************#######################
############# This script is used for estimation by using crown age.####################
#The R script is used for estimating initial and final speciation,extinction and diversification rates 
#as well as calculate the deviation and relative deviation.
#argument 1: row number
#argument 2: simulated rates file (sensitive to the colnames)
setwd("/home/u19/cdi/RPANDA_true_initial_run1/ccc")
#parameters1:input tree numbers, start form a ends at b
rm(list=ls())
# install packages. Remember to install the packages 
# on server-local place by running R script directly.
#install.packages("permute")
#install.packages("lattice")
#install.packages("vegan")
#install.packages("methods")#problem
#install.packages("nlme")
#install.packages("mgcv")
#install.packages("expm")
#install.packages("gee")
#install.packages("ape")
#install.packages("picante")
#install.packages("RPANDA")
#install.packages("pbapply")
#install.packages("Matrix")
library(permute)
library(lattice)
library(vegan)
library(methods)
library(Matrix)
library(nlme)
library(mgcv)
library(expm)
library(gee)
library(ape)
library(picante)
library(RPANDA)
library(pbapply)
#parameter2:the right model and model for simulation,cl+ll+el+ce+le+ee+cc+lc+ec
cc=1
lc=2
ec=3
cl=4
ll=5
el=6
ce=7
le=8
ee=9
#############################################################
#Need to be change manually for different sampling size. "fra" is the fraction of sampling species
##of total species. For example, fra=0.25 means there are 25% of species been sampled.
fra=1
args = commandArgs(trailingOnly=TRUE)
#*****need to be change for diffrent models.
#read trees according to the input row number and "num_tree" in input file.
#args[1]<-"1"
#args[2]<-"out.csv"
rownumber<-args[1]
rownumber<-as.numeric(rownumber)
simulate<-args[2]
simulate<-toString(simulate)
simulate<-read.csv(simulate)

#rownumber<-toString(rownumber)
simodel<-simulate$model_for_simulation[rownumber]
simodel<-toString(simodel)
rightmodel=get(simodel)

esmodel<-simulate$model_for_estimation[rownumber]
esmodel<-toString(esmodel)
r_esmodel<-esmodel
esmodel=get(esmodel)
#pay attention, still need to change some parameters in the RPANDA package
###################################################################
modelname<-c("cc","lc","ec","cl","ll","el","ce","le","ee")
#read arguments-inputs
treenumber<-simulate$num_tree[rownumber]
treenumber<-toString(treenumber)
intree<-paste0("tree_",simodel,"_",treenumber,".tre")
intree<-read.tree(intree)

### initial set up for RPANDA_function for speciation and extinction rates-------------------
#set intitial subobjects to null
s_tree<-NULL
#set initial parameters for RPANDA
#mu
f.mu_c<-function(t,y){y[1]}
f.mu_l<-function(t,y){y[1]+y[2]*t} 
f.mu_expo<-function(t,y){y[1]*exp(y[2]*t)}
#lamb
f.lamb_c<-function(t,y){y[1]}
f.lamb_l<-function(t,y){y[1]+y[2]*t} 
f.lamb_expo<-function(t,y){y[1]*exp(y[2]*t)}
#initial
lamb_par_init_c<-c(simulate$true_lamb_par1[rownumber])
mu_par_init_c<-c(simulate$true_mu_par1[rownumber])

lamb_par_init_linear<-c(simulate$true_lamb_par1[rownumber],simulate$true_lamb_linear_par2[rownumber])
mu_par_init_linear<-c(simulate$true_mu_par1[rownumber],simulate$true_mu_linear_par2[rownumber])

lamb_par_init_expo<-c(simulate$true_lamb_par1[rownumber],simulate$true_lamb_expo_par2[rownumber])
mu_par_init_expo<-c(simulate$true_mu_par1[rownumber],simulate$true_mu_expo_par2[rownumber])

#tot_time_crown
tot_time<-max(node.age(intree)$ages)
res_crown<-NULL
#records for initial rates
i_speci_crown<-NA
i_extin_crown<-NA
i_netdiv_crown<-NA
i_deviation_of_speciation_rate<-NA
i_deviation_of_extinction_rate<-NA
i_deviation_of_diversification_rate<-NA
i_relative_deviation_of_speciation_rate<-NA
i_relative_deviation_of_extinction_rate<-NA
i_relative_deviation_of_diversification_rate<-NA
#records for final rates
f_speci_crown<-NA
f_extin_crown<-NA
f_netdiv_crown<-NA
f_deviation_of_speciation_rate<-NA
f_deviation_of_extinction_rate<-NA
f_deviation_of_diversification_rate<-NA
f_relative_deviation_of_speciation_rate<-NA
f_relative_deviation_of_extinction_rate<-NA
f_relative_deviation_of_diversification_rate<-NA
#records of RPANDA parameters
lamb_par1<-NA
lamb_par2<-NA
mu_par1<-NA
mu_par2<-NA
tryCatch(
  { 
    ##############################################################################
    #######*******NEED TO BE CHANGE FOR DIFFERENT ESTIMATED MODEL
    if (esmodel==cc){
      #1 cc
      tryCatch(
        { res_crown<-fit_bd(intree, tot_time, f.lamb_c, f.mu_c, 
                            lamb_par_init_c, mu_par_init_c, f =fra,
                            # meth = "Nelder-Mead",
                            cst.lamb = TRUE, cst.mu = T,
                            expo.lamb = FALSE, expo.mu = FALSE, fix.mu = FALSE,
                            dt=0, cond = "crown")
        lamb_par1<-res_crown$lamb_par[1]
        mu_par1<-res_crown$mu_par[1]},
        error=function(e){})
      print("cc")} else if(esmodel==lc ){
        #2 lc
        tryCatch(
          { res_crown<-fit_bd(intree, tot_time, f.lamb_l, f.mu_c,
                              lamb_par_init_linear, mu_par_init_c, f =fra,
                              # meth = "Nelder-Mead",
                              cst.lamb = F, cst.mu = T,
                              expo.lamb = F, expo.mu = F, fix.mu = FALSE,
                              dt=0, cond = "crown")
          lamb_par1<-res_crown$lamb_par[1]
          mu_par1<-res_crown$mu_par[1]
          lamb_par2<-res_crown$lamb_par[2]},
          error=function(e){})
        print("lc")} else if (esmodel==ec){
          #3 ec
          tryCatch(
            { res_crown<-fit_bd(intree, tot_time, f.lamb_expo, f.mu_c,
                                lamb_par_init_expo, mu_par_init_c, f =fra,
                                # meth = "Nelder-Mead",
                                cst.lamb = F, cst.mu = T,
                                expo.lamb = T, expo.mu = F, fix.mu = FALSE,
                                dt=0, cond = "crown")
            lamb_par1<-res_crown$lamb_par[1]
            mu_par1<-res_crown$mu_par[1]
            lamb_par2<-res_crown$lamb_par[2]  },
            error=function(e){})
          print("3")} else if (esmodel==cl){
            #4 cl
            tryCatch(
              { res_crown<-fit_bd(intree, tot_time, f.lamb_c, f.mu_l, 
                                  lamb_par_init_c,mu_par_init_linear, f =fra,                     
                                  cst.lamb = TRUE, cst.mu = FALSE,
                                  expo.lamb = FALSE, expo.mu = FALSE, fix.mu = FALSE,
                                  dt=0, cond = "crown")
              lamb_par1<-res_crown$lamb_par[1]
              mu_par1<-res_crown$mu_par[1]
             # lamb_par2<-res_crown$lamb_par[2]
              mu_par2<-res_crown$mu_par[2]},
              error=function(e){})
            print("cl")} else if (esmodel==ll){
              #5 ll
              tryCatch(
                { res_crown<-fit_bd(intree, tot_time, f.lamb_l, f.mu_l, 
                                    lamb_par_init_linear,mu_par_init_linear, f =fra,
                                    # meth = "Nelder-Mead",
                                    cst.lamb = F, cst.mu = F,
                                    expo.lamb = F, expo.mu = F, fix.mu = FALSE,
                                    dt=0, cond = "crown")
                lamb_par1<-res_crown$lamb_par[1]
                mu_par1<-res_crown$mu_par[1]
                lamb_par2<-res_crown$lamb_par[2]
                mu_par2<-res_crown$mu_par[2]  },
                error=function(e){})
              print("ll")} else if (esmodel==el){
                #6 el
                tryCatch(
                  { res_crown<-fit_bd(intree, tot_time, f.lamb_expo, f.mu_l, 
                                      lamb_par_init_expo, mu_par_init_linear, f =fra,
                                      # meth = "Nelder-Mead",
                                      cst.lamb = F, cst.mu = F,
                                      expo.lamb = T, expo.mu = F, fix.mu = FALSE,
                                      dt=0, cond = "crown")
                  lamb_par1<-res_crown$lamb_par[1]
                  mu_par1<-res_crown$mu_par[1]
                  lamb_par2<-res_crown$lamb_par[2]
                  mu_par2<-res_crown$mu_par[2]  },
                  error=function(e){})
                print("el")} else if (esmodel==ce){
                  #7 ce
                  tryCatch(
                    { res_crown<-fit_bd(intree, tot_time, f.lamb_c, f.mu_expo, 
                                        lamb_par_init_c, mu_par_init_expo, f =fra,
                                        # meth = "Nelder-Mead",
                                        cst.lamb = TRUE, cst.mu = FALSE,
                                        expo.lamb = FALSE, expo.mu = T, fix.mu = FALSE,
                                        dt=0, cond = "crown")
                    lamb_par1<-res_crown$lamb_par[1]
                    mu_par1<-res_crown$mu_par[1]
                    # lamb_par2<-res_crown$lamb_par[2]
                    mu_par2<-res_crown$mu_par[2] },
                    error=function(e){})
                  print("ce")} else if (esmodel==le){
                    #8 le
                    tryCatch(
                      { res_crown<-fit_bd(intree, tot_time, f.lamb_l, f.mu_expo, 
                                          lamb_par_init_linear,mu_par_init_expo, f =fra,
                                          # meth = "Nelder-Mead",
                                          cst.lamb = F, cst.mu = F,
                                          expo.lamb = F, expo.mu = T, fix.mu = FALSE,
                                          dt=0, cond = "crown")
                      lamb_par1<-res_crown$lamb_par[1]
                      mu_par1<-res_crown$mu_par[1]
                      lamb_par2<-res_crown$lamb_par[2]
                      mu_par2<-res_crown$mu_par[2]  },
                      error=function(e){})
                    print("le")} else if (esmodel==ee){
                      #9 ee
                      tryCatch(
                        { res_crown<-fit_bd(intree, tot_time, f.lamb_expo, f.mu_expo, 
                                            lamb_par_init_expo, mu_par_init_expo, f =fra,
                                            # meth = "Nelder-Mead",
                                            cst.lamb = F, cst.mu = F,
                                            expo.lamb = T, expo.mu = T, fix.mu = FALSE,
                                            dt=0, cond = "crown")
                        lamb_par1<-res_crown$lamb_par[1]
                        mu_par1<-res_crown$mu_par[1]
                        lamb_par2<-res_crown$lamb_par[2]
                        mu_par2<-res_crown$mu_par[2]  },
                        error=function(e){})
                      print("ee")} else {
                        print("wrong input estimated model")
                      }
    ###############################################
    #######################################

    #initial
    inirates<-function(model,y1,y2,t){
      if (model=="c"){
        rate=y1
      }else if(model=="l"){
        rate=y1+y2*t
      }else if(model=="e"){
        rate=y1*exp(y2*t)
      }
      rate<-abs(rate)
      return(rate)
    }
    speciation_model<-substring(r_esmodel,1,1)
    i_speci_crown=
      inirates(speciation_model,lamb_par1,lamb_par2,tot_time)
    
    extinction_model<-substring(r_esmodel,2,2)
    i_extin_crown=
      inirates(extinction_model,mu_par1,mu_par2,tot_time)
    i_netdiv_crown<- i_speci_crown -  i_extin_crown
    
    i_deviation_of_speciation_rate<-i_speci_crown- simulate$simulated_initial_speciation_rate[rownumber]
    i_deviation_of_extinction_rate<- i_extin_crown-simulate$simulated_initial_extinction_rate[rownumber]
    i_deviation_of_diversification_rate<-i_netdiv_crown- simulate$simulated_initial_diversification_rate[rownumber]
    
    i_relative_deviation_of_speciation_rate<-100*(i_deviation_of_speciation_rate/(simulate$simulated_initial_speciation_rate[rownumber]))
    i_relative_deviation_of_extinction_rate<-100*(i_deviation_of_extinction_rate/(simulate$simulated_initial_extinction_rate[rownumber]))
    i_relative_deviation_of_diversification_rate<-100*(i_deviation_of_diversification_rate/(simulate$simulated_initial_diversification_rate[rownumber]))
    
    #final Use absolute values
    f_speci_crown<-abs(lamb_par1)
    f_extin_crown<-abs(mu_par1)
    f_netdiv_crown<- f_speci_crown -  f_extin_crown
    
    f_deviation_of_speciation_rate<-f_speci_crown- simulate$simulated_final_speciation_rate[rownumber]
    f_deviation_of_extinction_rate<- f_extin_crown-simulate$simulated_final_extinction_rate[rownumber]
    f_deviation_of_diversification_rate<-f_netdiv_crown- simulate$simulated_final_diversification_rate[rownumber]
    
    f_relative_deviation_of_speciation_rate<-100*(f_deviation_of_speciation_rate/(simulate$simulated_final_speciation_rate[rownumber]))
    f_relative_deviation_of_extinction_rate<-100*(f_deviation_of_extinction_rate/(simulate$simulated_final_extinction_rate[rownumber]))
    f_relative_deviation_of_diversification_rate<-100*(f_deviation_of_diversification_rate/(simulate$simulated_final_diversification_rate[rownumber]))
    s_tree<-c(
              res_crown$LH,
              res_crown$aicc,
              lamb_par1,
              lamb_par2,
              mu_par1,
              mu_par2,
              #1. initial rates
              i_speci_crown,
              i_extin_crown,
              i_netdiv_crown,

              
              i_deviation_of_speciation_rate,
              i_deviation_of_extinction_rate,
              i_deviation_of_diversification_rate,
              
              i_relative_deviation_of_speciation_rate,
              i_relative_deviation_of_extinction_rate,
              i_relative_deviation_of_diversification_rate,
              
              
              #2. final rates
              f_speci_crown,
              f_extin_crown,
              f_netdiv_crown,
              
              f_deviation_of_speciation_rate,
              f_deviation_of_extinction_rate,
              f_deviation_of_diversification_rate,
              
              f_relative_deviation_of_speciation_rate,
              f_relative_deviation_of_extinction_rate,
              f_relative_deviation_of_diversification_rate)
  },
  error=function(e){})

#write 9 models into sum_s 
s_tree<-rbind(s_tree)
colnames(s_tree)<-c(
                    "ti_maximum_log_likelihood_value",
                    "ti_aicc",
                    "ti_lamb_par1",
                    "ti_lamb_par2",
                    "ti_mu_par1",
                    "ti_mu_par2",
                    
                    
                    
                    "ti_estimated_initial_speciation_rate",
                    "ti_estimated_initial_extinction_rate",
                    "ti_estimated_initial_diversification_rate",

                    
                    "ti_deviation_of_initial_speciation_rate",
                    "ti_deviation_of_initial_extinction_rate",
                    "ti_deviation_of_initial_diversification_rate",
                    
                    "ti_relative_deviation_of_initial_speciation_rate",
                    "ti_relative_deviation_of_initial_extinction_rate",
                    "ti_relative_deviation_of_initial_diversification_rate",
                    
                    "ti_estimated_final_speciation_rate",
                    "ti_estimated_final_extinction_rate",
                    "ti_estimated_final_diversification_rate",

                    "ti_deviation_of_final_speciation_rate",
                    "ti_deviation_of_final_extinction_rate",
                    "ti_deviation_of_final_diversification_rate",
                    
                    "ti_relative_deviation_of_final_speciation_rate",
                    "ti_relative_deviation_of_final_extinction_rate",
                    "ti_relative_deviation_of_final_diversification_rate")

output<-cbind(simulate[rownumber,],s_tree)
write.csv(output,paste0("true_initial_",treenumber,"_",modelname[rightmodel],"_",modelname[esmodel],".csv"),row.names = F)


