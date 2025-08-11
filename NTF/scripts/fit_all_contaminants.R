library(tidyverse)
library(parallel)
library(BNPdensity)

##Data_load
source('scripts/load_rivm_db.R')

source('scripts/SSD_fit_functions.R')

source('scripts/fit_analysis_functions.R')

fit_and_save_one = function(c_name = 'Atrazine', cens = F, Nit = 10000, 
                            save_ = T, 
                            jittered = T, 
                            factor = 1, 
                            geomean = F, 
                            filt_hickey = T, fun){
  
  fit = c_name %>% 
    get_log_dat(cens = cens, 
                jittered = jittered, 
                factor = factor, 
                centred_scaled = F, 
                geomean = T, 
                filt_hickey = filt_hickey) %>% 
    fun(Nit = 10)
  
  censtxt = ifelse(cens, 'c','nc')
  
  fname = paste('NTF/saves/posterior_samples_rivm/', c_name, censtxt,'_',fit$qmethod, sep='') %>% 
    gsub(".","_",., fixed = T) %>% 
    paste(.,'.Rdata',sep='')
  
  dir.create(dirname(fname), recursive = TRUE, showWarnings = FALSE)


  if(file.exists(fname)){
    if(file.info(fname)$size>5*10**5){
      print(paste(fname, 'already computed'))
      return()
    }
    else{
      fit = fit$data_ %>% 
        fun(Nit = Nit)
      
      if(save_) {saveRDS(fit, file = fname)
        print(paste('saved', fname))}
    }
  }
  
  else{
    fit = fit$data_ %>% 
      fun(Nit = Nit)
    
    if(save_) {saveRDS(fit, file = fname)
      print(paste('saved', fname))}
  }
  
}

set.seed( as.integer((as.double(Sys.time())*1000+Sys.getpid()) %% 2^31) )

c_name_list = rivm_db_hickey_filter %>% 
  as.data.frame() %>% 
  dplyr::arrange(desc(n_noncensored_species_geomean)) %>% 
  .$CAS %>% 
  unique() %>% 
  as.character() %>% 
  .[1:179] # To have more than 8 species


funlist = list(function(dat, Nit) fit_BNP_mixture_ssd_on_log_data_uniform(ddat = dat, mu.pz0 = 0.1, sigma.pz0 = 1.5, Nit = Nit))

expand.grid(c_name_list, funlist, c(T), stringsAsFactors = F) %>% # Comment when testing
  # expand.grid(c_name_list[1:3], funlist, c(T), stringsAsFactors = F) %>% # For testing purposes, use only the first 3 contaminants
  sample_frac(size = 1) %>% #This to make sure that when running the program again, the chains are not started in the same order. Else, chains that get stuck might cripple the whole process
  (function(x){
    mcmapply(FUN = function(cid, fun, cens_){
      fit_and_save_one(c_name = cid, cens = cens_, save_ = T,
                       Nit = 12000, jittered = T, factor = 0.5, geomean = T, filt_hickey = T, fun = fun
      )
    }, x$Var1, x$Var2, x$Var3, mc.preschedule = F, mc.cores = detectCores())
  })


