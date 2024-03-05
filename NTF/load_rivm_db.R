library(tidyverse)

# ex = data.frame(conc.low = c(0.5,1,6,5,4),
#                 conc.upp = c(3,1,7,5.5,6),
#                 species = c(1,1,2,2,2),
#                 CAS = c(111,111,111,111,111))

geomean = function(x, na.rm = FALSE){
  if(any(x < 0, na.rm = TRUE)){
    return(NaN)
  }
  if(any(x == 0, na.rm = TRUE)){
    return(0)
  }
  else exp(mean(log(x), na.rm = na.rm))
}


deal_with_duplicates = function(db_grouped_by_CAS){
  db_grouped_by_CAS %>%
    mutate(conc.point = 0.5*(conc.low + conc.upp)) %>%
    group_by(CAS, chem.name, chem.name2, species) %>%
    summarise(conc.low = min(conc.low, na.rm = T),
              conc.upp = max(conc.upp, na.rm = T),
              conc.geomean = geomean(conc.point, na.rm = TRUE))
}



filter_according_to_Hickey2012 = function(rivm_db){
  ##############
  ## PREAMBLE ##
  ##############
  ## This code is produced by Peter Craig, Dept. of Mathematical Sciences,
  ## Durham University and Graeme Hickey, NIBHI, Manchester University.
  ## It is freely available for all purposes. However, the authors of the manuscript
  ## assume no responsibility for any possible errors in the code.
  ##
  ##
  # If you have any questions, please contact: P.S.Craig@durham.ac.uk.
  ##
  ## To being using the code, you will need to install the R statistical software
  ## program and the packages listed below.
  library(doBy)
  library(stringr)
  library(reshape)
  library(ggplot2)
  library(Matrix)
  # rivm = dget(file.choose()) # load a CSV version of the database from files
  rivm = rivm_db # I add this so it fits in the pipeline

  rivm = rivm[rivm$endpoint %in% c("EC50", "LC50", "NOEC"), ]
  rivm = droplevels(rivm)
  #----------------------------------------------------------------
  ######################
  ## DATA MANAGEMENT ##
  ######################

  # full robust acute (EC50/LC50)
  inds1 = with(rivm,
               (endpoint == "LC50" | endpoint == "EC50") &
                 (effect == "MOR" | effect == "IMM"))
  inds2 = with(rivm,
               (dur.low == 2) &
                 (major %in% c("CR", "IN")))
  inds3 = with(rivm,
               (dur.low == 4) &
                 !(major %in% c("CR", "IN")))
  inds4 = !(grepl(" sp$", rivm$species) | grepl(" sp.", rivm$species))
  inds5 = grepl(" ", rivm$species)
  inds6 = !(rivm$major == "MI")
  inds7 = !(rivm$conc.ind == "A")
  acute.r = rivm[which(inds1 & (inds2 | inds3) & inds4 & inds5 & inds6 & inds7), ]
  # acute.r = drop.levels(acute.r)
  acute.r = droplevels(acute.r)#Maybe there was a misspell

  #This is probably interesting but we are not making database-wide studies for the moment
  #   # Robust with n >= 5 distinct species *pointwise* measurements
  #   # (incl. censored data)
  #   n = by(
  #     acute.r,
  #     factor(acute.r$CAS),
  #     function(d) length(unique(d[d$conc.ind == "P", ]$species))
  #   )
  #   status = (n >= 5)
  #   acute.r2 = acute.r[acute.r$CAS %in% names(n)[status], ]
  #   acute.r2 = drop.levels(acute.r2)
  #   # Robust with n >= 5 distinct species *pointwise* measurements
  #   # (not incl. censored data)
  #   acute.r3 = acute.r2[acute.r2$conc.ind == "P", ]
  #   acute.r3 = drop.levels(acute.r3)

  return(acute.r)
}


Correct_some_stuff_in_the_database = function(db){
  db %>%
    mutate(chem.name=replace(chem.name, db$CAS==50066, '5-Ethyl-5-phenyl-2,4,6(1H,3H,5H)'),
           chem.name2=replace(chem.name2, db$CAS==50066, 'pyrimidinetrione')) %>%
    as.data.frame()
}


rivm_db_hickey_filter = get(load('rivm_db.Rdata')) %>%
  Correct_some_stuff_in_the_database %>%
  # subset(endpoint%in%c('LC50','EC50')) %>% #Selection some endpoints
  filter_according_to_Hickey2012 %>%
  deal_with_duplicates %>%
  group_by(CAS) %>% #Centring and scaling
  mutate(lconc.point = 0.5*(conc.low+conc.upp) %>% log10,
         lconc.geomean = log10(conc.geomean)) %>% #transforming interval censored data into pointwise data
  mutate(lconc.low = log10(conc.low), lconc.upp = log10(conc.upp) ) %>% #This was lost in the summarise operation
  mutate(centre = mean(lconc.point, na.rm = T), scale = sd(lconc.point, na.rm = T),
         centre_geo = mean(lconc.geomean, na.rm = T), scale_geo = sd(lconc.geomean, na.rm = T)) %>% #calculating the centre and sd of the distribution
  mutate(lconc.point.centred_scaled = (lconc.point - centre)/scale,
         lconc.geomean.centred_scaled = (lconc.geomean - centre_geo)/scale_geo,
         lconc.low.centred_scaled = (log10(conc.low) - centre)/scale,
         lconc.upp.centred_scaled = (log10(conc.upp) - centre)/scale) %>% #centring
  mutate(n_species = length(lconc.point),
         n_noncensored_species = lconc.point %>% na.omit %>% length,
         n_noncensored_species_geomean = lconc.geomean %>% na.omit %>% length) %>%
  mutate(shorter_name = ifelse(is.na(chem.name2), yes = chem.name %>% as.character(), no = chem.name2 %>% as.character()))


rivm_db = get(load('rivm_db.Rdata')) %>%
  Correct_some_stuff_in_the_database %>%
  subset(endpoint%in%c('LC50','EC50')) %>% #Selection some endpoints
  # filter_according_to_Hickey2012 %>%
  deal_with_duplicates %>%
  group_by(CAS) %>% #Centring and scaling
  mutate(lconc.point = 0.5*(conc.low+conc.upp) %>% log10,
         lconc.geomean = log10(conc.geomean)) %>% #transforming interval censored data into pointwise data
  mutate(lconc.low = log10(conc.low), lconc.upp = log10(conc.upp) ) %>% #This was lost in the summarise operation
  mutate(centre = mean(lconc.point, na.rm = T), scale = sd(lconc.point, na.rm = T),
         centre_geo = mean(lconc.geomean, na.rm = T), scale_geo = sd(lconc.geomean, na.rm = T)) %>% #calculating the centre and sd of the distribution
  mutate(lconc.point.centred_scaled = (lconc.point - centre)/scale,
         lconc.geomean.centred_scaled = (lconc.geomean - centre_geo)/scale_geo,
         lconc.low.centred_scaled = (log10(conc.low) - centre)/scale,
         lconc.upp.centred_scaled = (log10(conc.upp) - centre)/scale) %>% #centring
  mutate(n_species = length(lconc.point),
         n_noncensored_species = lconc.point %>% na.omit %>% length,
         n_noncensored_species_geomean = lconc.geomean %>% na.omit %>% length) %>%
  mutate(shorter_name = ifelse(is.na(chem.name2), yes = chem.name %>% as.character(), no = chem.name2 %>% as.character()))

create_CAS_shortname_converter = function(){
  cv = rivm_db %>%
    dplyr::select(CAS, shorter_name) %>%
    unique() %>%
    (function(x){
      c(x$CAS %>% as.character(),x$shorter_name) %>%
        setNames(c(x$shorter_name,x$CAS %>% as.character()))
    })

  function(x){
    x %>% sapply(function(y) cv[[y %>% as.character]])
  } %>% return

}

CAS_short_name_converter = create_CAS_shortname_converter()

is.censored = function(ddat){
  if(is.null(ncol(ddat))) FALSE
  else if(ncol(ddat)==1) FALSE
  else if(ncol(ddat)==2) TRUE
  else stop('Wrong type/dim of data input')
}

big_c_list = c('Atrazine','DDT','Endosulfan')
med_c_list = c('Zinc','Chlorine','Antimycin')
small_c_list = c('Ethanol','Leptophos','Hexane')

#For BAYSM
# big_c_list = c('Cadmium chloride','Zinc chloride','Lindane')
# med_c_list = c('Atrazine','Sodium cyanide','Antimycin')
# small_c_list = c('Captan','Propoxur','Alachlor')

#For the JASA paper, selected by Julyan
# big_c_list = c('122145','7733020','298000')
# med_c_list = c('7632000','654660','91203')
# small_c_list = c('9016459','61791262','1152028')

ex_c_list = c(big_c_list, med_c_list, small_c_list)

jitter_if_ = function(x, jittered, factor){
  x %>%
    (function(x){
      if(jittered) jitter(x, factor = factor)
      else x
    })
}

get_smallest_difference_censored_df = function(df){
  0.5*(df$left+df$right) %>%
    na.omit() %>%
    as.numeric() %>%
    unique() %>%
    combn(m = 2, simplify = F) %>% #produce all pairs
    lapply(FUN = function(pair) abs(pair[1]-pair[2])) %>%
    unlist() %>%
    min()
}

jitter_censored_df = function(df, factor){
  mindist = df %>% get_smallest_difference_censored_df

  rndjitter = runif(n = dim(df)[1], min = - mindist/5*factor, max = mindist/5*factor)

  df %>%
    mutate(left = left + rndjitter, right = right + rndjitter)
}

jitter_if_df = function(df, jittered, factor){
  df %>%
    (function(x){
      if(jittered) jitter_censored_df(x, factor = factor)
      else x
    })
}


# sort_after_jittering = function(xleft,xright){
#   if(is.na(xleft)|is.na(xright)) c(xleft,xright)
#   else c(xleft,xright) %>% sort
# }
#
# sort_df_after_jittering = function(df){
#   mapply(sort_after_jittering, df$left, df$right, SIMPLIFY = F) %>%
#     Reduce(rbind,.) %>%
#     data.frame() %>%
#     setNames(c('left','right'))
# }


add_name_if = function(df_, test, nm){
  df_ %>%
    (function(df){
      if(test) df %>% mutate(species = nm)
      else df
    })
}

as.numeric_keep_names = function(x){
  x %>%
    (function(y){
      y %>%
        as.numeric() %>%
        setNames(names(x))
    })
}

get_log_dat = function(contaminant, cens = FALSE,
                              centred_scaled = FALSE,
                              jittered = FALSE, factor = 1e-3,
                              geomean = F, filt_hickey = F, names = F){

    set.seed(0)

    if(filt_hickey) db = rivm_db_hickey_filter
    else db = rivm_db

      nm = db %>%
        subset(shorter_name == contaminant %>% as.character()| CAS == contaminant %>% as.character()) %>%
        .$species

    db %>%
      subset(shorter_name == contaminant %>% as.character()| CAS == contaminant %>% as.character()) %>%
      (function(x){
        if(cens){
          if(centred_scaled){
            data.frame(left = x$lconc.low.centred_scaled,
                       right = x$lconc.upp.centred_scaled ) %>%
              jitter_if_df(jittered, factor = factor) %>%
              add_name_if(names, nm)
          }
          else{
            data.frame(left = x$lconc.low,
                       right = x$lconc.upp) %>%
              jitter_if_df(jittered, factor = factor) %>%
              add_name_if(names, nm)
          }
        }
        else{
          if(centred_scaled){
            if(geomean){
              x$lconc.geomean.centred_scaled %>%
                setNames(nm) %>%
                na.omit %>%
                as.numeric_keep_names  %>%
                jitter_if_(jittered, factor)
            }
            else{
              x$lconc.point.centred_scaled %>%
                setNames(nm) %>%
                na.omit %>%
                as.numeric_keep_names  %>%
                jitter_if_(jittered, factor) }
          }
          else{
            if(geomean){
              x$lconc.geomean %>%
                setNames(nm) %>%
                na.omit %>%
                as.numeric_keep_names  %>%
                jitter_if_(jittered, factor)
            }
            else{
              x$lconc.point %>%
                setNames(nm) %>%
                na.omit %>%
                as.numeric_keep_names  %>%
                jitter_if_(jittered, factor)
            }
          }
        }
      })
}

get_log_dat_noname = function(contaminant, cens = FALSE,
                       centred_scaled = TRUE,
                       jittered = FALSE, factor = 1e-3,
                       geomean = F, filt_hickey = F){
  set.seed(0)

  if(filt_hickey) db = rivm_db_hickey_filter
  else db = rivm_db

  db %>%
    subset(shorter_name == contaminant %>% as.character()| CAS == contaminant %>% as.character()) %>%
    (function(x){
      if(cens){
        if(centred_scaled){
          data.frame(left = x$lconc.low.centred_scaled,
                     right = x$lconc.upp.centred_scaled ) %>%
            jitter_if_df(jittered, factor = factor)
        }

        else{
          data.frame(left = x$lconc.low,
                     right = x$lconc.upp) %>%
            jitter_if_df(jittered, factor = factor)
        }
      }
      else{
        if(centred_scaled){
          if(geomean){
            x$lconc.geomean.centred_scaled %>%
              na.omit %>%
              as.numeric  %>%
              jitter_if_(jittered, factor)
          }
          else{
            x$lconc.point.centred_scaled %>%
              na.omit %>%
              as.numeric  %>%
              jitter_if_(jittered, factor) }
        }
        else{
          if(geomean){
            x$lconc.geomean %>%
              na.omit %>%
              as.numeric  %>%
              jitter_if_(jittered, factor)
          }
          else{
            x$lconc.point %>%
              na.omit %>%
              as.numeric  %>%
              jitter_if_(jittered, factor)
          }
        }
      }
    })
}

saveRDS(rivm_db_hickey_filter, file = '~/Documents/Github/BNP_SSD/data/rivm_db_hickey_filter.RDS')


