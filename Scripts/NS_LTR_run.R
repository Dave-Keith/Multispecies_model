# add other packages here:
library(dplyr)
library(readr)
library(tibble)
library(csasdown)
# OK, so using the ICES assessments here's what we get for North Sea cod.
library(readxl)
library(tidyverse)
library(rio)
library(ggthemes)
library(cowplot)
library(gtable)
library(viridis)
library(arm)
library(mgcv)
library(lme4)
library(scales)
library(ggh4x)
library(gratia)
#library(legendary)
theme_set(theme_few(base_size = 14))


# The location of your data
loc <- 'D:/GitHub/Multispecies_model/'
# The parameters for the doubling time simulation
n.dt.sims <- 1000 # Increase to 1000 for final run, which takes like a full day or something 
n.sim.years <- 100

# What tunning model are you running or did you run!
what.2.tune <- 'fec_nm'

#source(paste0(loc,"/Scripts/functions/tuning_sim_fast.R"))

# Download the functions we'll need from github
funs <- c("https://raw.githubusercontent.com/dave-keith/ICM/main/Scripts/functions/simple_Lotka_r.r",
          "https://raw.githubusercontent.com/dave-keith/ICM/main/Scripts/functions/backwards_project.r",
          "https://raw.githubusercontent.com/dave-keith/ICM/main/Scripts/functions/forward_project.r",
          "https://raw.githubusercontent.com/dave-keith/ICM/main/Scripts/functions/tuning_sim_fast.r"
)
# Now run through a quick loop to load each one, just be sure that your working directory is read/write!
for(fun in funs) 
{
  download.file(fun,destfile = basename(fun))
  source(paste0(getwd(),"/",basename(fun)))
  file.remove(paste0(getwd(),"/",basename(fun)))
}

  ASR <- read_xlsx("D:/Github/Age-structured-marine-fish-database/MSM/MSM_data_NS_stocks.xlsx")
  
  datatypes <- unique(gsub(x = names(ASR), pattern = "[^a-zA-Z]", replacement=""))
  
  # we want:
  # Year, Num, WA, Catch, AM, NM, StockID, Management, Area, Order, Family, Genus, Species
  ASRdat <- ASR[,c(grep(x=names(ASR), "Num"),
                   grep(x=names(ASR), "WA"),
                   grep(x=names(ASR), "FM"),
                   grep(x=names(ASR), "AM"),
                   grep(x=names(ASR), "NM"))]
  # Making all the data numeric that should be numeric
  ASRdat <- apply(X = ASRdat, 2, as.numeric)
  # Getting the species info back
  ASRsp <- ASR[, which(!1:length(names(ASR)) %in% grep(x=names(ASR), ".", fixed=T))]
  # And binding it all back together
  ASR.trim <- cbind(ASRsp, ASRdat)
  # need a unique ID for stock
  #table(ASR.trim$Management, ASR.trim$Species)
  ASR.trim$Stock <- paste0(ASR.trim$Management, "_", ASR.trim$Area, "_", ASR.trim$Genus, "_", ASR.trim$Species)
  
  ASR.long <- ASR.trim %>%
    pivot_longer(!c("Management", "Area", "Order", "Family", "Genus", "Species", "Stock", "Year","Meeting_or_reference","Model","Case","Notes")) %>%
    separate(col=name, into=c("type", "age"), sep = "\\.")
  
  Stocks <- ASR.long %>%
                      dplyr::filter(!is.na(value)) %>%
                      dplyr::group_by(Stock, type) %>%
                      dplyr::summarize(count=length(unique(value))) %>%
                      dplyr::group_by(Stock) %>%
                      dplyr::summarize(types=length(unique(type))) %>%
                      dplyr::filter(types==5) %>%
                      dplyr::select(Stock)
  # The warning here is the tot field, not worried about that.
  ASR.stocks <- ASR.long %>%
                          dplyr::filter(Stock %in% Stocks$Stock) %>%
                          dplyr::filter(!is.na(value)) %>%
                          dplyr::arrange(Stock, Year, type, as.numeric(age))
  
  # Stocks <- ASR_stocks %>%
  #                       dplyr::group_by(Stock, Species, type) %>%
  #                       dplyr::summarize(ages=length(unique(age)),
  #                                 years=length(unique(Year))) %>%
  #                       dplyr::arrange(-years, -ages) %>%
  #                       dplyr::filter(!Species=="morhua") %>%
  #                       dplyr::distinct(Stock) %>%
  #                       dplyr::pull(Stock)
  # 
  # print(Stocks)
  Stocks <- Stocks$Stock
  
  # Going to remove Sebastes norvegicus because we can't get fecundity as the only data we have is total numbers.
  #Stocks <- Stocks[Stocks != "ICES-AFWG_DEEP1-2_Sebastes_norvegicus"]
  # Esmarki doesn't work great as it is assessed using a quarterly model.
  #Stocks <- Stocks[Stocks != "ICES-WGNSSK_NS 4-3aN_Trisopterus_esmarkii"]
  #i = 'ICES-AFWG_NEA1-2_Melanogrammus_aeglefinus'
years.tmp <- NULL
pnm.tmp <- NULL
waa.tmp <- NULL
ages.tmp <- NULL
rem.tmp <- NULL
mx.tmp <- NULL
NE.tmp <- NULL
vpa.tmp <- NULL
am.tmp <- NULL
mr.tmp <- NULL
all.abund.tmp <- NULL
full.abund.tmp <- NULL
full.rem.tmp <- NULL
FM.tmp <- NULL
abund.tmp <- NULL

for(i in Stocks)
{
  print(i)
  ASR.sub <- ASR.long %>%
    dplyr::select(Year, Stock, type, age, value) %>%
    dplyr::filter(Stock==i)

  
  age.mat <- ASR.sub %>% dplyr::filter(type=="AM") %>% dplyr::rename(AM=value) %>% dplyr::select(-type)
  # All of the age at maturities are NA in 0 and 1 year olds, so make them 0's, do that carefully just
  # in case that changes later...
  am.0s <- length(which((!is.na(age.mat$AM[age.mat$age ==0]))))
  am.1s <- length(which((!is.na(age.mat$AM[age.mat$age ==1]))))
  am.2s <- length(which((!is.na(age.mat$AM[age.mat$age ==2]))))
  am.3s <- length(which((!is.na(age.mat$AM[age.mat$age ==3]))))
  if(am.0s == 0) age.mat$AM[age.mat$age ==0] <- 0
  if(am.1s == 0) age.mat$AM[age.mat$age ==1] <- 0
  if(am.2s == 0) age.mat$AM[age.mat$age ==2] <- 0
  if(am.3s == 0) age.mat$AM[age.mat$age ==3] <- 0
  
   # Now we can get the abundance from the VPA models
  abund <- ASR.sub %>% dplyr::filter(type=="Num") %>% dplyr::rename(Num=value) %>% dplyr::select(-type)
  abund.no.tot <- abund[abund$age != 'tot',]
  removals <-   ASR.sub %>% dplyr::filter(type=="FM") %>% dplyr::rename(FM=value) %>% dplyr::select(-type)
  weight.age <- ASR.sub |> collapse::fsubset(type=="WA") |> dplyr::rename(WA=value) |> collapse::fselect(-type)
  removals$er <- 1-exp(-removals$FM)
  removals$catch <- abund.no.tot$Num*removals$er
  nat.mort <-  ASR.sub |> collapse::fsubset(type=="NM") |> dplyr::rename(NM=value) |> collapse::fselect(-type)
  
  #if(am.1s | am.0s | am.2s > 0) print("Stop, you need to check the age at maturity for either age 0 or 1 as there is data in there.")
  # Now do something similar for the Natural moralities, but in this case
  # we take the NM of the youngest age class we have information for and back that number up to year 0.
  # DK NOTE:  THIS IS A BIG CHANGE IN DIRECTION FOR THE LOTKA CALCS SO NEED TO NOTE THIS ONE!!
  
  data <- age.mat |>
    full_join(nat.mort) |>
    full_join(abund) |>
    full_join(weight.age) |>
    full_join(removals)
  
  # Aahh??  Is this NA dropping useful data, needed to quick fix the age at maturity data...
  data$available <- apply(is.na(data[, c("AM", "Num", "WA")]), 1, function(x) all(!x==T))
  data <- data[data$available==T,]
  
  # Tidy up the data for input...
  data$prop.nat.mort <- 1-exp(-data$NM)
  #prop.nat.mort[,-which(names(prop.nat.mort) %in% c("Year", "Stock", "type"))] <- 1-exp(-prop.nat.mort[,-which(names(prop.nat.mort) %in% c("Year", "Stock", "type"))])
  
  
  rem <- data |> collapse::fgroup_by(Year) |> collapse::fsummarize(rem=sum(catch,na.rm=T)) #|> dplyr::pull(rem)
  #rem <- data |> dplyr::group_by(Year,.drop=F) |> collapse::fsummarize(rem=sum(Catch,na.rm=T)) #|> dplyr::pull(rem)
  
  # Never run, so not doing anything.
  # missing_rem<- NULL
  # if(any(is.na(rem$rem))) 
  # {
  #   missing_rem <- unique(data$Year)[which(is.na(rem$rem))]
  #   rem$rem[is.na(rem$rem)] <- median(rem$rem, na.rm=T)
  # }
  
  #rowSums(removals[,-which(names(removals) %in% c("Year", "Stock", "type"))], na.rm=T)
  years <- data |> dplyr::pull(Year) |> unique() |> sort()
  
  N.end <- sum(data[data$Year==max(years),]$Num)
  vpa.abund <- data |> dplyr::group_by(Year) |> dplyr::summarize(vpa=sum(Num,na.rm=T)) |> pull(vpa)
  
  # The real mx matrix, recruits produced per individual in each age class... Not perfect as I need to offset recruits/ssb, but close enough for the moment..
  #minage
  minage <- min(as.numeric(data$age))
  #maxage
  maxage <- max(as.numeric(data$age))
  #recruits
  annual <- data.frame(Year=data$Year[data$age==minage], recruits=data$Num[data$age==minage])
  #ssn
  data$ssn <- data$Num * data$AM
  #ssb
  data$ssb <- data$ssn * data$WA #I need to figure out how to line this up with the number of recruits for the right years in the data object.... got it right for the overall r.p.ssb below.
  #tot.ssb
  annual <- data |> group_by(Year) |> summarize(tot.ssb = sum(ssb)) |> left_join(annual)
  
  #tst <- data |> group_by(Year) |> summarize(tot.ssn = sum(ssn)) |> left_join(annual)
  
  #r.p.ssb This needs to be offset by the age of recruits. So here's my lining this up logic.  For age 0, we get
  # the recruit and SSB estimates from the same survey, thus these age 0 recruits will have been produced by the SSB in the previous year...
  # So the R in 1970 were produced by SSB in 1969 when min age =0
  # If min.age = 1 then R in 1970 is produced by SSB in 1968... and so on.
  # FIX 1 "Recruit year" = SSB year + min.age + 1 (or min.age +2 as SSB.year is the first column)
  # Now for RPS there's another wrinkle because of this, because our model is N(1970) = lambda * N(1969) and lambda is made up of lx and mx
  # so if the recruits show up in 1970, the fecundity term that gets us there (irrespective of the age offset) has to be for 1969, becuase
  # that gets us to the babies that show up in 1970
  # FIX 2 So "RPS year, recs.per.age year, and mx year" = SSB year + min.age
  # I believe this to be correct, but we all need to discuss this to make sure the logic holds.
  annual$r.p.ssb <- c(rep(NA,minage),annual$recruits[(minage+2):nrow(annual)]/annual$tot.ssb[1:(nrow(annual)-minage-1)],NA)
  # recs.per.age
  data <- left_join(data, annual)
  # So I need to make an offset SSB in the 'data' object to line up with the correct r.p.ssb field.  This is probably gonna suck...
  tmp <- NULL
  for(j in 1:length(years))
  {
    tst <- data|> collapse::fsubset(Year == years[j]) |> collapse::fselect(ssb,ssn,Year,age)
    tst$Year <- tst$Year + minage # We want this to line up the the RPS year, so just adding minage here.
    names(tst) <- c("ssb.offset","ssn.offset","Year","age")
    tmp[[j]] <- tst
  }
  # Unpack the list
  ssb.off <- do.call('rbind',tmp)
  # And merge it with the data object
  data <- left_join(data,ssb.off,by=c("Year","age"))
  data$recs.per.age <- data$ssb.offset*data$r.p.ssb
  # mx
  # Some of the stocks are just the males and females, so these are already half the population.  Probably can just not include the male stocks in the end.
  # No longer dividing by 2 because we should really divide recruits and ssn by 2 to be female only, but that's pointless so we roll along with this.
  data$mx <- data$recs.per.age/data$ssn.offset # if(grepl("males",i)) 
  #if(!grepl("males",i))data$mx <- data$recs.per.age/data$ssn.offset/2 # Moms only! 
  data$mx[is.nan(data$mx)] <- 0 # if we don't have any spawners in an age class in a year their fecundity is 0
  # Easier to remove the years where we don't have recs.per.age...
  #data <- data |> collapse::fsubset(Year %in% years[(minage+1):length(years)])
  
  age.mat <- data |> collapse::fselect("Year", "age", "AM") |> 
                      pivot_wider(names_from=age, values_from = AM) |> collapse::fselect(-Year)
  
  prop.nat.mort <- data |> collapse::fselect("Year", "age", "prop.nat.mort") |> 
                      pivot_wider(names_from=age, values_from = prop.nat.mort) |> 
                      collapse::fselect(-Year) |> as.data.frame()
  
  weight.age <- data |> collapse::fselect("Year", "age", "WA") |> pivot_wider(names_from=age, values_from = WA) |> collapse::fselect(-Year)
  mx <- data |> collapse::fselect("Year", "age", "mx") |> 
                pivot_wider(names_from=age, values_from = mx) |> 
                collapse::fselect(-Year) |> as.data.frame()
  # Get the abundance matrix
  abund.age <- data |> collapse::fselect("Year", "age", "Num") |> 
    pivot_wider(names_from=age, values_from = Num) |> 
    collapse::fselect(-Year) |> as.data.frame()

  full.rems <- data |> collapse::fselect("Year", "age", "catch") |> pivot_wider(names_from=age, values_from = catch) |> collapse::fselect(-Year)
  full.Fs <-  data |> collapse::fselect("Year", "age", "FM") |> pivot_wider(names_from=age, values_from = FM) |> collapse::fselect(-Year)
  
  
  if(minage >0) 
  {
    rm <- -(1:minage)
    years <- years[rm]
    prop.nat.mort <- prop.nat.mort[rm,]
    weight.age <- weight.age[rm,]
    rem <- rem[rm,]
    mx <- mx[rm,]
    vpa.abund <- vpa.abund[rm]
    age.mat <- age.mat[rm,]
    abund.age <- abund.age[rm,]
    full.Fs <- full.Fs[rm,]
    full.rems <- full.rems[rm,]
  }
  #   mx.fill <- as.data.frame(matrix(rep(colMeans(mx,na.rm=T),minage),nrow=minage,byrow=T),colnames = names(mx))
  #   names(mx) <- names(mx.fill)
  #  # mx <- rbind(mx[(minage+1):nrow(mx),],mx.fill) # or is it
  #   mx <- rbind(mx.fill,mx[(minage+1):nrow(mx),]) # I think this is right!
  # }
  
  
  years.tmp[[i]] <- years
  pnm.tmp[[i]] <- prop.nat.mort
  waa.tmp[[i]] <- weight.age
  ages.tmp[[i]] <- minage:maxage
  rem.tmp[[i]] <- rem
  mx.tmp[[i]] <- mx
  NE.tmp[[i]] <- N.end
  vpa.tmp[[i]] <- vpa.abund
  am.tmp[[i]] <- age.mat
  #mr.tmp[[i]] <- missing_rem
  abund.tmp[[i]] <- abund.age
  FM.tmp[[i]] <- full.Fs
  full.rem.tmp[[i]] <- full.rems
} #end input data loop

save(years.tmp,pnm.tmp,waa.tmp,ages.tmp,rem.tmp,mx.tmp,NE.tmp,vpa.tmp,mr.tmp,am.tmp,Stocks,ASR.stocks,ASR.long,all.abund.tmp,full.abund.tmp,FM.tmp,
     full.rem.tmp,file =  paste0(loc,"/Results/model_inputs.Rdata"))

   
  
   
  
  

  
  ########################## Now run the tuning sims
  #
  
  load(file = paste0(loc,"/Results/model_inputs_all_ages.Rdata"))
  
  for.tune.all.res <- NULL
  for.tune.summary <- NULL
  for(i in Stocks)
  {
    years <- years.tmp[[i]]
    prop.nat.mort <- pnm.tmp[[i]] 
    #rem <- rem.tmp[[i]] 
    mx <- mx.tmp[[i]] 
    N.end <- NE.tmp[[i]] 
    ages <- ages.tmp[[i]]
    vpa.abund <- vpa.tmp[[i]] 
    #all.abund <- vpa.tmp[[i]]
    fs <-FM.tmp[[i]]
    full.rem <- full.rem.tmp[[i]]
    full.abund <- abund.tmp[[i]]
    
    #all.abund <- full.abund.tmp[[i]]
    #N.end <- vpa.abund[length(vpa.abund)]
    #Pick your step size
    ss.size <- 0.001
    # sensistive.stocks <- c("ICES-WGNSSK_NA_Trisopterus_esmarkii",
    #                        "ICES-HAWG_NS_Ammodytes_marinus",
    #                        "ICES-WGNSSK_NS NW_Gadus_morhua",
    #                        "ICES-WGNSSK_NS South_Gadus_morhua",
    #                        "ICES-WGNSSK_NS Viking_Gadus_morhua",
    #                        "ICES-WGNSSK_NS_Melanogrammus_aeglefinus",
    #                        "ICES-WGNSSK_NS_Merlangius_merlangus",
    #                        "ICES-HAWG_NS 1r_Ammodytes_marinus")
    # 
    # if(i %in% sensistive.stocks) ss.size = ss.size/10
  
    N.start <- vpa.abund[1]
     print(i)
     # The fecundities are 
     fast.for.tune <- fast.tunes(years,
                             tuner=what.2.tune,
                             step.size = ss.size, 
                             abund.ts = vpa.abund,
                             ages=0:max(ages),
                             nm = -(log(1-prop.nat.mort)),
                             fecund = mx,
                             fm = fs,
                             abund.age = full.abund, 
                             catch.age = full.rem, 
                             N.init = N.start,
                             direction= 'forwards'
                             )
     
     for.tune.all.res[[i]] <- fast.for.tune
     for.tune.tmp <- data.frame(fast.for.tune$res,Stock = i)
     for.tune.summary[[i]] <- for.tune.tmp
     
     N.end <- full.abund[length(full.abund)]

  }
  
  
  
  tst <- NULL
  tst1 <- NULL
  for(i in names(for.tune.all.res))
  {
    if(i != "ICES HAWG_NS-IV 3a,7d_Clupea_harengus") tst[[i]] <- for.tune.all.res[[i]]
    if(i != "ICES HAWG_NS-IV 3a,7d_Clupea_harengus") tst1[[i]] <- for.tune.summary[[i]]
  }
  
  for.tune.all.res <- tst
  for.tune.summary <- tst1
  #DFO_2J3KL_Gadus_morhua -27.7984654 2.815738e+08
  #ICES-WGCSE_IS6a-7b-7j_Dicentrarchus _labrax -21.6050608
  # ICES-WGCSE_IS6a-7b-7j_Dicentrarchus _labrax Lotka error
  # ICES-HAWG_CS 6a- 7b-7c_Clupea_harengus # lotka error
  
  
  #saveRDS(for.tune.summary,"D:/Github/ICM/Results/forwards_fully_tuned_summary_z.Rds")
  #saveRDS(for.tune.all.res,"D:/Github/ICM/Results/forwards_fully_tuned_all_res_z.Rds")
  # saveRDS(for.tune.summary,"D:/Github/ICM/Results/forwards_fully_tuned_summary_nm.Rds")
  # saveRDS(for.tune.all.res,"D:/Github/ICM/Results/forwards_fully_tuned_all_res_nm.Rds")
  #saveRDS(for.tune.summary,"D:/Github/ICM/Results/forwards_fully_tuned_summary_fm.Rds")
  #saveRDS(for.tune.all.res,"D:/Github/ICM/Results/forwards_fully_tuned_all_res_fm.Rds")
  saveRDS(for.tune.summary,paste0(loc,"/Results/forwards_fully_tuned_summary.Rds"))
  saveRDS(for.tune.all.res,paste0(loc,"/Results/forwards_fully_tuned_all_res.Rds"))
  # 


  # We are now looking at everything, lifetime reproductive success
  
  #for.tune.summary <- readRDS(paste0(loc,"/Results/forwards_fully_tuned_summary_",what.2.tune,".Rds"))
  #for.tune.all <- readRDS(paste0(loc,"/Results/forwards_fully_tuned_all_res_",what.2.tune,".Rds"))
  
# Get the no fish and fish info we need for later
for.tune.all <- for.tune.all.res
  
  res.fin <- NULL
  res.year.fin <- NULL
  res.lambda.fin <- NULL
  for(i in Stocks)
  {
    tst <- for.tune.all[[i]]
    meta.tmp <- i
    # Fishing mortality
    fm.org.tmp <- tst$fm.org
    fm.tmp <- tst$fm.opt
    # fecundity
    fec.tmp <- tst$fecund.opt
    fec.org.tmp <- tst$fecund.org
    # Natural mortality
    m.org.tmp <- exp(-tst$nm.org) # This is survivorship
    m.tmp <- exp(-tst$nm.opt) # This is survivorship
    # Total mortality
    z.tmp <- exp(-tst$z.opt) # This is survivorship
    z.org.tmp <- exp(-tst$nm.org -fm.org.tmp) # This is survivorship
    
    ages <- 0:(ncol(fec.tmp)-1)
    n.ages <- length(ages)
    years <- tst$res$year
    n.years <- length(years)
    # Now I need to grab the fecundities and nat morts by cohort
    if(n.years > n.ages)
    {
      res.cohort <- NULL
      res.year <- NULL
      res.lambda <- NULL
      for(y in 1:(n.years-n.ages+1))
      {
        count = 0
        fs <- NA
        # The lx is no fishing survivorship, zx I'm calling the survivorship including fishing
        lx <- lx.org <- lx.year <- zx <- zx.org <- zx.year <- 1
        fs.org <- NA
        
        for(a in 1:n.ages)
        {
          fs[a] <- fec.tmp[y+count,a]  
          fs.org[a] <- fec.org.tmp[y+count,a]  
          if(a > 1) 
          {
            # No fishing
            lx[a] <- lx[a-1] * m.tmp[y+count,a]  
            lx.org[a] <- lx.org[a-1]* m.org.tmp[y+count,a]  
            lx.year[a] <- lx.year[a-1] * m.tmp[y,a]
            # Including fishing
            zx[a] <- zx[a-1] * z.tmp[y+count,a]  
            zx.org[a] <- zx.org[a-1]* z.org.tmp[y+count,a]  
            zx.year[a] <- zx.year[a-1] * z.tmp[y,a]
          }  # end the a > 1 if
          count <- count + 1
        }# end the ages loop
        res.cohort[[y]] <- data.frame(mx.opt = fs,lx.opt = lx,mx.org = fs.org,lx.opt = lx,
                                      mx.lx.opt = fs*lx,mx.lx.org = fs.org*lx.org,
                                      x.mx.lx.opt = fs*lx*ages,x.mx.lx.org = fs.org*lx.org*ages,
                                      zx.opt = zx,zx.opt = zx,
                                      mx.zx.opt = fs*zx,mx.zx.org = fs.org*zx.org,
                                      x.mx.zx.opt = fs*zx*ages,x.mx.zx.org = fs.org*zx.org*ages,
                                      cohort=years[y],age = ages,stock = meta.tmp)
        
        res.year[[y]] <- data.frame(mx.opt = unlist(fec.tmp[y,]),lx.opt = lx.year,
                                    mx.lx.opt = unlist(fec.tmp[y,])*lx.year,
                                    x.mx.lx.opt = unlist(fec.tmp[y,])*lx.year*ages,
                                    zx.opt = zx.year,
                                    mx.zx.opt = unlist(fec.tmp[y,])*lx.year,
                                    x.mx.zx.opt = unlist(fec.tmp[y,])*lx.year*ages,
                                    age = ages,year=years[y],stock = meta.tmp)
        
      } # end the years loop
      
      # Do a years loop to get lambda with no fishing for each year, rather than by cohort...
      for(yy in 1:(n.years-1))
      {
        #if(yy == n.years) browser()
        nm <- as.numeric(tst$nm.opt[yy,])
        z <- as.numeric(tst$z.opt[yy,])
        fec <- as.numeric(tst$fecund.opt[yy,])
        nm.init <- as.numeric(tst$nm.org[yy,])
        fec.init <- as.numeric(tst$fecund.org[yy,])
        z.init <- as.numeric(tst$z.org[yy,])
        lam.no.fish <- exp(simple.lotka.r(mort = nm,fecund = fec,ages=ages)$res$r)
        lam.fish <- exp(simple.lotka.r(mort = z,fecund = fec,ages=ages)$res$r)
        lam.org.no.fish <- exp(simple.lotka.r(mort = nm.init,fecund = fec.init,ages=ages)$res$r)
        lam.org.fish <- exp(simple.lotka.r(mort =  z.init,fecund = fec.init,ages=ages)$res$r)
        
        # Now summarize these
        res.lambda[[yy]] <-  data.frame(lam.no.fish = lam.no.fish,lam.fish = lam.fish,
                                        lam.org.no.fish = lam.org.no.fish,
                                        lam.org.fish = lam.org.fish,
                                        year=years[yy],stock = meta.tmp)
      } # End the yy loop to get the lambdas...
      # unpack it...
      res.fin[[i]] <- do.call('rbind',res.cohort)
      res.year.fin[[i]] <- do.call('rbind',res.year)
      res.lambda.fin[[i]] <-  do.call('rbind',res.lambda)
    }  # end if to make sure we have enough data to get the cohort estimates.
    
  }# end the stock loop
  
  res.cohort.final <- do.call('rbind',res.fin)
  res.year.final <- do.call('rbind',res.year.fin)
  res.lambda.final <- do.call('rbind',res.lambda.fin)

  saveRDS(res.lambda.final,paste0(loc,"/Results/lambdas.Rds"))
  saveRDS(res.cohort.final,paste0(loc,"/Results/cohort_parameters.Rds"))
  saveRDS(res.year.final,paste0(loc,"/Results/annual_parameters.Rds"))

