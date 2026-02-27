# Here we develop a multi-species model for the North Sea.

#################  Section 1 Loading #################  Section 1 Loading #################  Section 1 Loading  ###############################################
library(tidyverse)
library(GGally)
library(cowplot)
library(ggthemes)
# Set the base plot theme
theme_set(theme_few(base_size = 14))
options(scipen = 999)
# Download the function to go from inla to sf
funs <- c("https://raw.githubusercontent.com/dave-keith/ICM/main/Scripts/functions/simple_Lotka_r.r",
          "https://raw.githubusercontent.com/dave-keith/ICM/main/Scripts/functions/simple_forward_sim.r",
          "https://raw.githubusercontent.com/dave-keith/ICM/main/Scripts/functions/forward_project.r"
          
)
# Now run through a quick loop to load each one, just be sure that your working directory is read/write!
for(fun in funs) 
{
  download.file(fun,destfile = basename(fun))
  source(paste0(getwd(),"/",basename(fun)))
  file.remove(paste0(getwd(),"/",basename(fun)))
}



#load(file = "D:/Github/ICM/Results/model_inputs.Rdata")
dat.loc <- 'D:/GitHub/ICM'
repo.loc <- "D:/GitHub/Multispecies_model/"
#loc <- "C:/Users/Owner/Documents/Github/ICM"
load(file = paste0(dat.loc,"/Results/all_cleaned_forward_tune_summaries_fec_nm.Rdata"))
load(file = paste0(dat.loc,"/Results/model_inputs.Rdata"))

########################### End Section 1 Loading  ###################################### End Section 1 Loading ###############################################

########################## Section 2 Parameters ########################## Section 2 Parameters ########################## Section 2 Parameters

n.yrs.proj <- 25 # How many years into the future we are going to project the stocks
n.sims <- 10 # The numbers of simulations to run, keeping low for testing...

# Get the right stocks
Stocks <- names(for.tune.all)
Stocks <- Stocks[grep("NS",Stocks)]
Stocks <- Stocks[Stocks != "ICES-WGHANSA_SP8abd_Sardina _pilchardus"]


# So here we are working to get the 'ecosystem' carrying capacity by looking at the total biomass for the NS stocks we have
# data for over the period of time we have data for all the stocks.
# So here we pull out the data we need to look at total abundance and total biomass in the system by year...
years.ns <- NULL
vpa.ns <- NULL
bm.ns <- NULL
num.ns <- NULL
waa.ns <- NULL
pnm.ns <- NULL
rem.ns <- NULL
mx.ns <- NULL
am.ns <- NULL
ages.ns <- NULL
for(i in  Stocks)
{
  years.ns[[i]] <- years.tmp[[i]]
  #vpa.ns[[i]] <- vpa.tmp[[i]]
  ages.ns[[i]] <- ages.tmp[[i]]
  num.ns[[i]] <- ASR_long |> collapse::fsubset(Stock == i & type == "Num")
  num.ns[[i]] <- num.ns[[i]] |> collapse::fsubset(age != "tot")
  waa.ns[[i]] <- ASR_long|> collapse::fsubset(Stock == i & type == "WA")
  if(i == "ICES-HAWG_NS_Ammodytes_dubius") waa.ns[[i]]$value <- waa.ns[[i]]$value/1000
  bm.ns[[i]] <- data.frame(Year = num.ns[[i]]$Year,Stock = num.ns[[i]]$Stock,age = num.ns[[i]]$age,
                           bm = num.ns[[i]]$value*waa.ns[[i]]$value,
                           num = num.ns[[i]]$value)
  pnm.ns[[i]] <- 1-exp(-for.tune.all[[i]]$nm.opt)
  mx.ns[[i]] <- for.tune.all[[i]]$fecund.opt
  vpa.ns[[i]] <- for.tune.all[[i]]$res$est.abund
  rem.ns[[i]] <- rem.tmp[[i]]
  rem.ns[[i]]$Stock <- i
  am.ns[[i]] <- am.tmp[[i]]
}
# Combine the biomass and abundance data into a dataframe
bm.tst <- do.call("rbind",bm.ns)
# Look at the biomass and abundance in the ecosystem
bm.tot <- bm.tst |> collapse::fgroup_by(Stock,Year) |> 
                    collapse::fsummarize(bm = sum(bm,na.rm=T),
                                         num = sum(num,na.rm=T))
# The 'ecosystem' biomass and numbers
eco.bm <- bm.tot |> collapse::fgroup_by(Year) |> 
                    collapse::fsummarize(num = sum(num),bm = sum(bm))


# Now we combine the ecosystem results with the stock biomass's
bm.final <- left_join(bm.tot,eco.bm,by=c("Year"))
names(bm.final) <- c("Stock","Year","bm.stock","num.stock","num.total","bm.total")
# Get the proportion of the total biomass each stock accounts for
bm.final <- bm.final |> collapse::fmutate(bm.prop = bm.stock/bm.total,
                                       num.prop = num.stock/num.total)
bm.final <- bm.final[bm.final$bm.stock > 0,]
bm.final <- as.data.frame(bm.final)
# This gets the average weight of individuals in each stock, we'll need this later to get an approximate exploitation rate
bm.final$avg.weight <- bm.final$bm.stock/bm.final$num.stock

# Now we subset to the years we have data for all the stocks
what.year <- bm.final |> collapse::fgroup_by(Stock) |> collapse::fsummarize(min = min(Year),
                                                                      max = max(Year))
# The years we have data for all stocks
first.year <- max(what.year$min)
last.year <- min(what.year$max)

# Now we subset the data to these years
bm.best <- bm.final |> collapse::fsubset(Year %in% first.year:last.year & Stock == Stocks[1], bm.total, num.total,Year) 

# Fix: This is not perfect way to get the past exploitation rates as the removals we have here are in numbers
# we'll need to think about whether using the avg.weight is appropriate it does get us 'close', but I think we can do better given the data we have.
rem.tst <- do.call("rbind",rem.ns)
fm.dat <- left_join(bm.final,rem.tst,by=c("Stock",'Year'))
# This is where we go from numbers to a biomass and get an exploitation rate in biomass.
fm.dat$exploit <- (fm.dat$rem*fm.dat$avg.weight)/fm.dat$bm.stock

# # Autocorrelation in abundance and biomass time series for the 'ecosystem'
K.cor <- pacf(log(bm.best$num.total))
K.cor.bm <- pacf(log(bm.best$bm.total))
# We see there is a decent correlation in the biomass time series, which is kinda nice, the ecosystem biomass isn't just jumping around randomly....
# visualizing the biomass time series, reasonable decline in the 'biomass in the early 2000s.
ggplot(bm.best,aes(x=Year,y=bm.total)) + geom_point() + geom_line()
# Numbers time series, one reason numbers doesn't work is the total domination by Sand Lance.
ggplot(bm.best,aes(x=Year,y=num.total)) + geom_point() + geom_line()


# So now we can 'statistically represent the abundance and biomass carrying capacity, note we do this on the log scale
N.target.sum <- bm.best |> collapse::fsummarise(sd = sd(log(num.total)),
                                            mn = mean(log(num.total)),
                                            med = median(log(num.total)))
# Working on implementing using the biomass as K rather than nubmers...
bm.target.sum <- bm.best |>   collapse::fsummarise(sd = sd(log(bm.total)),
                                               mn = mean(log(bm.total)),
                                               med = median(log(bm.total)))


# So here is where the fun comes in, given we have the autocorrelation, mean and standard deviation of the carrying capacity time series
# We can make a time series with this same characteristics, which we can use to simulation the future carrying capacity of the system.
# We can represent the carrying capacity in numbers or in biomass.  
# FIX: The ICM runs in numbers, carrying capacity is probably better thought of in biomass, think of justifiable ways to get ICM in biomass
# or logical ways to convert result of ICM to biomass (we have the weight at age, but for now doing it this way)
K.devs <- NULL
cors <- NULL
K.sims <- NULL
for(i in 1:n.sims) 
{
   K.devs[[i]] <- as.data.frame(arima.sim(list(order = c(1,0,0), ar = K.cor$acf[1]), 
                                                                 n = n.yrs.proj,
                                                                 sd = N.target.sum$sd))
  #K.devs[[i]] <- data.frame(x = rep(0,n.yrs.proj)) # Lets see what happens when it's fixed at a mean level.
  K.sims[[i]] <- data.frame(year = 2015:2039,sim = i,num = exp(N.target.sum$med + as.numeric(K.devs[[i]]$x)))
}

K.sims.4.plt <- data.frame(do.call('rbind',K.sims))
# Now does this work, lets visualize this...
# FIX: this is not bad, only issue I have is the starting point probably should be tied to the value in 2014, so will need to think about that
# The simulation was also skewing a bit high, switching from using the mean to the median when getting K.sims took care of that.
ggplot(K.sims.4.plt,aes(x=year,y=num)) + 
                        geom_point(aes(group=sim,color=sim,fill=sim)) + geom_line(aes(group=sim,color=sim,fill=sim)) + 
                        geom_line(data=bm.best,aes(x=Year,y=num.total),color='black',linewidth=2)

# FIX: This is a simplification that we'll need to consider, we might want to make this
# be where we figure out how to partition the biomass by trophic level.  But the moment it's ok.
# For the moment, lets just cut the world up so each stock gets a fixed % of K based on their biomass
prop.stock <- bm.final |> collapse::fgroup_by(Stock) |> collapse::fsummarize(prop = median(num.prop))
prop.stock.bm <- bm.final |> collapse::fgroup_by(Stock) |> collapse::fsummarize(prop = median(bm.prop))

# Now we use the total carrying capacity and the proportion of the abundance of each stock gets to get a carrying capacity of each stock
# over the time series.
K.stock <- NULL
fm.stock <- NULL
tmp <- NULL
for(i in 1:n.sims)
{
  tmp <- NULL
  for(s in Stocks)
  {
    # Pull out the proportion 
    prop <- prop.stock |> collapse::fsubset(Stock == s)
    tmp[[s]] <- prop$prop * as.numeric(K.sims[[i]]$num)
    # FIX: Really this should be in biomass too, but for now sticking with using abundance
    if(i == 1) fm.stock[[s]] <- fm.dat |> collapse::fsubset(Stock ==s) |> collapse::fsummarize(mn = median(exploit,na.rm=T),
                                                                                         sd = sd(log(exploit[exploit > 0]),na.rm=T))
  }
  K.stock[[i]] <- tmp
}

# Plots you are welcome to enjoy at your leisure
# ggplot(bm.best) + geom_line(aes(x=Year,y=bm.total),linewidth=2) + ylim(c(0,max(bm.final$bm.total))) + geom_hline(yintercept = mean(bm.best$bm.total))
# ggplot(bm.best) + geom_line(aes(x=Year,y=num.total)) + ylim(c(0,max(bm.final$num.total))) 
# ggplot(bm.best) + geom_line(aes(x=Year,y=bm.prop,group = Stock,color=Stock),linewidth=2) #+ ylim(c(0,0.5))
# ggplot(bm.best) + geom_bar(aes(x=Year,y=bm.prop,fill=Stock),position="stack",stat = 'identity') #+ ylim(c(0,0.5))
# prop.acfs <- NULL
# for(i in Stocks) prop.acfs[[i]] <- pacf(bm.best$bm.prop[bm.best$Stock == i])




############### Section 4 Multi-species model of North Sea ############### Section 4 Multi-species model of North Sea ############### Section 4 Multi-species model of North Sea

# Now we have our carrying capacity for each stock and we can get to business and running a model.
tmp.mx <- NULL
tmp.nm <- NULL
tmp.mat <- NULL
tmp.age <- NULL
tmp.waa <- NULL
mx.dev <- NULL
nm.dev <- NULL
ts.unpack <- NULL
r.unpack <- NULL
Ks <- NULL

# Get the year range, going from the 'last' year to n.yrs.proj in the future, note this will go 1 year less than your intuition because
# we want n.yrs of data, i.e., 20 years is 2000 to 2019, not 2020... )
years <- (last.year+1):(last.year+n.yrs.proj)

# So everything will need to get wrapped up in a simulation loop
for(j in 1:n.sims)
{
  st.time <- Sys.time()
  # This is getting the K for the stock for a particular simulation
  Kss <- K.stock[[j]]
  
  res.ts <- NULL
  res.r <- NULL
  for(s in Stocks)
  {
    # Obtain a sample of the correct size from the historic natural mortality and fecundities, Note
    # That in this set up this picks the fecundity and natural mortality from the 'same' year, which is 'fine', but lots of options here.
    # FIX: We could do this sample as a proper statistical sample rather than picking exactly what has been observed in the past.
    samp <- sample(nrow(mx.ns[[s]]),n.yrs.proj,replace=T)
    # Historical natural mortality estimates (using our retrospective simulations)
    prop.nat.mort <- pnm.ns[[s]] 
    prop.nat.mort <- prop.nat.mort[samp,] # Get the sample years.
    
    # maturity ogive
    #age.mat <- am.ns[[s]]
    # Historical fecundity estimates (using our retrospective simulations)
    mx <- mx.ns[[s]] 
    mx <- mx[samp,] # Get the sample years.
    # Ages, everything is calculated back to 'age 0', needed so lotka optimization 'knows' the age it's dealing with
    ages <- 0:(ncol(prop.nat.mort)-1)
    # Carrying capacity for the stocks, currently in numbers, need to think about it in biomass terms IMHO, but not there yet.
    Ks <- Kss[[s]]
    # The number of individuals (using our retrospective simulations)
    vpa.ns  <- bm.final$num.stock[bm.final$Stock == s]
    # Now get the final year abundance
    N.start <- vpa.ns[length(vpa.ns)]
    
    # Run the simulations
    tst <- simp.for.sim(years= years,
                        nm = -(log(1-prop.nat.mort)),
                        ages = ages,
                        rems =  list(fm.stock[[s]]$mn,fm.stock[[s]]$sd), #fm,
                        fecund = mx,
                        N.start = N.start,
                        pop.model = 'bounded_exp', 
                        sim= "project",
                        n.sims = 1,
                        K = Ks,
                        repo='preload')
    
# Tidy up the output, Abundance and population growth rates
res.ts[[s]] <- data.frame(tst$Pop[,-2],stock = s,sim= j)
res.r[[s]] <- data.frame(tst$r[,-3],stock=s,sim=j)

  } # end stock loop
  # Unpack the results
  ts.unpack[[j]] <- do.call('rbind',res.ts)
  r.unpack[[j]] <- do.call('rbind',res.r)
  
  # Pop a note when done each simulation
  timer <- Sys.time() - st.time
  print(paste("Simulation ", j))
  print(signif(timer,digits=2))
  
} # end n.sims

# Unpack all the results.
ts.final <- do.call("rbind",ts.unpack)
ts.final$fm <- ts.final$removals/ts.final$abund
r.final <- do.call("rbind",r.unpack)

quants <- ts.final |>  collapse::fgroup_by(years,stock) |> collapse::fsummarize(L.50 = quantile(abund,probs=c(0.25),na.rm=T),
                                                                          med = median(abund,na.rm=T),
                                                                          U.50 = quantile(abund,probs=c(0.75),na.rm=T),
                                                                          fml.50 = quantile(fm,probs=c(0.25),na.rm=T),
                                                                          fm = median(fm,na.rm=T),
                                                                          fmu.50 = quantile(fm,probs=c(0.75),na.rm=T))
# If happy save the 3 objects
saveRDS(object = ts.final,file = paste0(repo.loc,"/Results/NS_projections_",n.sims,"_sims_",min(years),"_to_",max(years),
                                         "_time_series_projections.Rds"))

saveRDS(object = quants,file = paste0(repo.loc,"/Results/NS_projections_",n.sims,"_sims_",min(years),"_to_",max(years),
                                      "_time_series_quantiles.Rds"))

saveRDS(object = r.final,file = paste0(repo.loc,"/Results/NS_projections_",n.sims,"_sims_",min(years),"_to_",max(years),
                                       "_r_projections.Rds"))


# Two simple plots. 
p.sims <- ggplot(ts.final ) + geom_line(aes(x=years,y=abund,group = sim,color=sim),alpha=0.8) +
  facet_wrap(~stock,scales = 'free_y') + ylim(c(0,NA)) + scale_x_continuous(breaks = seq(2010,2200,by=10)) 
#save_plot(paste0("D:/Github/ICM/Figures/NS_sims/NS_all_realizations_climate_starts_at_year_",c.effect,
#                 "_mx_decade_effect_",climate.mx.effect, "_nm_decade_effect_",climate.nm.effect,".png"),p.sims,base_height = 12,base_width = 20)

p.sims.quants <- ggplot(quants) + geom_line(aes(x=years,y=med)) + facet_wrap(~stock,scales = 'free_y') + ylim(c(0,NA)) +
  geom_ribbon(data=quants, aes(x=years,ymax=U.50,ymin = L.50),alpha=0.5,fill='blue',color='blue') 
#save_plot(paste0("D:/Github/ICM/Figures/NS_sims/NS_quantiles_climate_starts_at_year_",c.effect,
#                 "_mx_decade_effect_",climate.mx.effect, "_nm_decade_effect_",climate.nm.effect,".png"),p.sims.quants,base_height = 12,base_width = 20)



