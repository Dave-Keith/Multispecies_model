# Here we develop a multi-species model for the North Sea.


# DK Notes

#1: stocks:     # historical estimates, realized lambdas, and trophic levels for each stock in your ecosystem, should be a list
#                 with each stocks being a names list, with the name of the list being a unqiue stock name.
#2 lambdas:     # Lambda estimates from the LTR model run
#3: n.yrs.proj  # How many years into the future we are going to project the stocks
#4: n.sims      # The numbers of simulations to run, keeping low for testing...
#5: er.mn       # Average exploitation rate for the fishery for each stock, set up to be proportional 
#               # Should be the same length as the number of stocks. Defaults to NULL, which is no exploitation
#6: er.sd       # standard deviation of exploitation rate for the fishery for each stock, set up to be proportional 
#               # Should be the same length as the number of stocks. Defaults to NULL, which is no uncertainty
#7: repo.loc    # Location of the Github repo, defaults to "D:/GitHub/Multispecies_model/"


trophic.mod<-function(stocks = NULL,lambdas= NULL,n.yrs.proj = 50, n.sims = 20,
                      catch = list(catch = NULL,er.mn = NULL,er.sd = NULL),
                      repo.loc = "D:/GitHub/Multispecies_model",method = "not_sample")
{
stock.eco <- names(stocks)

library(tidyverse)
library(GGally)
library(cowplot)
library(ggthemes)
library(boot)
# Set the base plot theme
theme_set(theme_few(base_size = 22))
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


########################## Section 2 Parameters ########################## Section 2 Parameters ########################## Section 2 Parameters

# So here we are working to get the 'ecosystem' carrying capacity by looking at the total biomass for the NS stocks we have
# data for over the period of time we have data for all the stocks.
# So here we pull out the data we need to look at total abundance and total biomass in the system by year...
years <- NULL
bm <- NULL
num <- NULL
waa <- NULL
pnm <- NULL
mx <- NULL
ages <- NULL
rem.age <- NULL
for(s in  stock.eco)
{
  years[[s]] <- lambdas[[s]]$year
  #vpa[[s]] <- vpa.tmp[[s]]
  ages[[s]] <- lambdas[[s]]$age.min[1]:lambdas[[s]]$age.max[1]
  num[[s]] <- stocks[[s]] |> collapse::fsubset(type == "Num")
  num[[s]] <- num[[s]] |> collapse::fsubset(age != "tot")
  waa[[s]] <- stocks[[s]]|> collapse::fsubset(type == "WA")
  rem.age[[s]] <- stocks[[s]]|> collapse::fsubset(type == "Catch")
  tl <- rep(unique(stocks[[s]]$TL),nrow(rem.age[[s]]))
  if(s == "ICES-HAWG_NS_Ammodytes_tobianus") waa[[s]]$value <- waa[[s]]$value/1000
  bm[[s]] <- data.frame(Year = num[[s]]$Year,Stock = num[[s]]$Stock,age = num[[s]]$age,
                           bm = num[[s]]$value*waa[[s]]$value,
                           catch.num = rem.age[[s]]$value,
                           catch.bm = rem.age[[s]]$value*waa[[s]]$value,
                           num = num[[s]]$value,
                           trophic = tl,
                           troph.cat = as.character(floor(tl)),
                           Species = num[[s]]$Gen.Spec)
  #Need to clip out the years we don't have biomass data for...
  bm[[s]] <- bm[[s]] |> collapse::fsubset(Year %in% years[[s]])
  #pnm[[s]] <- 1-exp(-lambdas[[s]]$nm.opt)
  #mx[[s]] <- lambdas[[s]]$fecund.opt
  #vpa[[s]] <- lambdas[[s]]$res$est.abund
} # end for(s in  stock.eco)
# Combine the biomass and abundance data into a dataframe
bm.tst <- do.call("rbind",bm)


# Look at the biomass and abundance in the ecosystem
# FIX, about 1% of the catch biomasses are larger than the actual biomass observed, take a look
# and make sure that there isn't something mis-aligned for one of the stocks.
bm.tot <- bm.tst |> collapse::fgroup_by(Stock,Year,trophic,Species,troph.cat) |> 
                    collapse::fsummarize(bm = sum(bm,na.rm=T) + sum(catch.bm,na.rm=T),
                                         num = sum(num,na.rm=T)+ sum(catch.num,na.rm=T))


# The 'ecosystem' biomass and numbers
eco.tot.bm <- bm.tot |> collapse::fgroup_by(Year) |> 
                    collapse::fsummarize(num.eco = sum(num),bm.eco = sum(bm))
# Trophic level biomass and numbers.
trophic.bm <- bm.tot |> collapse::fgroup_by(Year,troph.cat) |> 
                    collapse::fsummarize(num.tl = sum(num),bm.tl = sum(bm))

# All the bm together
tl.eco.bm <- left_join(trophic.bm,eco.tot.bm,by="Year")
tl.eco.bm$prop.bm.tl <- tl.eco.bm$bm.tl/tl.eco.bm$bm.eco
tl.eco.bm$prop.num.tl <- tl.eco.bm$num.tl/tl.eco.bm$num.eco
# They years that are comparable
#tl.eco.bm.comp <- tl.eco.bm |> collapse::fsubset(Year %in% 1990:2014)

# So now take the bm.tot and merge that with the total biomass and the trophic level biomass so we can
# look at what the stock does within it's TL.


# Now we combine the ecosystem results with the stock biomass's
bm.final <- left_join(bm.tot,tl.eco.bm,by=c("Year","troph.cat"))
names(bm.final) <- c("Stock","Year","trophic","species","troph.cat","bm.stock","num.stock","num.tl","bm.tl",'num.eco','bm.eco',
                     'prop.bm.tl','prop.num.tl')
# Get the proportion of the total biomass each stock accounts for
bm.final <- bm.final |> collapse::fmutate(prop.bm.stock.eco = bm.stock/bm.eco,
                                          prop.num.stock.eco = num.stock/num.eco,
                                          prop.bm.stock.tl = bm.stock/bm.tl,
                                          prop.num.stock.tl = num.stock/num.tl)
# Remove 0s from the data
bm.final <- bm.final[bm.final$bm.stock > 0,]
bm.final <- as.data.frame(bm.final)
bm.final$troph.cat <- as.numeric(bm.final$troph.cat)
# This gets the average weight of individuals in each stock, we'll need this later to get an approximate exploitation rate
bm.final$avg.weight <- bm.final$bm.stock/bm.final$num.stock

# Now we subset to the years we have data for all the stocks
what.year <- bm.final |> collapse::fgroup_by(Stock) |> collapse::fsummarize(min = min(Year),
                                                                      max = max(Year))
# The years we have data for all stocks
first.year <- max(what.year$min)
last.year <- min(what.year$max)
n.years <- length(first.year:last.year)
# Now we subset the data to these years
bm.best <- bm.final |> collapse::fsubset(Year %in% first.year:last.year) 
# Now we subset the lambdas
lambdas.tmp <- NULL
for(s in stock.eco)   lambdas.tmp[[s]] <- lambdas[[s]][lambdas[[s]]$year %in% first.year:last.year,] 
lambdas <- lambdas.tmp

# Biomass by trophic level over time
bm.tl.plt <- ggplot(bm.best) + geom_line(aes(x=Year,y=bm.tl,group=as.character(troph.cat),color=as.character(troph.cat))) + 
  scale_color_manual(values = c("blue","red","darkgrey","lightgreen")) + scale_y_log10(name="Biomass")
save_plot(paste0(repo.loc,"/Figures/Biomass_by_trophic_level.png"),bm.tl.plt,base_height = 8,base_width = 11)
# This is real good now...
prop.bm.tl.plt <- ggplot(bm.best) + geom_line(aes(x=Year,y=prop.bm.tl,group=as.character(troph.cat),color=as.character(troph.cat))) + 
  scale_color_manual(values = c("blue","red","darkgrey","lightgreen")) + 
  scale_y_continuous(name="Proportion of Biomass")
save_plot(paste0(repo.loc,"/Figures/Prop_biomass_by_trophic_level.png"),prop.bm.tl.plt,base_height = 8,base_width = 11)

# The biomass for the ecosystem
bm.eco.plt <- ggplot(bm.best) + geom_line(aes(x=Year,y=bm.eco)) + 
                                scale_y_continuous(name="Biomass",limits = c(0,NA))
save_plot(paste0(repo.loc,"/Figures/Biomass_ns_ecosystem.png"),bm.eco.plt,base_height = 8,base_width = 11)


# The 'transfer efficiency' between our trophic levels
tl.3.to.4 <- bm.best$prop.bm.tl[bm.best$troph.cat==4][1:n.years]/bm.best$prop.bm.tl[bm.best$troph.cat==3][1:n.years]
tl.4.to.5 <- bm.best$prop.bm.tl[bm.best$troph.cat==5][1:n.years]/bm.best$prop.bm.tl[bm.best$troph.cat==4][1:n.years]
tl.3.to.5 <- bm.best$prop.bm.tl[bm.best$troph.cat==5][1:n.years]/bm.best$prop.bm.tl[bm.best$troph.cat==3][1:n.years]



# So now we want to look at stock level within a trophic level
# add some colors...
bm.best$color <- "black"
bm.best$color[bm.best$species %in% c("Clupea harengus","Melanogrammus aeglefinus","Gadus morhua")] <- "blue"
bm.best$color[bm.best$species %in% c("Pleuronectes platessa","Solea solea")] <- "green"
bm.best$color[bm.best$species %in% c("Trisopterus esmarkii","Merlangius merlangus")] <- "grey"
bm.best$color[bm.best$species %in% c("Pollachius virens")] <- "orange"

# Put in Species + trophic level
bm.best$spec.tl <- paste(bm.best$species,"(TL = ",bm.best$trophic,")")
# Pull out meta data
meta.dat <- bm.best |> dplyr::group_by(Stock,trophic,species,troph.cat,color,spec.tl) |> filter(row_number() >= (n() ))
meta.dat <- meta.dat[,c("Stock","trophic","species","troph.cat","color","spec.tl")]
meta.dat$troph.cat <- as.numeric(meta.dat$troph.cat)

colors <- distinct(bm.best, spec.tl, color)
pal <- colors$color
names(pal) <- colors$spec.tl
# Another color thing
colors2 <- distinct(bm.best, species, color)
pal2 <- colors$color
names(pal2) <- colors$species



stock.prop.bm.plt <- ggplot(bm.best) + geom_line(aes(x=Year,y=prop.bm.stock.tl,group = Stock,color=spec.tl),linewidth=2) + 
                  facet_wrap(~troph.cat) + guides(colour = guide_legend(nrow = 5)) + theme(legend.position = 'top') +
                  scale_y_log10(name= "Proportion of biomass",n.breaks=10) + scale_x_continuous(name="",labels = c(1990,2000,2010),breaks=c(1990,2000,2010))+
                  scale_color_manual(values=pal)
save_plot(paste0(repo.loc,"/Figures/Prop_Biomass_ns_by_stock.png"),stock.prop.bm.plt,base_height = 8,base_width = 15)

stock.bm.plt <- ggplot(bm.best) + geom_line(aes(x=Year,y=bm.stock,group = Stock,color=spec.tl),linewidth=2) + 
                     facet_wrap(~troph.cat) + scale_x_continuous(name="",labels = c(1990,2000,2010),breaks=c(1990,2000,2010))+
                     scale_y_log10(name = "Biomass",n.breaks=7) + theme(legend.position = 'top') +
                     guides(colour = guide_legend(nrow = 5)) + scale_color_manual(values=pal)
save_plot(paste0(repo.loc,"/Figures/Biomass_ns_by_stock.png"),stock.bm.plt,base_height = 8,base_width = 15)


# So Model 1: You're Basic
# OK, so within a TL each stock has it's own carrying capacity, that is nested within the trophic level carrying capacity
# so if the trophic level is below the carrying capacity each stock gets a bit of that K space for the logistic model. 
# The percentage of the K-space they get is contingent on their historic % of the carrying capacity the stock has had.
# Go with the logistic model too, but I need to build in some uncertainty to the logistic projection
# Initially I'm thinking I'll do (this comment will be outdated by the time you, dear reader, are reading this)
# Step 1: We have a total K for the ecosystem based on past K's, let it vary
# Step 2: We partition that to each trophic level, based on historic splits
# Step 3: We then partition that K to each stock, again based on historic proportion of the K, I wonder how 
#         this will work if a stock is over-fished, the others will be able to fill some of the K-space, but 
#         probably not all of it?
# Step 4: Run the logistic model with the K each stock gets apportioned and we have a trophic level multispecies model

# I think this should work, if we over-fish a stock everyone gets a bit of the free K space (including the overfished stock)
# probably means the trophic level K isn't entirely filled  which could cause problems
# So I can build the ecosystem to have a K that is based on the observed ecosystem biomass history and portion that out
# to each of the trophic levels appropriately, BUT, the population won't necessarily reach that K in any given year, but 
# I guess it should come close. So for base model we have the ecosystem biomass as our K, and 
# then we see if the model is able to get the populations to achieve that K. If we fish
# a bunch of stocks too hard, we have the K, but it'll never reach it. So, assumption that is could be
# a bit problematic, we assume the past ecosystem biomass is K for these stocks, but with this logic, if 
# we overfish we won't reach K, so we are assuming that these stocks were not overfished in totality and thus
# the historic B trend is K (but in reality K was probably > B).... that is if you belive in K in any way shape or form.
 

eco.tot.bm.best <- eco.tot.bm |> collapse::fsubset(Year %in% first.year:last.year)
trophic.bm.best <- trophic.bm |> collapse::fsubset(Year %in% first.year:last.year)

# The correlation in the ecosystem biomass trend, can see this is an AR1
K.cor <- pacf(eco.tot.bm.best$bm.eco,plot=F)
# The cross correlation between the ecosystem biomass trend and the trophic level biomasses
# All correlated, but strongest is unsurprisingly the link between the the ecosystem and the biomass in the
# lowest TL. I suspect this may structurally come out even without explicity building in a lot of
# correlation structure to the models.
K.tl.3.cor <- ccf(eco.tot.bm.best$bm.eco,trophic.bm.best$bm.tl[trophic.bm.best$troph.cat==3],plot = F)
K.tl.4.cor <- ccf(eco.tot.bm.best$bm.eco,trophic.bm.best$bm.tl[trophic.bm.best$troph.cat==4],plot = F)
K.tl.5.cor <- ccf(eco.tot.bm.best$bm.eco,trophic.bm.best$bm.tl[trophic.bm.best$troph.cat==5],plot = F)
# Within trophic levels...
# So these 3 mostly say if the biomass is up one TL, it is up in all TLs, tho there might be some negative between 3 and 4
# at Lag -1 (though that's not quite significant)
tl.3.4.cor <- ccf(trophic.bm.best$bm.tl[trophic.bm.best$troph.cat == 3],trophic.bm.best$bm.tl[trophic.bm.best$troph.cat == 4],plot = F)
tl.3.5.cor <- ccf(trophic.bm.best$bm.tl[trophic.bm.best$troph.cat == 3],trophic.bm.best$bm.tl[trophic.bm.best$troph.cat == 5],plot = F)
tl.4.5.cor <- ccf(trophic.bm.best$bm.tl[trophic.bm.best$troph.cat == 4],trophic.bm.best$bm.tl[trophic.bm.best$troph.cat == 5],plot = F)
# Looking at proportions, need to stew a bit on this because there is necessarily some 
# correlation built into proportions, but what is interesting is that
# the correlation strength is really really high between TL 3 and TL 4, it is weaker (more diffuse really) at TL 5
# and there is no correlation between 4 and 5
tl.3.4.prop.cor <- ccf(bm.best$prop.bm.tl[bm.best$troph.cat == 3][1:n.years],bm.best$prop.bm.tl[bm.best$troph.cat == 4][1:n.years],plot = F)
tl.3.5.prop.cor <- ccf(bm.best$prop.bm.tl[bm.best$troph.cat == 3][1:n.years],bm.best$prop.bm.tl[bm.best$troph.cat == 5][1:n.years],plot = F)
tl.4.5.prop.cor <- ccf(bm.best$prop.bm.tl[bm.best$troph.cat == 4][1:n.years],bm.best$prop.bm.tl[bm.best$troph.cat == 5][1:n.years],plot = F)


# So now really what I need to do first is make a quick simulation that gets me ecosystem K, trophic level K, and stock K
# once I have those then we just run the models :-)
mn.eco.bm <- mean(eco.tot.bm.best$bm.eco)
start.eco.sim <- eco.tot.bm.best$bm.eco[length(eco.tot.bm.best$bm.eco)]
sd.eco.bm <- sd(eco.tot.bm.best$bm.eco)
# trophic level biomass and proportions... for the proportion will probably wanna sample from a beta distro
# So not sure how to do that nicely...
# 


# First get the ecosystem biomass in a correlated time series, there are a whole lot of ways one could do this, this
# is one of many different ideas. I think we could get the 4 and 5 correlations better another way, but
# For a first pass I'm ok with this.
# Ok, duh, use the mean of the time series then the arima gives us the deviations from that mean and we get a nice time series.

# Used for simulations to get good time series for the K for TL3,4, and 5 
tl.3.prop.bm.ts <- bm.best$prop.bm.tl[bm.best$troph.cat==3][1:n.years]
# Extract the frist two components from the pacf to get the two AR components from the model.
tl.3.prop.pacf <- pacf(tl.3.prop.bm.ts,plot = F)
tl.3.prop.bm.lag.1 <- tl.3.prop.pacf$acf[1]
tl.3.prop.bm.lag.2 <- tl.3.prop.pacf$acf[2]
# TL 4 and 5 splits historically
tl.4.5.prop.bm <- bm.best$bm.tl[bm.best$troph.cat == 5][1:n.years]/(bm.best$bm.tl[bm.best$troph.cat == 4][1:n.years]+bm.best$bm.tl[bm.best$troph.cat == 5][1:n.years])
# This is the correlation between 4 and 5
tl.4.5.prop.4.5.bm <- pacf(tl.4.5.prop.bm,plot = F)


troph.levels <- sort(unique(bm.best$troph.cat))

# 
sim.K.stock <- NULL
sim.Ks <- NULL
sim.eco.bm <- NULL
bm.trophic.Ks <- NULL

# Get necessary data on logit scale
tl.3.logit <- logit(tl.3.prop.bm.ts)
tl.4.5.logit <- logit(tl.4.5.prop.bm)

# Starting values for the ecosystem and the proportions, logit needed for arima models with the proportions
start.eco.sim <- eco.tot.bm.best$bm.eco[nrow(eco.tot.bm.best)]
start.tl.3.prop.bm <- tl.3.prop.bm.ts[length(tl.3.prop.bm.ts)]
start.tl.3.logit <- tl.3.logit[length(tl.3.logit)]
start.tl.4.5.prop.bm <- tl.4.5.prop.bm[length(tl.4.5.prop.bm)]
start.tl.4.5.logit <- tl.4.5.logit[length(tl.4.5.logit)]

# Mean values for the trophic levels
mn.tl.3.prop.bm <- mean(tl.3.prop.bm.ts)
mn.tl.3.logit <- mean(tl.3.logit)
mn.tl.4.5.prop.bm <- mean(tl.4.5.prop.bm)
mn.tl.4.5.logit <- mean(tl.4.5.logit)
# We could instead use the most recent year as the mean...
mn.tl.3.prop.bm <- start.tl.3.prop.bm
mn.tl.3.logit <- start.tl.3.logit
mn.tl.4.5.prop.bm <- start.tl.4.5.prop.bm
mn.tl.4.5.logit <- start.tl.4.5.logit

# Difference between starting value an mean (if we use the most recent year as the 'mean', then this is 0)
start.eco.diff = start.eco.sim- mn.eco.bm
start.tl.3.diff <- start.tl.3.logit - mn.tl.3.logit
start.tl.4.5.diff <- start.tl.4.5.logit - mn.tl.4.5.logit

# the standard deviations
sd.tl.3.logit <- sd(tl.3.logit)
sd.tl.4.5.logit <- sd(tl.4.5.logit)
# Lag for the Arima model
tl.4.5.prop.bm.lag.1 <- tl.4.5.prop.4.5.bm$acf[1]
# convert to logit scale for the arima models


for(i in 1:n.sims) 
{
 # The ecosystem K, using the mean of the ecosystem with the correlation observed of the time series.
 # This starts the time series at the last value of the time series, then moves it to the mean value, bam!!  This will be done for each of these arima sims.
  sim.eco.bm[[i]] <- data.frame(bm = c(arima.sim(model =list(ar = K.cor$acf[1]),n = n.yrs.proj,n.start=1,start.innov = start.eco.diff/K.cor$acf[1],
                                                 innov = c(0,rnorm(n.yrs.proj-1,0,sd.eco.bm))) + mn.eco.bm),
                                Years = 1:n.yrs.proj,sim = i) 
  #pacf(sim.eco.bm[[i]]$bm) # looks good

  # So then from my simulated ecosystem I want each trophic level to get it's cut of the biomass, 
  # FIX: I am using the AR2, but I know the start innovation is slightly incorrect, but it make almost no difference for the NS case so I'll stick with it
  # so probably should figure out how to specify that right as it just works by luck here I think, if the difference was larger
  # or correlations different it wouldn't do so well (e.g., it isn't nice for the stock level ones.)
  sim.tl.3.prop.bm <-inv.logit(mn.tl.3.logit + 
                                 arima.sim(model =list(ar = c(tl.3.prop.bm.lag.1,tl.3.prop.bm.lag.2)),n = n.yrs.proj,
                                           n.start =2, start.innov = c(start.tl.3.diff/tl.3.prop.bm.lag.1,start.tl.3.diff/tl.3.prop.bm.lag.1), 
                                           innov = c(0,rnorm(n.yrs.proj-1,0,sd.tl.3.logit))))
  
  bm.sim.3 <- sim.tl.3.prop.bm * sim.eco.bm[[i]]$bm
  # So this is what is left for 3 and 4
  bm.left.4.5<- sim.eco.bm[[i]]$bm - bm.sim.3
  # So then we use the historical split between 4 and 5 can see 5 gets about 1/3-1-5 of 3
   # so then simulate this split
  sim.tl.5.4.prop.bm <- inv.logit(mn.tl.4.5.logit + 
                                    arima.sim(model =list(ar = tl.4.5.prop.bm.lag.1),n = n.yrs.proj,
                                              n.start =1, start.innov = c(start.tl.4.5.diff/tl.4.5.prop.bm.lag.1), 
                                              innov = c(0,rnorm(n.yrs.proj-1,0,sd.tl.3.logit))))
  # And now TL 5 gets this proportion of the 4 and 5 biomass
  bm.sim.5 <- bm.left.4.5 * sim.tl.5.4.prop.bm
  # And TL4 gets the rest, and so the ecosystem biomass is a portion of the whole biomass
  bm.sim.4 <- sim.eco.bm[[i]]$bm - bm.sim.3-bm.sim.5
  
  bm.trophic.Ks[[i]] <- data.frame(Years = rep(1:n.yrs.proj,3), sim =i,
                                   bm.tl = c(bm.sim.3,bm.sim.4,bm.sim.5),troph.cat = as.factor(sort(rep(c(3,4,5),n.yrs.proj))),
                                   bm.eco = rep(sim.eco.bm[[i]]$bm,3))
  bm.trophic.Ks[[i]]$prop.bm.tl <- bm.trophic.Ks[[i]]$bm.tl/bm.trophic.Ks[[i]]$bm.eco
  
  # OK, so now we have the trophic level K values simulated in a 'nice' way. Next how do we partition these to the stocks
  # Give each stock a proportion of the K in it's ecosystem based on their historical cuts of the K, and include the time series correlation in that.
  # I'm going to build in correlation to their K time series (this could 100% be fishery induced correlation), could also put in 
  # cross correlation for species with multiple stocks, but for now, let's just do the AR1/2 thing with this for the proportion of the trophic level 
  # biomass each stock gets.
  
  for(tl in troph.levels)
  {
    tl.stocks <- unique(bm.best$Stock[bm.best$troph.cat==tl])
    n.stock.tl <- length(tl.stocks)
    count =0
    for(s in tl.stocks)
    {
      count = count+1
      # Now get the time series for each stock...
      if(count == 1 ||  n.stock.tl != 2)
      {
        tmp.dat <- bm.best[bm.best$Stock ==s,]
        tmp.cor <- pacf(tmp.dat$prop.bm.tl,plot=F) # Get the correlation, use AR1 and AR2 but no more.
        tmp.cor.lag.1 <- tmp.cor$acf[1]
        #tmp.cor.lag.2 <- tmp.cor$acf[2]
        #tmp.beta <- estBetaParams(mean(tmp.dat$prop.bm.tl),sd(tmp.dat$prop.bm.tl)^2)
        # Logit tranform the proportions and do the ARIMA on the logits
        bm.logit <- logit(tmp.dat$prop.bm.stock.tl)
        start.bm.logit <- bm.logit[length(bm.logit)]
        mn.bm.logit <- mean(bm.logit)
        # Use the most recent bm on logit scale as the 'mean' value for the simulation
        mn.bm.logit <- start.bm.logit
        # And the standard deviation
        sd.bm.logit <- sd(bm.logit)
        diff.bm.logit <- start.bm.logit - mn.bm.logit
        
        # Then backtransform and everything will stay positive! Just using the AR1 term for these
        # FIX: SEE above comment for where I'm using the AR2, here using the AR2 would give some poor starting values
        # So I'm not comfy doing that (it works by luck in the above for the NS IMHO.)
        tmp.prop.bm <- c(inv.logit(arima.sim(model =list(ar = c(tmp.cor.lag.1)),
                                             n.start = 1, start.innov = c(diff.bm.logit/tmp.cor.lag.1),
                                             n = n.yrs.proj,innov = c(0,rnorm(n.yrs.proj-1,0,sd.bm.logit))) + mn.bm.logit))
        sim.Ks[[s]] <- data.frame(Years = 1:n.yrs.proj, sim = i,
                                       Stock = s, troph.cat = tl,
                                       bm.stock = tmp.prop.bm*bm.trophic.Ks[[i]]$bm.tl[bm.trophic.Ks[[i]]$troph.cat==tl])
      } # end the if(count == 1 ||  n.stock.tl != 2)
      # If there are only 2 stocks in a trophic level, then the second stock get the rest of the trophic levels biomass
      
      if(count == 2 & n.stock.tl == 2) 
      {
        sim.Ks[[s]] <-  data.frame(Years = 1:n.yrs.proj, sim=i,
                                   Stock = s, troph.cat = as.numeric(tl),
                                   bm.stock = (1-tmp.prop.bm)*bm.trophic.Ks[[i]]$bm.tl[bm.trophic.Ks[[i]]$troph.cat==tl])
      } # end the case of just 2 stocks
    } # end the stocks loop
  } # end the trophic level loop
  sim.K.stock[[i]] <- do.call("rbind",sim.Ks)
  
} # end the simulation loop

sim.K.stocks <- do.call("rbind",sim.K.stock)
sim.troph.K <- do.call("rbind",bm.trophic.Ks)
sim.eco.K <- do.call("rbind",sim.eco.bm)
# Wrap up the K time series for each simulation
sim.K.stocks$Species <- substr(sim.K.stocks$Stock,14,100)
sim.stock.K.plt <- ggplot(sim.K.stocks[sim.K.stocks$sim==1,]) + geom_line(aes(x=Years,y=bm.stock,group=Species,color=Species),linewidth=2) + 
                             facet_wrap(~troph.cat) + scale_y_log10(name="Biomass") + theme(legend.position = 'top') +
                             guides(colour = guide_legend(nrow = 7))
save_plot(filename = paste0(repo.loc,"/Figures/Simulation_stock_K.png"),sim.stock.K.plt,base_height = 8,base_width = 11)

sim.tl.K.plt <- ggplot(sim.troph.K) + geom_line(aes(x=Years,y=bm.tl,group=as.factor(sim),color=as.factor(sim))) + 
                      facet_wrap(~troph.cat) + theme(legend.position='none') + 
                      scale_y_log10(name="Biomass")

save_plot(filename = paste0(repo.loc,"/Figures/Simulation_trophic_K.png"),sim.tl.K.plt,base_height = 8,base_width = 11)
sim.eco.K.plt <- ggplot(sim.eco.K) + geom_line(aes(x=Years,y=bm,group=as.factor(sim),color=as.factor(sim))) +
                                 theme(legend.position = 'none')
save_plot(filename = paste0(repo.loc,"/Figures/Simulation_eco_K.png"),sim.eco.K.plt,base_height = 8,base_width = 11)

# Comparing TL and ecosystem K going stock by stock with the trophic level and ecosystem K's that I originally made up
# And it's not perfect, but I think for a first pass this work, they keep the characteristics we want in terms of
# correlation and the K's are quite similar to the original ones. For TL3 it is perfect, for 4 and 5 it can be slightly off
# because I aimed to keep the trophic level having the correlation over focusing on getting the K exactly right
# There could be ways to do both I haven't thought of, but think this is ok for now.
 #tst <- sim.K.stocks |> collapse::fgroup_by(troph.cat,Years,sim) |> collapse::fsummarise(tot.bm = sum(bm.stock))
# tst2 <- left_join(tst,sim.troph.K,by=c("Years","troph.cat","sim"))
# tst2 <- tst2 |> collapse::fgroup_by(c('Years','sim')) |>collapse::fmutate(eco.bm.new = sum(tot.bm)) |> as.data.frame()
# # So they definitely differ stock by stock from the original trophic level splits, but most of the time
# # it is within 20% of the original
# tst2$per.diff <- 100*(tst2$tot.bm - tst2$bm.tl) / tst2$bm.tl
# hist(tst2$per.diff[tst2$troph.cat == 5])
# # They do retain the time series characteristics tho
# pacf(tst2$tot.bm[tst2$troph.cat == 3 & tst2$sim == 1])
# # Ecosystem is within 5% way more than 75% of the time
# summary(100*(tst2$eco.bm.new - tst2$bm.eco)/tst2$bm.eco)
# # And still has characteristics we want
# pacf(tst2$eco.bm.new[tst2$sim == 1][1:n.yrs.proj])
# # OK I think we're good to here, have simulated biomass time series :-)
#ggplot(sim.K.stocks |> collapse::fsubset(sim == 1)) + geom_line(aes(x=Years,y=bm.stock,group=Stock,color=troph.cat)) + scale_y_log10()
# ggplot(bm.best) + geom_line(aes(x=Year,y=bm.stock,group=Stock,color=troph.cat))+ scale_y_log10()

# Fix: This is not perfect way to get the past exploitation rates as the removals we have here are in numbers
# Give we have the database with the age specific removals and age specific weights, this
# should be tweaked to use that data. That said, for the moment this should be 'good enough' 
#rem.tst <- do.call("rbind",rem)
#fm.dat <- left_join(bm.best,rem.tst,by=c("Stock",'Year'))
# This is where we go from numbers to a biomass and get an exploitation rate in biomass.
#fm.dat$exploit <- (fm.dat$rem*fm.dat$avg.weight)/fm.dat$bm.stock
# Extract the LTR results.


############### Section 4 Multi-species model of North Sea ############### Section 4 Multi-species model of North Sea ############### Section 4 Multi-species model of North Sea

# Now we have our carrying capacity for each stock and we can get to business and running a model.
# Initialize some things, or maybe no things
res.ts <- NULL
ts.unpack <- NULL
# Get the year range, going from the 'last' year to n.yrs.proj in the future, note this will go 1 year less than your intuition because
# we want n.yrs of data, i.e., 20 years is 2000 to 2019, not 2020... )
#browser()
years <- (last.year+1):(last.year+n.yrs.proj)
# Take the biomass data for the north sea and subset it to the years we have data
bm.mod.yrs <- bm.tst |> collapse::fsubset(Year %in% first.year:last.year)
for(s in Stocks)
bm.start.year <- bm.mod.yrs |> collapse::fsubset(Year == last.year) |> 
                         collapse::fgroup_by(Stock,trophic,troph.cat,Species) |> 
                         collapse::fsummarise(bm.tot = sum(bm,na.rm=T))
# Get the initial ecosystem biomass..
init.eco.bm <- sum(bm.start.year$bm.tot)
init.tl.bm <- bm.start.year |> collapse::fgroup_by(troph.cat) |> collapse::fsummarise(bm.tl = sum(bm.tot))
init.stock.bm <- bm.start.year
# Get the average weight of the fish in the stocks so we can go from biomass to abundance for the model
# FIX: NOT SURE I NEED THIS ANYMORE This could definitely be done more sophisisticatedly!
av.wgt <- bm.best |> collapse::fgroup_by(Stock,troph.cat) |> collapse::fsummarise(mn.wgt = mean(avg.weight,na.rm=T))
# FIX: NOT SURE I NEED THIS ANYMORE Let's try getting the most recent year weight to go from biomass to numbers as average may be somewhat misleading
# So here the idea is that the most recent years 
av.wgt <- bm.best |> dplyr::group_by(Stock,troph.cat) |> filter(row_number() >= (n() ))
av.wgt <- data.frame(Stock = av.wgt$Stock,troph.cat = as.numeric(av.wgt$troph.cat),mn.wgt = av.wgt$avg.weight)
# For some debugging, if still here you can delete I'm sure
#count = 0

# So everything will need to get wrapped up in a simulation loop
for(j in 1:n.sims)
{
  st.time <- Sys.time()
  
  for(t in 1:n.yrs.proj)
  {
    # Get some starting points. These are for the current year
    base.eco.K.tmp <- sim.eco.K |> collapse::fsubset(sim == j & Years ==t)
    base.tl.K.tmp <- sim.troph.K |> collapse::fsubset(sim == j & Years ==t)
    base.stock.K.tmp <- sim.K.stocks |> collapse::fsubset(sim == j & Years ==t)
    # Now get the stock biomass from last year.
    if(t ==1)
    {
      stock.bm.last <- init.stock.bm
      stock.bm.last <- stock.bm.last[order(stock.bm.last$troph.cat),]
      eco.bm.last <- init.eco.bm
      tl.bm.last <- init.tl.bm
    }
    # Then we'll need to get these from the model simulations.
    if(t > 1)
    {
      # Use the handy av.wgt data.frame I made above
      bm.stocks <- data.frame(bm = NA,meta.dat)
      for(s in stock.eco) bm.stocks$bm[bm.stocks$Stock == s] <- res.ts[[s]]$bm[res.ts[[s]]$Years == t-1]
      bm.stocks$bm <- bm.stocks$bm
      stock.bm.last <- bm.stocks
      eco.bm.last <- sum(bm.stocks$bm)
      tl.bm.last <- bm.stocks |> collapse::fgroup_by(troph.cat) |> collapse::fsummarise(bm.tl = sum(bm))
    }  
   # Now we need to figure out what K space is available for each stock within the trophic level.
   # First is our Trophic level above the K we have for it.
    # So this is the K space available in a given trophic level in a year
    base.tl.K.tmp$prop.K.space <- base.tl.K.tmp$bm.tl/tl.bm.last$bm.tl
    # We can then adjust the stock K's by the available K space in each stock
    
    base.stock.K.tmp$K.space <- NA
    base.stock.K.tmp$K.space[base.stock.K.tmp$troph.cat ==3] <- base.stock.K.tmp$bm.stock[base.stock.K.tmp$troph.cat ==3] * 
                                                        (base.tl.K.tmp$prop.K.space[base.tl.K.tmp$troph.cat ==3]-1)
    base.stock.K.tmp$K.space[base.stock.K.tmp$troph.cat ==4] <- base.stock.K.tmp$bm.stock[base.stock.K.tmp$troph.cat ==4] * 
      (base.tl.K.tmp$prop.K.space[base.tl.K.tmp$troph.cat ==4]-1)
    base.stock.K.tmp$K.space[base.stock.K.tmp$troph.cat ==5] <- base.stock.K.tmp$bm.stock[base.stock.K.tmp$troph.cat ==5] * 
      (base.tl.K.tmp$prop.K.space[base.tl.K.tmp$troph.cat ==5]-1)
    
    base.stock.K.tmp$adj.K <- base.stock.K.tmp$bm.stock + base.stock.K.tmp$K.space
    
    # So now I have Carrying Capacities that take up (or lose) any available K space.
    # Now we can convert these to numbers using the historic 'average weight' of the stocks, to avoid complication
    # I'm just using the average of the average weight for each stock...
    base.stock.K.tmp <- left_join(base.stock.K.tmp,av.wgt,by=c("Stock","troph.cat"))
    # And now we can get a K in numbers....
    base.stock.K.tmp$adj.K.num <- base.stock.K.tmp$adj.K/base.stock.K.tmp$mn.wgt
    # Since I have Years and sim recorded, I should just be able to recursivly rbind this...
    if(t ==1 & j == 1) 
    {
      base.stock.K <- base.stock.K.tmp
      base.tl.K <- base.tl.K.tmp
      base.eco.K <- base.eco.K.tmp
    } else {
            base.stock.K <- rbind(base.stock.K,base.stock.K.tmp)
            base.tl.K <- rbind(base.tl.K,base.tl.K.tmp)
            base.eco.K <- rbind(base.eco.K,base.eco.K.tmp)
            } # end the else...

  for(s in stock.eco)
  {
      # Reset samples
      #browser()
      stock.lambdas <- lambdas[[s]] 
      tmp.bm.last <- stock.bm.last |> collapse::fsubset(Stock == s)
      tmp.stock.K <- base.stock.K.tmp |> collapse::fsubset(Stock == s)
      bm.ts.stock <- bm.final[bm.final$Stock == s & bm.final$Year %in% first.year:last.year,]  
      #a=s
      #browser()
      # Now get the final year bm
      if(t == 1) 
      { 
        
        bm.start <- bm.ts.stock$bm.stock[bm.ts.stock$Year == last.year]
        res.ts[[s]] <- data.frame(bm = bm.start,removals = NA,ex.rate = NA,
                                  Stock = s,sim= j,lambda = NA,Years=t-1,
                                  troph.cat = as.numeric(unique(bm.tot$troph.cat[bm.tot$Stock ==s])),
                                  K.bm = NA)
        
      } else{ bm.start <- res.ts[[s]]$bm[res.ts[[s]]$Years == t-1]}
      
      # Sort out which of the years are low or high bm
      # I'm using 0.5 as the cut off, other options are valid (0.4 is my fav...)
      low.vs.high <- 0.5
      low.vs.high.bm <- low.vs.high * max(bm.ts.stock$bm.stock)
      # Have to drop the final year because we don't have a lambda estimate for the final year
      low.bm.years <- which(bm.ts.stock$bm.stock[-nrow(bm.ts.stock)] < low.vs.high.bm)
      if(length(low.bm.years) == 0) low.years <- F else low.years <- T
      high.bm.years <- which(bm.ts.stock$bm.stock[-nrow(bm.ts.stock)] >= low.vs.high.bm)
      # ANd get what our carrying capacity is at this moment
      cur.K <- tmp.stock.K$adj.K
      
      # So first, get a sample 
      method <- "not_sample"
      if(method == 'sample')
      {
        # Pick one of these to sample if that's how we want to roll, if we have low biomass years (as
        # defined by our cut off low.vs.high)
        if(low.years == T)
        {
        if(bm.start < low.vs.high.bm) samp <- sample(low.bm.years,1)
        if(bm.start >= low.vs.high.bm) samp <- sample(high.bm.years,1)
        } # end If we have low years
        if(low.years == F) samp <- sample(high.bm.years,1)
        # The simple way to do it is just to sample from the natural mortality distribution
        # Now using the right lambda, go look at trends from the stocks that are declining to see what's up there.
        lam.samp <- stock.lambdas$lam.no.fish[samp] # Get the sample years.  
      } # end the sample method.
      
      # Or do it the fun way...
      if(method != "sample")
      {
      
        # The fun way to do it is to do something multivariate! Note these are instantaneous now!!
        if(bm.start < low.vs.high.bm) 
        {
          if(length(low.bm.years) >0)
          {
          lam.mn <- mean(stock.lambdas$lam.no.fish[low.bm.years],na.rm=T)
          lam.sd <- sd(log(stock.lambdas$lam.no.fish[low.bm.years]),na.rm=T)
          if(length(low.bm.years) == 1) lam.sd <- 0.2 # In case there is just one low biomass year
          }
          if(length(low.bm.years) ==0)
          {
            lam.mn <- mean(stock.lambdas$lam.no.fish,na.rm=T)
            lam.sd <- sd(log(stock.lambdas$lam.no.fish),na.rm=T)
          }
          lam.samp <- rlnorm(1,log(lam.mn),lam.sd)
          #if(is.na(lam.samp)) browser()
        } # end if(bm.start < low.vs.high.bm) 
        
        if(bm.start >= low.vs.high.bm & bm.start < cur.K) 
        {
          lam.mn <- mean(stock.lambdas$lam.no.fish[high.bm.years],na.rm=T)
          lam.sd <- sd(log(stock.lambdas$lam.no.fish[high.bm.years]),na.rm=T)
          lam.samp <- rlnorm(1,log(lam.mn),lam.sd)
          
        } # end if(bm.start < low.vs.high.bm) 
        
      } # end if(method != "sample")
      #if(is.na(lam.samp)) browser()
      while(is.na(lam.samp)) lam.samp <- rlnorm(1,log(lam.mn),lam.sd)
      # Final one, if we are above the K, we are just doing the high biomass scenario for now.
      # Solution, sample from the lambdas at high biomass, but only take lambdas that are <= 1
      if(bm.start >= cur.K) 
      {
        lam.mn <- mean(stock.lambdas$lam.no.fish[high.bm.years],na.rm=T)
        lam.sd <- sd(log(stock.lambdas$lam.no.fish[high.bm.years]),na.rm=T)
        lam.samp <- rlnorm(1,log(lam.mn),lam.sd)
        #if(is.na(lam.samp)) browser()
        while(lam.samp >1) lam.samp <- rlnorm(1,log(lam.mn),lam.sd)
      }
      #while(is.na(lam.samp)) lam.samp <- rlnorm(1,log(lam.mn),lam.sd)
      
# Simple way to include removals

if(is.null(catch$er.mn)) {ex.rate = 0; ex.sd = 0}
if(!is.null(catch$er.mn))
{
  er.mn <- catch$er.mn
  if(!is.null(catch$er.sd)) er.sd <- catch$er.sd
  if(is.null(catch$er.sd)) er.sd <- data.frame(er.sd=0,Stock =s)
  ex.rate = er.mn$ex.mn[er.mn$Stock == s]
  ex.sd = er.sd$ex.sd[er.sd$Stock == s]
} # end if(!is.null(catch$er.mn))
if(is.null(catch$catch)) 
{
  #print(ex.rate)
  # Convert proportion to F
  ex.rate <- -log(1-ex.rate)
  er <- rlnorm(1,log(ex.rate),ex.sd)
  # Go from F back to proportion
  er <- 1-exp(-er)
  removals <- bm.start*er
} # end if(is.null(catch$catch)) 
      
# If you have a catch estimate for the stock
  if(!is.null(catch$catch))
  {
    limit.er <- 0.4
    removals.tmp <- catch$catch$catch[catch$catch$Stock ==s]
    er <- removals.tmp/(bm.start+removals.tmp)
    if(er > limit.er) 
    {
      removals.tmp <- limit.er*bm.start
      er <- limit.er
    }
    removals <- removals.tmp
  } # end if(!is.null(catch$catch))


#print(er)

#lam.samp <- rlnorm(1,log(1),0.1)
#tst.res <- (lam.samp)*bm.start - removals
# We want to grow after removals because otherwise we can get negative values given exploitation was
# calculated using the initial biomass
tst.res <- lam.samp*(bm.start - removals)
     
# Because I'm only doing this one year at a time, there's something in here I need to mess around with to get the output tidy...

res.ts[[s]] <- rbind(res.ts[[s]] ,data.frame(bm = tst.res,removals =removals,ex.rate = er,
                                             Stock = s,sim= j,lambda = lam.samp,Years=t,
                                             troph.cat = as.numeric(unique(bm.tot$troph.cat[bm.tot$Stock ==s])),
                                             K.bm = tmp.stock.K$adj.K))
#if(tst$r$r[1] > 14) stop("WTF")
#abund.new[[s]] <- data.frame(abund = tst$Pop$abund[2])
#res.r[[s]] <- data.frame(tst$r[1,-2],stock=s,sim=j)

  } # end stock loop
  #res.tst[[t]] <- do.call("rbind",res.ts)
  } # end the t looping through each year.
  
#ggplot(base.stock.K) + geom_line(aes(x= Years,y=bm.stock,group=Stock,color=Stock)) + facet_wrap(~troph.cat) + scale_y_log10()
  
  # Unpack the results
  ts.unpack[[j]] <- do.call('rbind',res.ts)
  
#ggplot(ts.unpack[[j]]) + geom_line(aes(x= Years,y=abund,group=Stock,color=Stock)) + facet_wrap(~troph.cat) + scale_y_log10()
  
  
  # Pop a note when done each simulation
  timer <- Sys.time() - st.time
  print(paste("Simulation ", j))
  print(signif(timer,digits=2))
  
} # end n.sims

# Unpack all the results.
ts.final <- do.call("rbind",ts.unpack)

#ggplot(ts.final) + geom_line(aes(x= Years,y=abund,group=sim,color=sim)) + facet_wrap(~Stock) + scale_y_log10()


#ts.final$fm <- ts.final$removals/ts.final$abund
#av.wgt$troph.cat <- as.numeric(av.wgt$troph.cat)
ts.final <- left_join(ts.final,av.wgt,by=c("Stock","troph.cat"))
#ts.final$biomass <- ts.final$abund*ts.final$mn.wgt
#r.final <- do.call("rbind",r.unpack)

quants <- ts.final |>  collapse::fgroup_by(Years,Stock,troph.cat) |> collapse::fsummarize(L.50 = quantile(bm,probs=c(0.25),na.rm=T),
                                                                          med = median(bm,na.rm=T),
                                                                          U.50 = quantile(bm,probs=c(0.75),na.rm=T))#,
                                                                          #fml.50 = quantile(fm,probs=c(0.25),na.rm=T),
                                                                          #fm = median(fm,na.rm=T),
                                                                          #fmu.50 = quantile(fm,probs=c(0.75),na.rm=T))

ts.final <- left_join(ts.final,meta.dat,by = c("Stock","troph.cat"))
quants <- left_join(quants,meta.dat,by = c("Stock","troph.cat"))


# If happy save the 2 objects
saveRDS(object = ts.final,file = paste0(repo.loc,"/Results/NS_projections_",n.sims,"_sims_",min(years),"_to_",max(years),
                                         "_time_series_projections.Rds"))

saveRDS(object = quants,file = paste0(repo.loc,"/Results/NS_projections_",n.sims,"_sims_",min(years),"_to_",max(years),
                                       "_time_series_quantiles.Rds"))
# 
# saveRDS(object = r.final,file = paste0(repo.loc,"/Results/NS_projections_",n.sims,"_sims_",min(years),"_to_",max(years),
#                                        "_r_projections.Rds"))


# Two simple plots. 
p.sims <- ggplot(ts.final ) + geom_line(aes(x=Years,y=bm,group = sim,color=sim),alpha=0.8) +
  facet_wrap(~Stock,scales = 'free_y') + 
  scale_x_continuous(breaks = seq(1,50,by=49),labels=c(2015,2065)) +
  scale_y_log10(name = "Biomass") + 
  theme(legend.position = 'none') 

save_plot(paste0(repo.loc,"/Figures/biomass_trends.png"),p.sims,base_height = 12,base_width = 24)



p.sims.quants <- ggplot(quants) + geom_line(aes(x=Years,y=med,group=Stock,color=spec.tl)) + 
  facet_wrap(~troph.cat,scales = 'free_y') +  scale_y_log10(name="Biomass") +   theme(legend.position = 'top') +
  guides(colour = guide_legend(nrow = 7)) +
  scale_x_continuous(breaks = seq(1,50,by=49),labels=c(2015,2065)) 
  #geom_ribbon(data=quants, aes(x=Years,ymax=U.50,ymin = L.50),alpha=0.5,fill='blue',color='blue') 
save_plot(paste0(repo.loc,"/Figures/Quantile_abundance_trends.png"),p.sims.quants,base_height = 8,base_width = 16)


colors <- distinct(bm.best, spec.tl, color)
pal <- colors$color
names(pal) <- colors$spec.tl

p.sims.quants <- ggplot(quants) + geom_line(aes(x=Years,y=med,group=Stock,color=spec.tl),linewidth=2) + 
  facet_wrap(~troph.cat,scales = 'free_y') +  scale_y_log10(name="Biomass") +   theme(legend.position = 'top') +
  guides(colour = guide_legend(nrow = 5)) + scale_color_manual(values=pal) +
  scale_x_continuous(name="",breaks = seq(1,50,by=49),labels=c(2015,2065)) 
#geom_ribbon(data=quants, aes(x=Years,ymax=U.50,ymin = L.50),alpha=0.5,fill='blue',color='blue') 
save_plot(paste0(repo.loc,"/Figures/Quantile_biomass_trends.png"),p.sims.quants,base_height = 8,base_width = 16)



# ggplot(bm.best) + geom_line(aes(x=Year,y=prop.bm.stock.tl,group = Stock,color=spec.tl),linewidth=2) + 
#   facet_wrap(~troph.cat) + guides(colour = guide_legend(nrow = 5)) + theme(legend.position = 'top') +
#   scale_y_log10(name= "Proportion of biomass",n.breaks=10) + scale_x_continuous(name="",labels = c(1990,2000,2010),breaks=c(1990,2000,2010))+
#   scale_color_manual(values=pal)

return(list(sim.quantiles = quants,
            sim.ts = ts.final,
            past.bm = bm.best,
            sim.K.stocks = sim.K.stocks,
            sim.troph.K = sim.troph.K,
            sim.eco.K = sim.eco.K))
} # end function