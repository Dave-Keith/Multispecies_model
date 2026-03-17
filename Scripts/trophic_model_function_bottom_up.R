# Here we develop a multi-species model for the North Sea.


# DK Notes

#1: stocks:     # historical estimates, realized lambdas, and trophic levels for each stock in your ecosystem, should be a list
#                 with each stocks being a names list, with the name of the list being a unqiue stock name.
#2 lambdas:     # Lambda estimates from the LTR model run
#3: n.yrs.proj  # How many years into the future we are going to project the stocks
#4: n.sims      # The numbers of simulations to run, keeping low for testing...
#5: explolit    # List of the stock, remmovals and exploitation rate history of the stock
#6: manage      # Management scenario to test. If NULL (default) management # If you don't do anything, the default is to fish it using the 
                # expoitation from the time series as the RR, the LRP is 40% of median, and the USR is the median
                # No fishing below the LRP, and fishing turns down linearily between the USR and LRP.           # Should be the same length as the number of stocks. Defaults to NULL, which is no uncertainty
#7: repo.loc    # Location of the Github repo, defaults to "D:/GitHub/Multispecies_model/"




trophic.mod<-function(stocks = NULL,lambdas= NULL,n.yrs.proj = 50, n.sims = 20,
                      exploit = NULL,
                      repo.loc = "D:/GitHub/Multispecies_model",method = "not_sample")
{
stock.eco <- names(stocks)

library(tidyverse)
library(GGally)
library(cowplot)
library(ggthemes)
library(boot)

if(is.null(exploit)) stop("Nope, not happening, you need exploitation data")
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
#ages <- NULL
rem.age <- NULL
for(s in  stock.eco)
{
  years[[s]] <- lambdas[[s]]$year
  #vpa[[s]] <- vpa.tmp[[s]]
  #ages[[s]] <- lambdas[[s]]$age.min[1]:lambdas[[s]]$age.max[1]
  num[[s]] <- stocks[[s]] |> collapse::fsubset(type == "Num")
  num[[s]] <- num[[s]] |> collapse::fsubset(age != "tot")
  waa[[s]] <- stocks[[s]]|> collapse::fsubset(type == "WA")
  #rem.age[[s]] <- stocks[[s]]|> collapse::fsubset(type == "catch")
  tl <- rep(unique(stocks[[s]]$TL),nrow(waa[[s]]))
  tc <- rep(unique(stocks[[s]]$troph.cat),nrow(waa[[s]]))
  
  #if(s == "ICES-HAWG_NS_Ammodytes_tobianus") waa[[s]]$value <- waa[[s]]$value/1000
  #browser()
  bm[[s]] <- data.frame(Year = num[[s]]$Year,Stock = num[[s]]$Stock,age = num[[s]]$age,
                           bm = num[[s]]$value*waa[[s]]$value,
                           #catch.num = rem.age[[s]]$value,
                           #catch.bm = rem.age[[s]]$value*waa[[s]]$value,
                           num = num[[s]]$value,
                           trophic = tl,
                           troph.cat = tc,
                           Species = num[[s]]$Gen.Spec,
                           Stock.short = num[[s]]$Stock.short,
                           common = num[[s]]$common)
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
bm.tot <- bm.tst |> collapse::fgroup_by(Stock,Year,trophic,Species,troph.cat,Stock.short,common) |> 
                    collapse::fsummarize(bm = sum(bm,na.rm=T), #+ sum(catch.bm,na.rm=T),
                                         num = sum(num,na.rm=T))# sum(catch.num,na.rm=T),
                                         #catch = sum(catch.bm,na.rm=T))

#bm.tot$er <- bm.tot$catch/bm.tot$bm # I could add catch back into the denominator of this equation...


# The 'ecosystem' biomass and numbers
eco.tot.bm <- bm.tot |> collapse::fgroup_by(Year) |> 
                    collapse::fsummarize(num.eco = sum(num),bm.eco = sum(bm))
# Trophic level biomass and numbers.
trophic.bm <- bm.tot |> collapse::fgroup_by(Year,troph.cat) |> 
                    collapse::fsummarize(num.tl = sum(num),bm.tl = sum(bm))
trophic.bm$troph.cat <- factor(trophic.bm$troph.cat,levels = c("≤ 4.0","4.1-4.9","≥ 5.0"))

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
names(bm.final) <- c("Stock","Year","trophic","species","troph.cat","Stock.short",'common',"bm.stock","num.stock","num.tl","bm.tl",'num.eco','bm.eco',
                     'prop.bm.tl','prop.num.tl')
# Get the proportion of the total biomass each stock accounts for
bm.final <- bm.final |> collapse::fmutate(prop.bm.stock.eco = bm.stock/bm.eco,
                                          prop.num.stock.eco = num.stock/num.eco,
                                          prop.bm.stock.tl = bm.stock/bm.tl,
                                          prop.num.stock.tl = num.stock/num.tl)
# Remove 0s from the data
bm.final <- bm.final[bm.final$bm.stock > 0,]
bm.final <- as.data.frame(bm.final)
bm.final$troph.cat <- factor(bm.final$troph.cat,levels = c("≤ 4.0","4.1-4.9","≥ 5.0"))
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
bm.tl.plt <- ggplot(bm.best) + geom_line(aes(x=Year,y=bm.tl,group=troph.cat,color=troph.cat)) + 
  scale_color_manual(values = c("blue","red","darkgrey","lightgreen")) + scale_y_log10(name="Biomass Pool (tonnes)") + theme(legend.title = element_blank()) 
save_plot(paste0(repo.loc,"/Figures/BU/Historic_Biomass_by_trophic_level.png"),bm.tl.plt,base_height = 8,base_width = 11)
# This is real good now...
prop.bm.tl.plt <- ggplot(bm.best) + geom_line(aes(x=Year,y=prop.bm.tl,group=troph.cat,color=troph.cat)) + 
  scale_color_manual(values = c("blue","red","darkgrey","lightgreen")) + 
  scale_y_continuous(name="Proportion of Biomass Pool") + theme(legend.title = element_blank()) 
save_plot(paste0(repo.loc,"/Figures/BU/Historic_Prop_biomass_by_trophic_level.png"),prop.bm.tl.plt,base_height = 8,base_width = 11)

# The biomass for the ecosystem
bm.eco.plt <- ggplot(bm.best) + geom_line(aes(x=Year,y=bm.eco)) + 
                                scale_y_continuous(name="Biomass Pool (tonnes)",limits = c(0,NA))
save_plot(paste0(repo.loc,"/Figures/BU/Historic_Biomass_ns_ecosystem.png"),bm.eco.plt,base_height = 8,base_width = 11)

# The 'transfer efficiency' between our trophic levels
tl.3.to.4 <- bm.best$prop.bm.tl[bm.best$troph.cat=="4.1-4.9"][1:n.years]/bm.best$prop.bm.tl[bm.best$troph.cat=="≤ 4.0"][1:n.years]
tl.4.to.5 <- bm.best$prop.bm.tl[bm.best$troph.cat=="≥ 5.0"][1:n.years]/bm.best$prop.bm.tl[bm.best$troph.cat=="4.1-4.9"][1:n.years]
tl.3.to.5 <- bm.best$prop.bm.tl[bm.best$troph.cat=="≥ 5.0"][1:n.years]/bm.best$prop.bm.tl[bm.best$troph.cat=="≤ 4.0"][1:n.years]

#browser()

# So now we want to look at stock level within a trophic level
# add some colors... this is so clunky....
tl3s <- unique(bm.best$species[bm.best$troph.cat=="≤ 4.0"])
tl4s <- unique(bm.best$species[bm.best$troph.cat=="4.1-4.9"])
tl5s <- unique(bm.best$species[bm.best$troph.cat=="≥ 5.0"])
bm.best$color <- "black"
count=1
for(c in tl3s) 
{
  if(count == 2) bm.best$color[bm.best$species == tl3s[count]] <- "blue"
  if(count == 3) bm.best$color[bm.best$species == tl3s[count]] <- "green"
  if(count == 4) bm.best$color[bm.best$species == tl3s[count]] <- "grey"
  if(count == 5) bm.best$color[bm.best$species == tl3s[count]] <- "orange"
  if(count == 6) bm.best$color[bm.best$species == tl3s[count]] <- "firebrick2"
  count = count + 1
}
# TL4 colors
count=1
for(c in tl4s) 
{
  if(count == 2) bm.best$color[bm.best$species == tl4s[count]] <- "blue"
  if(count == 3) bm.best$color[bm.best$species == tl4s[count]] <- "green"
  if(count == 4) bm.best$color[bm.best$species == tl4s[count]] <- "grey"
  if(count == 5) bm.best$color[bm.best$species == tl4s[count]] <- "orange"
  if(count == 6) bm.best$color[bm.best$species == tl4s[count]] <- "firebrick2"
  count = count + 1
}
# TL5 colors
count=1
for(c in tl5s) 
{
  if(count == 2) bm.best$color[bm.best$species == tl5s[count]] <- "blue"
  if(count == 3) bm.best$color[bm.best$species == tl5s[count]] <- "green"
  if(count == 4) bm.best$color[bm.best$species == tl5s[count]] <- "grey"
  if(count == 5) bm.best$color[bm.best$species == tl5s[count]] <- "orange"
  if(count == 6) bm.best$color[bm.best$species == tl5s[count]] <- "firebrick2"
  count = count + 1
}
#browser()
# Put in Species + trophic level
bm.best$spec.tl <- paste(bm.best$species,"(TL is ",bm.best$trophic,")")
# Pull out meta data
meta.dat <- bm.best |> dplyr::group_by(Stock,trophic,species,troph.cat,color,spec.tl) |> filter(row_number() >= (n() ))
meta.dat <- meta.dat[,c("Stock","trophic","species","troph.cat","color","spec.tl","Stock.short",'common')]
#meta.dat$troph.cat <- as.numeric(meta.dat$troph.cat)

colors <- distinct(bm.best, spec.tl, color)
pal <- colors$color
names(pal) <- colors$spec.tl

stock.prop.bm.plt <- ggplot(bm.best) + geom_line(aes(x=Year,y=prop.bm.stock.tl,group = Stock,color=spec.tl),linewidth=2) + 
                  facet_wrap(~troph.cat) + guides(colour = guide_legend(nrow = 5)) + theme(legend.position = 'top',legend.title = element_blank()) +
                  scale_y_log10(name= "Proportion of Biomass Pool",n.breaks=10) + scale_x_continuous(name="",labels = c(1990,2000,2010),breaks=c(1990,2000,2010))+
                  scale_color_manual(values=pal)
save_plot(paste0(repo.loc,"/Figures/BU/Historic_Prop_Biomass_ns_by_stock.png"),stock.prop.bm.plt,base_height = 8,base_width = 15)

stock.bm.plt <- ggplot(bm.best) + geom_line(aes(x=Year,y=bm.stock,group = Stock,color=spec.tl),linewidth=2) + 
                     facet_wrap(~troph.cat) + scale_x_continuous(name="",labels = c(1990,2000,2010),breaks=c(1990,2000,2010))+
                     scale_y_log10(name = "Biomass Pool (tonnes)",n.breaks=7) + theme(legend.position = 'top',legend.title = element_blank()) +
                     guides(colour = guide_legend(nrow = 5)) + scale_color_manual(values=pal)
save_plot(paste0(repo.loc,"/Figures/BU/Historic_Biomass_ns_by_stock.png"),stock.bm.plt,base_height = 8,base_width = 15)


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
K.tl.3.cor <- ccf(eco.tot.bm.best$bm.eco,trophic.bm.best$bm.tl[trophic.bm.best$troph.cat=="≤ 4.0"],plot = F)
K.tl.4.cor <- ccf(eco.tot.bm.best$bm.eco,trophic.bm.best$bm.tl[trophic.bm.best$troph.cat=="4.1-4.9"],plot = F)
K.tl.5.cor <- ccf(eco.tot.bm.best$bm.eco,trophic.bm.best$bm.tl[trophic.bm.best$troph.cat=="≥ 5.0"],plot = F)
# Within trophic levels...
# So these 3 mostly say if the biomass is up one TL, it is up in all TLs, tho there might be some negative between 3 and 4
# at Lag -1 (though that's not quite significant)
tl.3.4.cor <- ccf(trophic.bm.best$bm.tl[trophic.bm.best$troph.cat =="≤ 4.0"],trophic.bm.best$bm.tl[trophic.bm.best$troph.cat =="4.1-4.9"],plot = F)
tl.3.5.cor <- ccf(trophic.bm.best$bm.tl[trophic.bm.best$troph.cat =="≤ 4.0"],trophic.bm.best$bm.tl[trophic.bm.best$troph.cat =="≥ 5.0"],plot = F)
tl.4.5.cor <- ccf(trophic.bm.best$bm.tl[trophic.bm.best$troph.cat =="4.1-4.9"],trophic.bm.best$bm.tl[trophic.bm.best$troph.cat =="≥ 5.0"],plot = F)
# Looking at proportions, need to stew a bit on this because there is necessarily some 
# correlation built into proportions, but what is interesting is that
# the correlation strength is really really high between TL 3 and TL 4, it is weaker (more diffuse really) at TL 5
# and there is no correlation between 4 and 5
tl.3.4.prop.cor <- ccf(bm.best$prop.bm.tl[bm.best$troph.cat =="≤ 4.0"][1:n.years],bm.best$prop.bm.tl[bm.best$troph.cat =="4.1-4.9"][1:n.years],plot = F)
tl.3.5.prop.cor <- ccf(bm.best$prop.bm.tl[bm.best$troph.cat =="≤ 4.0"][1:n.years],bm.best$prop.bm.tl[bm.best$troph.cat =="≥ 5.0"][1:n.years],plot = F)
tl.4.5.prop.cor <- ccf(bm.best$prop.bm.tl[bm.best$troph.cat =="4.1-4.9"][1:n.years],bm.best$prop.bm.tl[bm.best$troph.cat =="≥ 5.0"][1:n.years],plot = F)


# So now really what I need to do first is make a quick simulation that gets me ecosystem K, trophic level K, and stock K
# once I have those then we just run the models :-)
mn.eco.bm <- mean(eco.tot.bm.best$bm.eco)
#start.eco.sim <- eco.tot.bm.best$bm.eco[length(eco.tot.bm.best$bm.eco)]
sd.eco.bm <- sd(eco.tot.bm.best$bm.eco)
# trophic level biomass and proportions... for the proportion will probably wanna sample from a beta distro
# So not sure how to do that nicely...
# 


# First get the ecosystem biomass in a correlated time series, there are a whole lot of ways one could do this, this
# is one of many different ideas. I think we could get the 4 and 5 correlations better another way, but
# For a first pass I'm ok with this.
# Ok, duh, use the mean of the time series then the arima gives us the deviations from that mean and we get a nice time series.

# Used for simulations to get good time series for the K for TL3,4, and 5 
tl.3.prop.bm.ts <- bm.best$prop.bm.tl[bm.best$troph.cat=="≤ 4.0"][1:n.years]
# Extract the frist two components from the pacf to get the two AR components from the model.
tl.3.prop.pacf <- pacf(tl.3.prop.bm.ts,plot = F)
tl.3.prop.bm.lag.1 <- tl.3.prop.pacf$acf[1]
tl.3.prop.bm.lag.2 <- tl.3.prop.pacf$acf[2]
# TL 4 and 5 splits historically
tl.4.5.prop.bm <- bm.best$bm.tl[bm.best$troph.cat =="≥ 5.0"][1:n.years]/(bm.best$bm.tl[bm.best$troph.cat =="4.1-4.9"][1:n.years]+bm.best$bm.tl[bm.best$troph.cat =="≥ 5.0"][1:n.years])
# This is the correlation between 4 and 5
tl.4.5.prop.4.5.bm <- pacf(tl.4.5.prop.bm,plot = F)


troph.levels <- sort(unique(bm.best$troph.cat))

# 
sim.K.stock.lst <- NULL
sim.Ks <- NULL
sim.eco.bm <- NULL
bm.trophic.Ks <- NULL
bm.sim.5 <- NULL
bm.sim.4 <- NULL
bm.sim.3 <- NULL
sim.tl.3.to.4 <- NULL
sim.tl.4.to.5 <- NULL

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
#mn.tl.3.prop.bm <- start.tl.3.prop.bm
#mn.tl.3.logit <- start.tl.3.logit
#mn.tl.4.5.prop.bm <- start.tl.4.5.prop.bm
#mn.tl.4.5.logit <- start.tl.4.5.logit

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
#browser()
# NEED TO SORT OUT THE INDEXING ON THE Bm.sim.tl objects and how I use these later on!!
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
  # sim.tl.3.prop.bm <-inv.logit(mn.tl.3.logit + 
  #                                arima.sim(model =list(ar = c(tl.3.prop.bm.lag.1,tl.3.prop.bm.lag.2)),n = n.yrs.proj,
  #                                          n.start =2, start.innov = c(start.tl.3.diff/tl.3.prop.bm.lag.1,start.tl.3.diff/tl.3.prop.bm.lag.1), 
  #                                          innov = c(0,rnorm(n.yrs.proj-1,0,sd.tl.3.logit))))
  # 
  sim.tl.3.prop.bm <-inv.logit(mn.tl.3.logit + 
                                 arima.sim(model =list(ar = c(tl.3.prop.bm.lag.1)),n = n.yrs.proj,
                                           n.start =1, start.innov = c(start.tl.3.diff/tl.3.prop.bm.lag.1), 
                                           innov = c(0,rnorm(n.yrs.proj-1,0,sd.tl.3.logit))))
  
  bm.sim.3[[i]] <- sim.tl.3.prop.bm * sim.eco.bm[[i]]$bm
  # So this is what is left for TL 4 and 5
  bm.left.4.5<- sim.eco.bm[[i]]$bm - bm.sim.3[[i]]
  # So then we use the historical split between 4 and 5 can see 5 gets about 1/3-1-5 of 3
   # so then simulate this split
  sim.tl.5.4.prop.bm <- inv.logit(mn.tl.4.5.logit + 
                                    arima.sim(model =list(ar = tl.4.5.prop.bm.lag.1),n = n.yrs.proj,
                                              n.start =1, start.innov = c(start.tl.4.5.diff/tl.4.5.prop.bm.lag.1), 
                                              innov = c(0,rnorm(n.yrs.proj-1,0,sd.tl.3.logit))))
  # And now TL 5 gets this proportion of the 4 and 5 biomass
  bm.sim.5[[i]] <- bm.left.4.5 * sim.tl.5.4.prop.bm
  # And TL4 gets the rest, and so the ecosystem biomass is a portion of the whole biomass
  bm.sim.4[[i]] <- sim.eco.bm[[i]]$bm - bm.sim.3[[i]]-bm.sim.5[[i]]
  #browser()
  bm.trophic.Ks[[i]] <- data.frame(Years = rep(1:n.yrs.proj,3), sim =i,
                                   bm.tl = c(bm.sim.3[[i]],bm.sim.4[[i]],bm.sim.5[[i]]),troph.cat = sort(rep(troph.levels,n.yrs.proj)),
                                   bm.eco = rep(sim.eco.bm[[i]]$bm,3))
  bm.trophic.Ks[[i]]$prop.bm.tl <- bm.trophic.Ks[[i]]$bm.tl/bm.trophic.Ks[[i]]$bm.eco
  
  sim.tl.3.to.4[[i]] <- bm.sim.4[[i]]/bm.sim.3[[i]]
  sim.tl.4.to.5[[i]] <- bm.sim.5[[i]]/bm.sim.4[[i]]
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
      
        tmp.dat <- bm.best[bm.best$Stock ==s,]
        tmp.cor <- pacf(tmp.dat$prop.bm.tl,plot=F) # Get the correlation, use AR1 and AR2 but no more.
        tmp.cor.lag.1 <- tmp.cor$acf[1]
        #tmp.cor.lag.2 <- tmp.cor$acf[2]
        #tmp.beta <- estBetaParams(mean(tmp.dat$prop.bm.tl),sd(tmp.dat$prop.bm.tl)^2)
        # Logit tranform the proportions and do the ARIMA on the logits
        bm.logit <- logit(tmp.dat$prop.bm.stock.tl)
        start.bm.logit <- bm.logit[length(bm.logit)]
        mn.bm.logit <- mean(bm.logit)
        #DK Note, I thought about using the most recent bm on logit scale as the 'mean' value for the simulation to 
        # start where we finished, but I don't like that behaviour in TL3, so going to use the mean which means the stocks
        # will want to go back to an average value of the 'sharing' of biomass. Shit-canning this will make
        # some of the below unnecessarily complicated.
        #mn.bm.logit <- start.bm.logit
        # And the standard deviation
        sd.bm.logit <- sd(bm.logit)
        diff.bm.logit <- start.bm.logit - mn.bm.logit
        
        # Then backtransform and everything will stay positive! Just using the AR1 term for these
        # FIX: SEE above comment for where I'm using the AR2, here using the AR2 would give some poor starting values
        # So I'm not comfy doing that (it works by luck in the above for the NS IMHO.)
        # These can sum to be > 1 within a trophic level, thanks for catching that HB
        # so we'll need to tidy that up
        tmp.prop.bm <- c(inv.logit(arima.sim(model =list(ar = c(tmp.cor.lag.1)),
                                             n.start = 1, start.innov = c(diff.bm.logit/tmp.cor.lag.1),
                                             n = n.yrs.proj,innov = c(0,rnorm(n.yrs.proj-1,0,sd.bm.logit))) + mn.bm.logit))
       
        sim.Ks[[s]] <- data.frame(Years = 1:n.yrs.proj, sim = i,
                                       Stock = s, troph.cat = tl,
                                       prop.bm.stock = tmp.prop.bm,
                                       cor.prop.bm = NA,
                                       bm.stock = bm.trophic.Ks[[i]]$bm.tl[bm.trophic.Ks[[i]]$troph.cat==tl])
      
      
 
      
      # If there are only 2 stocks in a trophic level, then the second stock get the rest of the trophic levels biomass
     
    } # end the stocks loop
    # Now we need to get the proportions summing to 1
    #if(tl != 3)
    #{
      #browser()
      tl.stock.list <- as.data.frame(do.call('rbind',sim.Ks[tl.stocks]))
      tl.stock.list <- tl.stock.list |> dplyr::group_by(Years) |> dplyr::mutate(cor.prop.bm = prop.bm.stock/sum(prop.bm.stock))
      # Now remake the sim.Ks thing so it works with the below... clunky yes...
      for(s in tl.stocks) 
      {
        tl.stock.list.tmp <- tl.stock.list[tl.stock.list$Stock ==s,]
        sim.Ks[[s]] <- data.frame(Years = 1:n.yrs.proj, sim = i,
                                  Stock = s, troph.cat = tl,
                                  prop.bm.stock = tl.stock.list.tmp$prop.bm.stock,
                                  cor.prop.bm = tl.stock.list.tmp$cor.prop.bm,
                                  bm.stock = tl.stock.list.tmp$cor.prop.bm*bm.trophic.Ks[[i]]$bm.tl[bm.trophic.Ks[[i]]$troph.cat==tl])
      }
    #}# end the tl if
  } # end the trophic level loop

  sim.K.stock.lst[[i]] <- do.call("rbind",sim.Ks)
  
} # end the simulation loop

#browser()
sim.K.stocks.tmp <- do.call("rbind",sim.K.stock.lst)
sim.troph.K <- do.call("rbind",bm.trophic.Ks)
sim.eco.K <- do.call("rbind",sim.eco.bm)

# Get the meta data into the K stuff

sim.K.stocks <- left_join(sim.K.stocks.tmp,meta.dat,by=c("Stock",'troph.cat'))

#browser()
quants.K.stocks <- sim.K.stocks |> dplyr::group_by(Stock,common,species,Stock.short,Years,troph.cat) |> dplyr::summarize(mn = mean(bm.stock,na.rm=T),
                                                                                                                         log.mn = mean(log(bm.stock),na.rm=T),
                                                                                                                         med = median(bm.stock,na.rm=T),
                                                                                                                         sd = sd(log(bm.stock),na.rm=T))
quants.K.stocks$UCI <- exp(quants.K.stocks$log.mn + quants.K.stocks$sd)
quants.K.stocks$LCI <- exp(quants.K.stocks$log.mn - quants.K.stocks$sd)
# If this CI goes negative than set it to 1.  probably should do this on the log scale...
quants.K.stocks$LCI[quants.K.stocks$LCI < 0] <- 1

quants.K.tl <- sim.troph.K |> collapse::fgroup_by(troph.cat,Years) |> collapse::fsummarize(mn = mean(bm.tl,na.rm=T),
                                                                                           log.mn = mean(log(bm.tl),na.rm=T),
                                                                                           med = median(bm.tl,na.rm=T),
                                                                                           sd = sd(log(bm.tl),na.rm=T))
quants.K.tl$UCI <- exp(quants.K.tl$log.mn + quants.K.tl$sd)
quants.K.tl$LCI <- exp(quants.K.tl$log.mn - quants.K.tl$sd)

quants.K.eco <- sim.eco.K |> collapse::fgroup_by(Years) |> collapse::fsummarize(mn = mean(bm,na.rm=T),
                                                                                log.mn = mean(log(bm),na.rm=T),
                                                                                med = median(bm,na.rm=T),
                                                                                sd = sd(log(bm),na.rm=T))
quants.K.eco$UCI <- exp(quants.K.eco$log.mn + quants.K.eco$sd)
quants.K.eco$LCI <- exp(quants.K.eco$log.mn - quants.K.eco$sd)


# Wrap up the K time series for each simulation
#sim.K.stocks$Species <- substr(sim.K.stocks$Stock,14,100)
sim.stock.K.plt <- ggplot(sim.K.stocks) + geom_line(aes(x=Years+max(bm.best$Year)-1,y=bm.stock,group=sim),linewidth=2,alpha=0.2) + 
                                          geom_line(data=bm.best,aes(x=Year,y=bm.stock)) +
                                          facet_wrap(~troph.cat+Stock.short,scales='free_y') + scale_x_continuous(name='')+
                                          scale_y_continuous(name="Biomass Pool (tonnes)") + theme(legend.position = 'none') #+
                                          #guides(colour = guide_legend(nrow = 7))

save_plot(filename = paste0(repo.loc,"/Figures/BU/Simulation_stock_K.png"),sim.stock.K.plt,base_height = 10,base_width = 20)




sim.tl.K.plt <- ggplot(sim.troph.K) + geom_line(aes(x=Years+max(bm.best$Year)-1,y=bm.tl,group=as.factor(sim),color=as.factor(sim))) + 
                                      geom_line(data=bm.best,aes(x=Year,y=bm.tl)) +
                                      facet_wrap(~troph.cat) + theme(legend.position='none') + scale_x_continuous(name='')+ 
                                      scale_y_log10(name="Biomass Pool (tonnes)")
save_plot(filename = paste0(repo.loc,"/Figures/BU/Simulation_trophic_K.png"),sim.tl.K.plt,base_height = 10,base_width = 20)

sim.eco.K.plt <- ggplot(sim.eco.K) + geom_line(aes(x=Years+max(bm.best$Year)-1,y=bm,group=as.factor(sim),color=as.factor(sim))) +
                                      geom_line(data=bm.best,aes(x=Year,y=bm.eco)) +
                                    theme(legend.position = 'none')
save_plot(filename = paste0(repo.loc,"/Figures/BU/Simulation_eco_K.png"),sim.eco.K.plt,base_height = 8,base_width = 11)

# Now get the quantile plots with some uncertainty
#browser()
sim.stock.K.quant.plt <- ggplot(quants.K.stocks) + 
                          geom_line(aes(x=Years+max(bm.best$Year)-1,y=mn),linewidth=1,alpha=0.2) + 
                          geom_ribbon(aes(x=Years+max(bm.best$Year)-1,ymax=UCI,ymin=LCI),linewidth=1,alpha=0.2,fill='blue') + 
                          geom_line(data=bm.best,aes(x=Year,y=bm.stock)) +
                          facet_wrap(~Stock.short,scales='free_y') + scale_x_continuous(name='')+
                          scale_y_continuous(name="Biomass Pool (tonnes)") + theme(legend.position = 'none') #+
#guides(colour = guide_legend(nrow = 7))

save_plot(filename = paste0(repo.loc,"/Figures/BU/Simulation_quantile_stock_K.png"),sim.stock.K.quant.plt,base_height = 10,base_width = 20)

# Trophic level quantile plots

sim.TL.K.quant.plt <- ggplot(quants.K.tl) + 
  geom_line(aes(x=Years+max(trophic.bm.best$Year)-1,y=mn),linewidth=1.5,alpha=0.2) + 
  geom_ribbon(aes(x=Years+max(bm.best$Year)-1,ymax=UCI,ymin=LCI),linewidth=2,alpha=0.2,fill='blue') + 
  geom_line(data=trophic.bm.best,aes(x=Year,y=bm.tl)) +
  facet_wrap(~troph.cat) +
  scale_x_continuous(name='')+
  scale_y_continuous(name="Biomass Pool (tonnes)") + theme(legend.position = 'none') #+
#guides(colour = guide_legend(nrow = 7))

save_plot(filename = paste0(repo.loc,"/Figures/BU/Simulation_quantile_TL_K.png"),sim.TL.K.quant.plt,base_height = 10,base_width = 20)


# Ecosystem quantile plots
sim.eco.K.quant.plt <- ggplot(quants.K.eco) + 
  geom_line(aes(x=Years+max(eco.tot.bm.best$Year)-1,y=mn),linewidth=1.5,alpha=0.2) + 
  geom_ribbon(aes(x=Years+max(eco.tot.bm.best$Year)-1,ymax=UCI,ymin=LCI),linewidth=2,alpha=0.2,fill='blue') + 
  geom_line(data=eco.tot.bm.best,aes(x=Year,y=bm.eco)) +
  scale_x_continuous(name='')+
  scale_y_continuous(name="Biomass Pool (tonnes)") + theme(legend.position = 'none') #+
#guides(colour = guide_legend(nrow = 7))

save_plot(filename = paste0(repo.loc,"/Figures/BU/Simulation_quantile_eco_K.png"),sim.eco.K.quant.plt,base_height = 10,base_width = 20)



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

# Get the year range, going from the 'last' year to n.yrs.proj in the future, note this will go 1 year less than your intuition because
# we want n.yrs of data, i.e., 20 years is 2000 to 2019, not 2020... )
#browser()
years <- (last.year+1):(last.year+n.yrs.proj)
# Take the biomass data for the north sea and subset it to the years we have data
bm.mod.yrs <- bm.tst |> collapse::fsubset(Year %in% first.year:last.year)
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
av.wgt <- data.frame(Stock = av.wgt$Stock,troph.cat = av.wgt$troph.cat,mn.wgt = av.wgt$avg.weight)
# For some debugging, if still here you can delete I'm sure
#count = 0


#fake.tl.K <- 1e6



# So everything will need to get wrapped up in a simulation loop
res.ts <- NULL

ts.unpack <- NULL
results <- NULL
# define the trophic levels we have.

for(j in 1:n.sims)
{
  if(j == 1) count <- 0 
  Ks <- NULL
  st.time <- Sys.time()
  
  for(t in 1:n.yrs.proj)
  {
    # Now we need to run this through each trophic level, for bottom up we go 3, then 4, then 5....
    for(tl in troph.levels)
    {
      if(tl =="≤ 4.0") tl.K <- bm.sim.3[[j]][t]
      # This gets the Ks
      if(tl != "≤ 4.0") tl.K.space.prop <- Ks[[j]]$k.space[Ks[[j]]$tl==tl & Ks[[j]]$Years==t]
      
      bm.stash <- NULL
      tl.stocks <- unique(sim.K.stocks$Stock[sim.K.stocks$troph.cat ==tl])
      
      # Now get the stock biomass from last year.
      if(t ==1)  tl.bm.last <- init.tl.bm[init.tl.bm$troph.cat==tl,]
      
      # Then we'll need to get these from the model simulations.
      #if(t==2) browser()
      if(t > 1)
      {
        # Use the handy av.wgt data.frame I made above
        bm.stocks <- data.frame(bm = NA,meta.dat)
        for(s in stock.eco) bm.stocks$bm[bm.stocks$Stock == s] <-  res.ts[[t-1]]$bm[res.ts[[t-1]]$Stock ==s]
        bm.stocks$bm <- bm.stocks$bm
        stock.bm.last <- bm.stocks
        eco.bm.last <- sum(bm.stocks$bm)
        tl.bm.last <- bm.stocks |> collapse::fsubset(troph.cat ==tl) |> collapse::fsummarise(bm.tl = sum(bm))
      }  
      # Now we need to figure out what K space is available for each stock within the trophic level.
      # I'm going to base the K space for t+1 on the biomass available in the higher trophic level
      # in year t along with the 'transfer efficiency', which is currently the proportion of the total 
      # ecosystem biomass that the trophic level gets.
      # We can then adjust the stock K's by the available K space in each stock
      # So this is the K space available in a given trophic level in a year
      if(tl =="≤ 4.0") 
      {
        #browser()
        base.eco.K.tmp <- sim.eco.K |> collapse::fsubset(sim == j & Years ==t)
        # Get the trophic level 3 only since we are going bottom up
        base.tl.K.tmp <- sim.troph.K |> collapse::fsubset(sim == j & Years ==t & troph.cat==tl)
        base.stock.K.tmp <- sim.K.stocks |> collapse::fsubset(sim == j & Years ==t & troph.cat==tl)
        #base.stock.K.tmp$prop.bm.stock <- base.stock.K.tmp
        base.stock.K.tmp$tl.K <- tl.K 
        base.tl.K.tmp$prop.K.space <- base.tl.K.tmp$bm.tl/tl.bm.last$bm.tl
        # We can then adjust the stock K's by the available K space in each stock
        base.stock.K.tmp$K.space <- NA
        base.stock.K.tmp$K.space <- base.stock.K.tmp$bm.stock[base.stock.K.tmp$troph.cat ==tl] * 
          (base.tl.K.tmp$prop.K.space-1)
        base.stock.K.tmp$adj.K <- base.stock.K.tmp$bm.stock + base.stock.K.tmp$K.space
      } # end if tl =="≤ 4.0"
      
      if(tl !="≤ 4.0") 
      {
        #browser()
        base.stock.K.tmp <- sim.K.stocks |> collapse::fsubset(sim == j & Years ==t & troph.cat == tl)
        base.tl.K.tmp <- sim.troph.K |> collapse::fsubset(sim == j & Years ==t & troph.cat == tl)
        
        # Get the new trophic level values right
        base.tl.K.tmp$prop.K.space <- tl.K.space.prop
        base.tl.K.tmp$adj.K <- base.tl.K.tmp$prop.K.space * base.tl.K.tmp$bm.tl
        base.tl.K.tmp$K.space <- base.tl.K.tmp$adj.K - base.tl.K.tmp$bm.tl
        
        # Now get the stock right, all we have to do is multiply the bm.stock
        # by the prop.K.space (i.e. what the proportinonal chance in the avilable K-space is)
        base.stock.K.tmp$adj.K <- base.tl.K.tmp$prop.K.space * base.stock.K.tmp$bm.stock
        base.stock.K.tmp$K.space <- base.stock.K.tmp$adj.K - base.stock.K.tmp$bm.stock
        
      }
      # So now I have Carrying Capacities that take up (or lose) any available K space.
      # Now we can convert these to numbers using the historic 'average weight' of the stocks, to avoid complication
      # I'm just using the average of the average weight for each stock...
      #browser()
      base.stock.K.tmp <- left_join(base.stock.K.tmp,av.wgt,by=c("Stock","troph.cat"))
      # And now we can get a K in numbers....
      base.stock.K.tmp$adj.K.num <- base.stock.K.tmp$adj.K/base.stock.K.tmp$mn.wgt
      
      for(s in tl.stocks)
      {
        count <- count + 1
        stock.lambdas <- lambdas[[s]] 
        bm.ts.stock <- bm.final[bm.final$Stock == s & bm.final$Year %in% first.year:last.year,]  
      
      if(t == 1) 
      { 
        bm.start <- bm.ts.stock$bm.stock[bm.ts.stock$Year == last.year]
        results[[s]] <- data.frame(bm = bm.start,removals = NA,ex.rate = NA,
                                  Stock = s,sim= j,lambda = NA,Years=t-1,troph.cat = tl,
                                  K.bm = NA)
        
      } else{ bm.start <- res.ts[[t-1]]$bm[res.ts[[t-1]]$Stock ==s]}
      
      
     
      # Since I have Years and sim recorded, I should just be able to recursivly rbind this...
      
      # This should do the trick for creating this for the very first stop and then updating thereafter.
      # if(t ==1 & j ==1 & count ==1) 
      # {
      #   base.stock.K <- base.stock.K.tmp
      #   base.tl.K <- base.tl.K.tmp
      #   base.eco.K <- base.eco.K.tmp
      # } else {
      #   base.stock.K <- rbind(base.stock.K,base.stock.K.tmp)
      #   base.tl.K <- rbind(base.tl.K,base.tl.K.tmp)
      #   base.eco.K <- rbind(base.eco.K,base.eco.K.tmp)
        
      #} # end the else...
      #browser()
      l.v.h <- 0.4
      #if(s == "ICES-HAWG_ NS-IV 3a,7d_Clupea_harengus")  l.v.h <- 0.6 # DK Note, using 0.6 for herring stock didn't decline below 50% in this time period.
      #if(s == "ICES-WGNSSK_NS4 _Scopthalmus_maximus")  l.v.h <- 0.9 # # DK Note, using 0.9 for this stock because it only declined to 66% of max in time period.
      #if(t==25) 
      cur.K <- base.stock.K.tmp$adj.K[base.stock.K.tmp$Stock ==s]
      init.K <- base.stock.K.tmp$bm.stock[base.stock.K.tmp$Stock ==s]
      
      # Here is where we insert the management to set the catch.
      # Note this is using all the data, not just the data in recent years to calculate the RPs and RR...
      #If you don't do anything, the default is to fish it using the 
      # expoitation from the time series as the RR
      # the LRP is 40% of median
      # and the USR is the median
      # No fishing below the LRP, and fishing turns down linearily between the USR and LRP.
      if(is.null(exploit))
      {
        bm.n.er.hist <- bm.best[bm.best$Stock == s,]
        ex.dat <- data.frame(lrp = NA, urp = NA, rr= NA,er.mn = NA,er.below.lrp = NA)
        # If no data, we'll make the lrp be 40% of the median historic biomass
        ex.dat$lrp <- 0.4*median(bm.n.er.hist$bm.stock,na.rm=T)
        # If no data, we'll make the urp be 100% of the median historic biomass
        ex.dat$urp <- median(bm.n.er.hist$bm.stock,na.rm=T)
          # If no data, we'll make the rr be the median historic exploitation rate
        ex.dat$rr <- median(exploit[[s]]$er,na.rm=T)
        # If we didn't set the rr, we'll make it be the rr
        ex.dat$er.mn <- exploit[[s]]$rr
        # If we didn't set the exploitation rate below the lrp, we'll make it 0
        ex.dat$er.below.lrp <- 0
        ex.dat$er.sd <- 0.1
      } #end if exploit is null.
      #if(s == "ICES-HAWG_NS_Ammodytes_tobianus") browser()
      #ex.dat <- exploit[exploit$stock == s,]
      ## Now calculate the exploitation rate, which will be based on lrp and usr.
     # if (s == "ICES-WGNSSK_NS4_Solea_solea") browser()

      er <- proj.catch.eqn(dat = ex.dat,bm = bm.start)
      
      res <- pop.dam(lambda = stock.lambdas,stock.K = cur.K,bm.start = bm.start, catch = er, stock = s,
                     bm.stock = bm.ts.stock,low.vs.high = l.v.h,method="not_sample")
      #if(t==2) browser()
      results[[s]] <- data.frame(bm = res$tst.res,removals =res$removals,ex.rate = res$ex.rate,
                                 Stock = s,sim= j,lambda = res$lambda,Years=t,
                                 troph.cat = tl,K.bm = res$K)
      #browser()
      bm.stash[[s]] <- data.frame(k.stock = res$K,bm.stock = res$tst.res,k.space = res$K - res$tst.res,K.init =init.K)
      } # end stock loop by trophic level
      # How this is set up we have the biomass after fishing impacting the K in that year
      # We might want to have this impact the following years K.
      # K for the next trophic level
      if(tl =="≤ 4.0") 
      {
        #browser()
        bms <- do.call('rbind', bm.stash)
        # This way takes the proportion of k-space available (or missing) and removes it from the trophic levels below
        # i.e. if the lower trophic levels isn't a K, then the higher trophic levels are missing food
        # so 90% of their initially allocated K (I have been debating using the K-space, but I think is problematic 
        # in ways that don't make biological sense to propogate up the chain, I think
        # as we move up the chain means everyone below is hurting too.
        # Doing it this way to ensure that the K's remain positive (by adding actual space you can get negatives....)
        next.tl.K.space.prop <- (sum(bms$bm.stock)/sum(bms$K.init))#/mean(sim.tl.4.to.5[[j]]) # Old way this was used to give them the tl efficicnecy, worked badly...
        #next.tl.K <- sum(do.call('rbind',bm.stash))*sim.tl.3.to.4[[j]][t]
      }
      if(tl =="4.1-4.9") 
      {
        #if(t%in% c(1,10,20,50)) browser()
        bms <- do.call('rbind', bm.stash)
        # This way takes the proportion of k-space available (or missing) and removes it from the trophic levels below
        # i.e. if the lower trophic levels isn't a K, then the higher trophic levels are missing food
        # so 90% of their initially allocated K (I have been debating using the reallocated K-space, but I think is problematic as the collapse of the lower levels makes the K space collapse
        # in ways that don't make biological sense to propogate up the chain, I think
        # as we move up the chain means everyone below is hurting too.
        # Doing it this way to ensure that the K's remain positive (by adding actual space you can get negatives....)
        next.tl.K.space.prop <- (sum(bms$bm.stock)/sum(bms$K.init))#/mean(sim.tl.4.to.5[[j]]) # Old way this was used to give them the tl efficicnecy, worked badly...
        #next.tl.K <- sum(do.call('rbind',bm.stash))*sim.tl.3.to.4[[j]][t]
        
      }
      
      if(tl != "≥ 5.0") 
      {
        mv.to.next <- which(troph.levels == tl)+1
        Ks[[j]] <- rbind(Ks[[j]],data.frame(k.space.prop=next.tl.K.space.prop,tl = troph.levels[mv.to.next],Years=t,sim=j))
      }
    }# End trophic level loop.
    
  #browser()
  res.ts[[t]] <- do.call("rbind",results)
  } # end the t looping through each year.
  
#ggplot(base.stock.K) + geom_line(aes(x= Years,y=bm.stock,group=Stock,color=Stock)) + facet_wrap(~troph.cat) + scale_y_log10()
  #browser()
  # Unpack the results
  #st.time <- Sys.time()
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
saveRDS(object = ts.final,file = paste0(repo.loc,"/Results/BU/NS_projections_",n.sims,"_sims_",min(years),"_to_",max(years),
                                         "_time_series_projections.Rds"))

saveRDS(object = quants,file = paste0(repo.loc,"/Results/BU/NS_projections_",n.sims,"_sims_",min(years),"_to_",max(years),
                                       "_time_series_quantiles.Rds"))
# 
# saveRDS(object = r.final,file = paste0(repo.loc,"/Results/BU/NS_projections_",n.sims,"_sims_",min(years),"_to_",max(years),
#                                        "_r_projections.Rds"))


n.breaks <- 6
# Two simple plots. 
p.sims <- ggplot(ts.final ) + geom_line(aes(x=Years,y=bm,group = sim,color=sim),alpha=0.8) +
  facet_wrap(~Stock) + 
  scale_x_continuous(name='',breaks = seq(0,n.yrs.proj,length=n.breaks)) +
  scale_y_log10(name = "Biomass Pool (tonnes)") + 
  theme(legend.position = 'none') 

save_plot(paste0(repo.loc,"/Figures/BU/biomass_trends.png"),p.sims,base_height = 12,base_width = 24)


colors <- distinct(bm.best, spec.tl, color)
pal <- colors$color
names(pal) <- colors$spec.tl

# p.sims.quants <- ggplot(quants) + geom_line(aes(x=Years,y=med,group=Stock,color=spec.tl)) + 
#   facet_wrap(~troph.cat,scales = 'free_y') +  scale_y_log10(name="Biomass Pool (tonnes)") +   theme(legend.position = 'top') +
#   guides(colour = guide_legend(nrow = 7)) +
#   scale_x_continuous(breaks = seq(1,50,by=49),labels=c(2015,2065)) 
#   #geom_ribbon(data=quants, aes(x=Years,ymax=U.50,ymin = L.50),alpha=0.5,fill='blue',color='blue') 
# save_plot(paste0(repo.loc,"/Figures/BU/Quantile_abundance_trends.png"),p.sims.quants,base_height = 8,base_width = 16)



#

p.sims.quants <- ggplot(quants) + geom_line(aes(x=Years,y=med,group=Stock,color=spec.tl),linewidth=2) + 
  facet_wrap(~troph.cat) +  scale_y_log10(name="Biomass Pool (tonnes)") +   theme(legend.position = 'top') +
  guides(colour = guide_legend(nrow = 5)) + scale_color_manual(values=pal) +
  scale_x_continuous(name="",breaks = ) 
#geom_ribbon(data=quants, aes(x=Years,ymax=U.50,ymin = L.50),alpha=0.5,fill='blue',color='blue') 
save_plot(paste0(repo.loc,"/Figures/BU/Quantile_biomass_trends.png"),p.sims.quants,base_height = 8,base_width = 16)



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