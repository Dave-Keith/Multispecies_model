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




trophic.mod<-function(dat=bm.best,n.yrs.proj = 50, n.sims = 20,
                      manage = NULL, method = 'log_linear', mod.pred = NULL,
                      repo.loc = "D:/GitHub/Multispecies_model")
{
  
set.seed(1)


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
          "https://raw.githubusercontent.com/dave-keith/ICM/main/Scripts/functions/forward_project.r",
          "https://raw.githubusercontent.com/dave-keith/Multispecies_model/main/Scripts/NS_caatch_function.R"
)
# Now run through a quick loop to load each one, just be sure that your working directory is read/write!
for(fun in funs) 
{
  download.file(fun,destfile = basename(fun))
  source(paste0(getwd(),"/",basename(fun)))
  file.remove(paste0(getwd(),"/",basename(fun)))
}


########################## Section 2 Parameters ########################## Section 2 Parameters ########################## Section 2 Parameters
stock.eco <- unique(bm.best$Stock)
bm.best <- dat

# The 'ecosystem' biomass and numbers
com.tot.bm <- bm.best |> collapse::fgroup_by(Year) |> 
  collapse::fsummarize(num.eco = sum(num.stock),bm.eco = sum(bm.stock))
# Trophic level biomass and numbers.
trophic.bm <- bm.best |> collapse::fgroup_by(Year,troph.cat) |> 
  collapse::fsummarize(num.tl = sum(num.stock),bm.tl = sum(bm.stock))
#trophic.bm$troph.cat <- factor(trophic.bm$troph.cat,levels = c("Low","Medium","High"))

# All the bm together
tl.com.bm <- left_join(trophic.bm,com.tot.bm,by="Year")
tl.com.bm$prop.bm.tl <- tl.com.bm$bm.tl/tl.com.bm$bm.eco
tl.com.bm$prop.num.tl <- tl.com.bm$num.tl/tl.com.bm$num.eco
# So here we are working to get the 'ecosystem' carrying capacity by looking at the total biomass for the NS stocks we have
# data for over the period of time we have data for all the stocks.
# So here we pull out the data we need to look at total abundance and total biomass in the system by year...

# The 'transfer efficiency' between our trophic levels
tl.3.to.4 <- bm.best$prop.bm.tl[bm.best$troph.cat=="Medium"][1:n.years]/bm.best$prop.bm.tl[bm.best$troph.cat=="Low"][1:n.years]
tl.4.to.5 <- bm.best$prop.bm.tl[bm.best$troph.cat=="High"][1:n.years]/bm.best$prop.bm.tl[bm.best$troph.cat=="Medium"][1:n.years]
tl.3.to.5 <- bm.best$prop.bm.tl[bm.best$troph.cat=="High"][1:n.years]/bm.best$prop.bm.tl[bm.best$troph.cat=="Low"][1:n.years]


# So now we want to look at stock level within a trophic level
# add some colors... this is so clunky....
tl3s <- unique(bm.best$Species[bm.best$troph.cat=="Low"])
tl4s <- unique(bm.best$Species[bm.best$troph.cat=="Medium"])
tl5s <- unique(bm.best$Species[bm.best$troph.cat=="High"])
bm.best$color <- "black"
count=1
for(c in tl3s) 
{
  if(count == 2) bm.best$color[bm.best$Species == tl3s[count]] <- "blue"
  if(count == 3) bm.best$color[bm.best$Species == tl3s[count]] <- "green"
  if(count == 4) bm.best$color[bm.best$Species == tl3s[count]] <- "grey"
  if(count == 5) bm.best$color[bm.best$Species == tl3s[count]] <- "orange"
  if(count == 6) bm.best$color[bm.best$Species == tl3s[count]] <- "firebrick2"
  count = count + 1
}
# TL4 colors
count=1
for(c in tl4s) 
{
  if(count == 2) bm.best$color[bm.best$Species == tl4s[count]] <- "blue"
  if(count == 3) bm.best$color[bm.best$Species == tl4s[count]] <- "green"
  if(count == 4) bm.best$color[bm.best$Species == tl4s[count]] <- "grey"
  if(count == 5) bm.best$color[bm.best$Species == tl4s[count]] <- "orange"
  if(count == 6) bm.best$color[bm.best$Species == tl4s[count]] <- "firebrick2"
  count = count + 1
}
# TL5 colors
count=1
for(c in tl5s) 
{
  if(count == 2) bm.best$color[bm.best$Species == tl5s[count]] <- "blue"
  if(count == 3) bm.best$color[bm.best$Species == tl5s[count]] <- "green"
  if(count == 4) bm.best$color[bm.best$Species == tl5s[count]] <- "grey"
  if(count == 5) bm.best$color[bm.best$Species == tl5s[count]] <- "orange"
  if(count == 6) bm.best$color[bm.best$Species == tl5s[count]] <- "firebrick2"
  count = count + 1
}

# Put in Species + trophic level
bm.best$spec.tl <- paste(bm.best$Species,"(TL is ",bm.best$trophic,")")
# Pull out meta data
meta.dat <- bm.best |> dplyr::group_by(Stock,trophic,Species,troph.cat,color,spec.tl) |> filter(row_number() >= (n() ))
meta.dat <- meta.dat[,c("Stock","trophic","Species","troph.cat","color","spec.tl","Stock.short",'common')]
#meta.dat$troph.cat <- as.numeric(meta.dat$troph.cat)

colors <- distinct(bm.best, spec.tl, color)
pal <- colors$color
names(pal) <- colors$spec.tl

stock.prop.bm.plt <-  ggplot(bm.best) + geom_line(aes(x=Year,y=prop.bm.stock.tl,group = Stock,color=spec.tl),linewidth=2) + 
                                        facet_wrap(~troph.cat) + guides(colour = guide_legend(nrow = 5)) + theme(legend.position = 'top',legend.title = element_blank()) +
                                        scale_y_log10(name= "Proportion of biomass Pool",n.breaks=10) + scale_x_continuous(name="",breaks=seq(1970,2200,by=5))+
                                        scale_color_manual(values=pal)
save_plot(paste0(repo.loc,"/Figures/BU/Historic_Prop_Biomass_ns_by_stock.png"),stock.prop.bm.plt,base_height = 8,base_width = 15)

stock.bm.plt <-ggplot(bm.best) + geom_line(aes(x=Year,y=bm.stock/1e3,group = Stock,color=spec.tl),linewidth=2) + 
                                  facet_wrap(~troph.cat) + scale_x_continuous(name="",breaks=seq(1970,2200,by=5))+
                                  scale_y_log10(name = "Biomass Pool (thousands of tonnes)",n.breaks=7) + theme(legend.position = 'top',legend.title = element_blank()) +
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
 

com.tot.bm.best <- com.tot.bm |> collapse::fsubset(Year %in% first.year:last.year)
trophic.bm.best <- trophic.bm |> collapse::fsubset(Year %in% first.year:last.year)

# The correlation in the ecosystem biomass trend, can see this is an AR1
K.cor <- pacf(com.tot.bm.best$bm.eco,plot=F)
# The cross correlation between the ecosystem biomass trend and the trophic level biomasses
# All correlated, but strongest is unsurprisingly the link between the the ecosystem and the biomass in the
# lowest TL. I suspect this may structurally come out even without explicity building in a lot of
# correlation structure to the models.
K.tl.3.cor <- ccf(com.tot.bm.best$bm.eco,trophic.bm.best$bm.tl[trophic.bm.best$troph.cat=="Low"],plot = F)
K.tl.4.cor <- ccf(com.tot.bm.best$bm.eco,trophic.bm.best$bm.tl[trophic.bm.best$troph.cat=="Medium"],plot = F)
K.tl.5.cor <- ccf(com.tot.bm.best$bm.eco,trophic.bm.best$bm.tl[trophic.bm.best$troph.cat=="High"],plot = F)
# Within trophic levels...
# So these 3 mostly say if the biomass is up one TL, it is up in all TLs, tho there might be some negative between 3 and 4
# at Lag -1 (though that's not quite significant)
tl.3.4.cor <- ccf(trophic.bm.best$bm.tl[trophic.bm.best$troph.cat =="Low"],trophic.bm.best$bm.tl[trophic.bm.best$troph.cat =="Medium"],plot = F)
tl.3.5.cor <- ccf(trophic.bm.best$bm.tl[trophic.bm.best$troph.cat =="Low"],trophic.bm.best$bm.tl[trophic.bm.best$troph.cat =="High"],plot = F)
tl.4.5.cor <- ccf(trophic.bm.best$bm.tl[trophic.bm.best$troph.cat =="Medium"],trophic.bm.best$bm.tl[trophic.bm.best$troph.cat =="High"],plot = F)
# Looking at proportions, need to stew a bit on this because there is necessarily some 
# correlation built into proportions, but what is interesting is that
# the correlation strength is really really high between TL 3 and TL 4, it is weaker (more diffuse really) at TL 5
# and there is no correlation between 4 and 5
tl.3.4.prop.cor <- ccf(bm.best$prop.bm.tl[bm.best$troph.cat =="Low"][1:n.years],bm.best$prop.bm.tl[bm.best$troph.cat =="Medium"][1:n.years],plot = F)
tl.3.5.prop.cor <- ccf(bm.best$prop.bm.tl[bm.best$troph.cat =="Low"][1:n.years],bm.best$prop.bm.tl[bm.best$troph.cat =="High"][1:n.years],plot = F)
tl.4.5.prop.cor <- ccf(bm.best$prop.bm.tl[bm.best$troph.cat =="Medium"][1:n.years],bm.best$prop.bm.tl[bm.best$troph.cat =="High"][1:n.years],plot = F)

# So now really what I need to do first is make a quick simulation that gets me ecosystem K, trophic level K, and stock K
# once I have those then we just run the models :-)
# I am saying the actual ecosystem biomass is the maximum observered instead of the mean, this will give some
# spaces between the biomass and the Pool sizes that we didn't get when using the mean
max.com.bm <- max(com.tot.bm.best$bm.eco,na.rm=T) 
#start.com.sim <- com.tot.bm.best$bm.eco[length(com.tot.bm.best$bm.eco)]
sd.com.bm <- sd(com.tot.bm.best$bm.eco)
# trophic level biomass and proportions... for the proportion will probably wanna sample from a beta distro
# So not sure how to do that nicely...
# 


# First get the ecosystem biomass in a correlated time series, there are a whole lot of ways one could do this, this
# is one of many different ideas. I think we could get the 4 and 5 correlations better another way, but
# For a first pass I'm ok with this.
# Ok, duh, use the mean of the time series then the arima gives us the deviations from that mean and we get a nice time series.

# Used for simulations to get good time series for the K for TL3,4, and 5 
tl.3.prop.bm.ts <- bm.best$prop.bm.tl[bm.best$troph.cat=="Low"][1:n.years]
# Extract the frist two components from the pacf to get the two AR components from the model.
tl.3.prop.pacf <- pacf(tl.3.prop.bm.ts,plot = F)
tl.3.prop.bm.lag.1 <- tl.3.prop.pacf$acf[1]
tl.3.prop.bm.lag.2 <- tl.3.prop.pacf$acf[2]
# TL 4 and 5 splits historically
# This is the correlation between 4 and 5# This is the correlation between 4 and 5
tl.4.5.prop.bm <- bm.best$bm.tl[bm.best$troph.cat =="High"][1:n.years]/(bm.best$bm.tl[bm.best$troph.cat =="Medium"][1:n.years]+
                                                                          bm.best$bm.tl[bm.best$troph.cat =="High"][1:n.years])
tl.4.5.prop.4.5.bm <- pacf(tl.4.5.prop.bm,plot = F)


troph.levels <- sort(unique(bm.best$troph.cat))

# 
sim.pool.stock.lst <- NULL
sim.pools <- NULL
sim.com.bm <- NULL
bm.trophic.pools <- NULL
bm.sim.5 <- NULL
bm.sim.4 <- NULL
bm.sim.3 <- NULL
sim.tl.3.to.4 <- NULL
sim.tl.4.to.5 <- NULL

# Get necessary data on logit scale
tl.3.logit <- logit(tl.3.prop.bm.ts)
tl.4.5.logit <- logit(tl.4.5.prop.bm)

# Starting values for the ecosystem and the proportions, logit needed for arima models with the proportions
# For the community I am starting it at the maximum and running from there
start.com.sim <- max.com.bm#com.tot.bm.best$bm.eco[nrow(com.tot.bm.best)]
# Tor the rest I am starting it at where they ended.
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
start.com.diff = start.com.sim- max.com.bm
start.tl.3.diff <- start.tl.3.logit - mn.tl.3.logit
start.tl.4.5.diff <- start.tl.4.5.logit - mn.tl.4.5.logit

# the standard deviations
sd.tl.3.logit <- sd(tl.3.logit)
sd.tl.4.5.logit <- sd(tl.4.5.logit)
# Lag for the Arima model
tl.4.5.prop.bm.lag.1 <- tl.4.5.prop.4.5.bm$acf[1]
# convert to logit scale for the arima models

# NEED TO SORT OUT THE INDEXING ON THE Bm.sim.tl objects and how I use these later on!!
for(i in 1:n.sims) 
{
 # The ecosystem K, using the mean of the ecosystem with the correlation observed of the time series.
 # This starts the time series at the last value of the time series, then moves it to the mean value, bam!!  This will be done for each of these arima sims.
  sim.com.bm[[i]] <- data.frame(bm = c(arima.sim(model =list(ar = K.cor$acf[1]),n = n.yrs.proj,n.start=1,start.innov = start.com.diff/K.cor$acf[1],
                                                 innov = c(0,rnorm(n.yrs.proj-1,0,sd.com.bm))) + max.com.bm),
                                                 Years = 1:n.yrs.proj,sim = i) 
  #pacf(sim.com.bm[[i]]$bm) # looks good

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
  
  bm.sim.3[[i]] <- sim.tl.3.prop.bm * sim.com.bm[[i]]$bm
  # So this is what is left for TL 4 and 5
  bm.left.4.5<- sim.com.bm[[i]]$bm - bm.sim.3[[i]]
  # So then we use the historical split between 4 and 5 can see 5 gets about 1/3-1-5 of 3
   # so then simulate this split
  sim.tl.5.4.prop.bm <- inv.logit(mn.tl.4.5.logit + 
                                    arima.sim(model =list(ar = tl.4.5.prop.bm.lag.1),n = n.yrs.proj,
                                              n.start =1, start.innov = c(start.tl.4.5.diff/tl.4.5.prop.bm.lag.1), 
                                              innov = c(0,rnorm(n.yrs.proj-1,0,sd.tl.3.logit))))
  # And now TL 5 gets this proportion of the 4 and 5 biomass
  bm.sim.5[[i]] <- bm.left.4.5 * sim.tl.5.4.prop.bm
  # And TL4 gets the rest, and so the ecosystem biomass is a portion of the whole biomass
  bm.sim.4[[i]] <- sim.com.bm[[i]]$bm - bm.sim.3[[i]]-bm.sim.5[[i]]
  
  bm.trophic.pools[[i]] <- data.frame(Years = rep(1:n.yrs.proj,3), sim =i,
                                   bm.tl = c(bm.sim.3[[i]],bm.sim.4[[i]],bm.sim.5[[i]]),troph.cat = sort(rep(troph.levels,n.yrs.proj)),
                                   bm.eco = rep(sim.com.bm[[i]]$bm,3))
  bm.trophic.pools[[i]]$prop.bm.tl <- bm.trophic.pools[[i]]$bm.tl/bm.trophic.pools[[i]]$bm.eco
  
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
        tmp.cor <- pacf(tmp.dat$prop.bm.tl,plot=F) # Get the correlation, use AR1 and no more.
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
        #browser()
        tmp.prop.bm <- c(inv.logit(arima.sim(model =list(ar = c(tmp.cor.lag.1)),
                                             n.start = 1, start.innov = c(diff.bm.logit/tmp.cor.lag.1),
                                             n = n.yrs.proj+1,innov = c(0,rnorm(n.yrs.proj,0,sd.bm.logit))) + mn.bm.logit))[2:(n.yrs.proj+1)]
       
        sim.pools[[s]] <- data.frame(Years = 1:n.yrs.proj, sim = i,
                                       Stock = s, troph.cat = tl,
                                       prop.bm.stock = tmp.prop.bm,
                                       cor.prop.bm = NA,
                                       bm.stock = bm.trophic.pools[[i]]$bm.tl[bm.trophic.pools[[i]]$troph.cat==tl])
      
      
 
      
      # If there are only 2 stocks in a trophic level, then the second stock get the rest of the trophic levels biomass
     
    } # end the stocks loop
    # Now we need to get the proportions summing to 1
    #if(tl != 3)
    #{
      #
      tl.stock.list <- as.data.frame(do.call('rbind',sim.pools[tl.stocks]))
      tl.stock.list <- tl.stock.list |> dplyr::group_by(Years) |> dplyr::mutate(cor.prop.bm = prop.bm.stock/sum(prop.bm.stock))
      # Now remake the sim.pools thing so it works with the below... clunky yes...
      for(s in tl.stocks) 
      {
        tl.stock.list.tmp <- tl.stock.list[tl.stock.list$Stock ==s,]
        sim.pools[[s]] <- data.frame(Years = 1:n.yrs.proj, sim = i,
                                  Stock = s, troph.cat = tl,
                                  prop.bm.stock = tl.stock.list.tmp$prop.bm.stock,
                                  cor.prop.bm = tl.stock.list.tmp$cor.prop.bm,
                                  bm.stock = tl.stock.list.tmp$cor.prop.bm*bm.trophic.pools[[i]]$bm.tl[bm.trophic.pools[[i]]$troph.cat==tl])
      }
    #}# end the tl if
  } # end the trophic level loop

  sim.pool.stock.lst[[i]] <- do.call("rbind",sim.pools)
  
} # end the simulation loop

#
sim.pool.stocks.tmp <- do.call("rbind",sim.pool.stock.lst)
sim.troph.pool <- do.call("rbind",bm.trophic.pools)
sim.com.pool <- do.call("rbind",sim.com.bm)

# Get the meta data into the K stuff

sim.pool.stocks <- left_join(sim.pool.stocks.tmp,meta.dat,by=c("Stock",'troph.cat'))

#
quants.pool.stocks <- sim.pool.stocks |> dplyr::group_by(Stock,common,Species,Stock.short,Years,troph.cat) |> dplyr::summarize(mn = mean(bm.stock,na.rm=T),
                                                                                                                         log.mn = mean(log(bm.stock),na.rm=T),
                                                                                                                         med = median(bm.stock,na.rm=T),
                                                                                                                         sd = sd(log(bm.stock),na.rm=T))
quants.pool.stocks$UCI <- exp(quants.pool.stocks$log.mn + quants.pool.stocks$sd)
quants.pool.stocks$LCI <- exp(quants.pool.stocks$log.mn - quants.pool.stocks$sd)
# If this CI goes negative than set it to 1.  probably should do this on the log scale...
quants.pool.stocks$LCI[quants.pool.stocks$LCI < 0] <- 1

quants.pool.tl <- sim.troph.pool |> collapse::fgroup_by(troph.cat,Years) |> collapse::fsummarize(mn = mean(bm.tl,na.rm=T),
                                                                                           log.mn = mean(log(bm.tl),na.rm=T),
                                                                                           med = median(bm.tl,na.rm=T),
                                                                                           sd = sd(log(bm.tl),na.rm=T))
quants.pool.tl$UCI <- exp(quants.pool.tl$log.mn + quants.pool.tl$sd)
quants.pool.tl$LCI <- exp(quants.pool.tl$log.mn - quants.pool.tl$sd)

quants.pool.eco <- sim.com.pool |> collapse::fgroup_by(Years) |> collapse::fsummarize(mn = mean(bm,na.rm=T),
                                                                                log.mn = mean(log(bm),na.rm=T),
                                                                                med = median(bm,na.rm=T),
                                                                                sd = sd(log(bm),na.rm=T))
quants.pool.eco$UCI <- exp(quants.pool.eco$log.mn + quants.pool.eco$sd)
quants.pool.eco$LCI <- exp(quants.pool.eco$log.mn - quants.pool.eco$sd)

# Number of beaks for the x-axis of the figures

n.breaks <- 6
start.year <- 10*floor(min(bm.best$Year/10))
end.year <- n.yrs.proj+10*ceiling(max(bm.best$Year/10))

few.breaks <- floor(seq(start.year,end.year,length=n.breaks))
lotsa.breaks <- floor(seq(start.year,end.year,length=2*n.breaks))
# Wrap up the K time series for each simulation
#sim.pool.stocks$Species <- substr(sim.pool.stocks$Stock,14,100)
sim.stock.pool.plt <- ggplot(sim.pool.stocks) + geom_line(aes(x=Years+max(bm.best$Year)-1,y=bm.stock/1e3,group=sim),linewidth=2,alpha=0.2,color='grey50') + 
                                            geom_line(data=bm.best,aes(x=Year,y=bm.stock/1e3),color='black') +
                                            facet_wrap(~troph.cat+Stock.short,scales='free_y') +  scale_x_continuous(name='',breaks = few.breaks)+ 
                                            scale_y_continuous(name="Biomass Pool (thousands of tonnes)",limits=c(0,NA)) + theme(legend.position = 'none') #+

save_plot(filename = paste0(repo.loc,"/Figures/BU/Simulation_stock_K.png"),sim.stock.pool.plt,base_height = 10,base_width = 20)


sim.tl.pool.plt <- ggplot(sim.troph.pool) + geom_line(aes(x=Years+max(bm.best$Year)-1,y=bm.tl/1e6,group=as.factor(sim),color=as.factor(sim))) + 
                                      geom_line(data=bm.best,aes(x=Year,y=bm.tl/1e6)) +
                                      facet_wrap(~troph.cat,scales='free_y') + theme(legend.position='none') + scale_x_continuous(name='',breaks = few.breaks) + 
                                      scale_y_continuous(name="Biomass Pool (millions of tonnes)",limits=c(0,NA))  
save_plot(filename = paste0(repo.loc,"/Figures/BU/Simulation_trophic_K.png"),sim.tl.pool.plt,base_height = 10,base_width = 20)

sim.com.pool.plt <- ggplot(sim.com.pool) + geom_line(aes(x=Years+max(bm.best$Year)-1,y=bm/1e6,group=as.factor(sim),color=as.factor(sim))) +
                                      geom_line(data=bm.best,aes(x=Year,y=bm.eco/1e6)) + scale_x_continuous(name='',breaks = lotsa.breaks)+ 
                                      scale_y_continuous(name="Biomass Pool (millions of tonnes)",limits=c(0,NA)) + theme(legend.position = 'none')
save_plot(filename = paste0(repo.loc,"/Figures/BU/Simulation_eco_K.png"),sim.com.pool.plt,base_height = 8,base_width = 11)

# Now get the quantile plots with some uncertainty
#
sim.stock.pool.quant.plt <-  ggplot(quants.pool.stocks) + geom_line(aes(x=Years+max(bm.best$Year)-1,y=mn/1e3),linewidth=1,alpha=0.2) + 
                                                    geom_ribbon(aes(x=Years+max(bm.best$Year)-1,ymax=UCI/1e3,ymin=LCI/1e3),linewidth=1,alpha=0.2,fill='blue') + 
                                                    geom_line(data=bm.best,aes(x=Year,y=bm.stock/1e3)) +
                                                    facet_wrap(~Stock.short,scales='free_y') + scale_x_continuous(name='',breaks = few.breaks)+ 
                                                    scale_y_continuous(name="Biomass Pool (thousands of tonnes)",limits=c(0,NA)) + theme(legend.position = 'none') #+
#guides(colour = guide_legend(nrow = 7))

save_plot(filename = paste0(repo.loc,"/Figures/BU/Simulation_quantile_stock_K.png"),sim.stock.pool.quant.plt,base_height = 10,base_width = 20)

# Trophic level quantile plots

sim.TL.pool.quant.plt <- ggplot(quants.pool.tl) + geom_line(aes(x=Years+max(trophic.bm.best$Year)-1,y=mn/1e6),linewidth=1.5,alpha=0.2) + 
                                            geom_ribbon(aes(x=Years+max(bm.best$Year)-1,ymax=UCI/1e6,ymin=LCI/1e6),linewidth=2,alpha=0.2,fill='blue') + 
                                            geom_line(data=trophic.bm.best,aes(x=Year,y=bm.tl/1e6)) +
                                            facet_wrap(~troph.cat,scales='free_y') +
                                            scale_x_continuous(name='',breaks = few.breaks)+ 
                                            scale_y_continuous(name="Biomass Pool (millions of tonnes)",limits=c(0,NA)) + theme(legend.position = 'none') #
#guides(colour = guide_legend(nrow = 7))

save_plot(filename = paste0(repo.loc,"/Figures/BU/Simulation_quantile_TL_K.png"),sim.TL.pool.quant.plt,base_height = 10,base_width = 20)


# Ecosystem quantile plots
sim.com.pool.quant.plt <- ggplot(quants.pool.eco) +  geom_line(aes(x=Years+max(com.tot.bm.best$Year)-1,y=mn/1e6),linewidth=1.5,alpha=0.2) + 
                                               geom_ribbon(aes(x=Years+max(com.tot.bm.best$Year)-1,ymax=UCI/1e6,ymin=LCI/1e6),linewidth=2,alpha=0.2,fill='blue') + 
                                               geom_line(data=com.tot.bm.best,aes(x=Year,y=bm.eco/1e6)) +
                                               scale_x_continuous(name='',breaks = lotsa.breaks)+ 
                                               scale_y_continuous(name="Biomass Pool (millions of tonnes)",limits=c(0,NA)) + theme(legend.position = 'none') #
#guides(colour = guide_legend(nrow = 7))

save_plot(filename = paste0(repo.loc,"/Figures/BU/Simulation_quantile_eco_K.png"),sim.com.pool.quant.plt,base_height = 10,base_width = 20)



# Comparing TL and ecosystem K going stock by stock with the trophic level and ecosystem K's that I originally made up
# And it's not perfect, but I think for a first pass this work, they keep the characteristics we want in terms of
# correlation and the K's are quite similar to the original ones. For TL3 it is perfect, for 4 and 5 it can be slightly off
# because I aimed to keep the trophic level having the correlation over focusing on getting the K exactly right
# There could be ways to do both I haven't thought of, but think this is ok for now.
 #tst <- sim.pool.stocks |> collapse::fgroup_by(troph.cat,Years,sim) |> collapse::fsummarise(tot.bm = sum(bm.stock))
# tst2 <- left_join(tst,sim.troph.pool,by=c("Years","troph.cat","sim"))
# tst2 <- tst2 |> collapse::fgroup_by(c('Years','sim')) |>collapse::fmutate(com.bm.new = sum(tot.bm)) |> as.data.frame()
# # So they definitely differ stock by stock from the original trophic level splits, but most of the time
# # it is within 20% of the original
# tst2$per.diff <- 100*(tst2$tot.bm - tst2$bm.tl) / tst2$bm.tl
# hist(tst2$per.diff[tst2$troph.cat == 5])
# # They do retain the time series characteristics tho
# pacf(tst2$tot.bm[tst2$troph.cat == 3 & tst2$sim == 1])
# # Ecosystem is within 5% way more than 75% of the time
# summary(100*(tst2$com.bm.new - tst2$bm.eco)/tst2$bm.eco)
# # And still has characteristics we want
# pacf(tst2$com.bm.new[tst2$sim == 1][1:n.yrs.proj])
# # OK I think we're good to here, have simulated biomass time series :-)
#ggplot(sim.pool.stocks |> collapse::fsubset(sim == 1)) + geom_line(aes(x=Years,y=bm.stock,group=Stock,color=troph.cat)) + scale_y_log10()
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
#
first.year <- min(bm.best$Year,na.rm=T)
last.year <- 2022
n.years <- length(first.year:last.year)
years <- (last.year+1):(last.year+n.yrs.proj)
# Take the biomass data for the north sea and subset it to the years we have data
init.stock.bm <- bm.best |> collapse::fsubset(Year == last.year) 
# Get the initial ecosystem biomass..
init.com.bm <- sum(init.stock.bm$bm.stock)
init.tl.bm <- init.stock.bm |> collapse::fgroup_by(troph.cat) |> collapse::fsummarise(bm.tl = sum(bm.stock))
 
# Get the average weight of the fish in the stocks so we can go from biomass to abundance for the model
# FIX: NOT SURE I NEED THIS ANYMORE This could definitely be done more sophisisticatedly!
av.wgt <- bm.best |> collapse::fgroup_by(Stock,troph.cat) |> collapse::fsummarise(mn.wgt = mean(avg.weight,na.rm=T))
# FIX: NOT SURE I NEED THIS ANYMORE Let's try getting the most recent year weight to go from biomass to numbers as average may be somewhat misleading
# So here the idea is that the most recent years 
av.wgt <- bm.best |> dplyr::group_by(Stock,troph.cat) |> filter(row_number() >= (n() ))
av.wgt <- data.frame(Stock = av.wgt$Stock,troph.cat = av.wgt$troph.cat,mn.wgt = av.wgt$avg.weight)
# For some debugging, if still here you can delete I'm sure
#count = 0


# Before jumping into the simulation, get the management strategy for each stock set....

# Here is where we insert the management to set the catch.
# If you don't do anything, the default is to fish it using the 
# expoitation from the time series as the RR
# the LRP is 40% of median
# and the USR is the median
# No fishing below the LRP, and fishing turns down linearily between the USR and LRP.
ex.dat <- NULL

for(ss in stock.eco)
{
  if(!is.null(manage))
  {
    
    ms <- manage[manage$Stock == ss,]
    bm.n.er.hist <- bm.best[bm.best$Stock == ss,]
    ex.dat[[ss]] <- data.frame(Stock=ss,lrp = NA, urp = NA, rr= NA,er.mn = NA,er.below.lrp = ms$er.below.lrp,er.sd = ms$er.sd,use.hcr=ms$use.hcr)
    ex.dat[[ss]]$fm.below.lrp <- 1-exp(-ex.dat[[ss]]$er.below.lrp)
    # Get the average eploitation, this is a jumble of logic...
    # If we have entered an exploitation rate we multiply that by the relative exploitation to get the average exploitation
    # These could in theory be set too high (since er is a proportion) and crash everything...
    
    if(!is.na(as.numeric(ms$er))) 
    {
      ex.dat[[ss]]$fm.mn <- ms$relative.er * (1-exp(-as.numeric(ms$er)))
      ex.dat[[ss]]$er.mn <- -log(1-ex.dat[[ss]]$fm.mn)
    }
    # If we wanted the relative exploitation but didn't set the er, then we are setting the  exploitation to the long term median exploitation observed
    # and adjusting that up or down using the relative er. This with a relative er of 1 is the same er as the NULL scenario
    
    if(is.na(as.numeric(ms$er))) {
      ex.dat[[ss]]$fm.mn <- ms$relative.er * (1-exp(-match.fun(ms$er.fun)(bm.n.er.hist$er,na.rm=T)))
      ex.dat[[ss]]$er.mn <- -log(1-ex.dat[[ss]]$fm.mn)
    }
    
    # Are you using a harvest control rule for this stock? Note the removal reference point is set as instantaneous on expoitation rate
    # doing so because it makes sure that exploitation never goes above 1.
    if(ms$use.hcr== F) 
    {
      ex.dat[[ss]]$lrp <- 0
      ex.dat[[ss]]$urp <- 0
      # Note here that if you have your relative er set > 1, the removal reference will be higher than otherwise
      # which makes sense since the fm.mn is the target F for the fishery above the urp.
      ex.dat[[ss]]$rr <- ex.dat[[ss]]$fm.mn
    }
    
    if(ms$use.hcr== T) 
    {
      ex.dat[[ss]]$lrp <- ms$lrp * match.fun(ms$rp.fun)(bm.n.er.hist$bm.stock,na.rm=T)
      ex.dat[[ss]]$urp <- ms$urp * match.fun(ms$rp.fun)(bm.n.er.hist$bm.stock,na.rm=T)
      ex.dat[[ss]]$rr <- ex.dat[[ss]]$fm.mn
    }
    
  } # end if(is.null(exploit))
  # If management is blank, then set it up with these defaults (vaguely Canadian PA approach).
  if(is.null(manage))
  {
    #stk.ex.dat <- exploit[exploit$Stock == s,] 
    #bm.n.er.hist <- bm.best[bm.best$Stock == s,]
    ex.dat[[ss]] <- data.frame(Stock == ss,lrp = NA, urp = NA, rr= NA,er.mn = NA,er.below.lrp = NA)
    # If no data, we'll make the lrp be 40% of the median historic biomass
    ex.dat[[ss]]$lrp <- 0.4*median(bm.n.er.hist$bm.stock,na.rm=T)
    # If no data, we'll make the urp be 100% of the median historic biomass
    ex.dat[[ss]]$urp <- median(bm.n.er.hist$bm.stock,na.rm=T)
    # If no data, we'll make the rr be the median historic fishing mortality, ths is set as instantaneous!
    ex.dat[[ss]]$rr <- 1-exp(-median(bm.best$er[bm.best$Stock == ss],na.rm=T))
    # If we didn't set the rr, we'll make it be the rr
    ex.dat[[ss]]$er.mn <- -log(1-ex.dat[[ss]]$rr)
    ex.dat[[ss]]$fm.mn <- ex.dat[[ss]]$rr
    # If we didn't set the exploitation rate below the lrp, we'll make it 0
    ex.dat[[ss]]$er.below.lrp <- 0
    ex.dat[[ss]]$er.sd <- 0.1
  } # end if(is.null(exploit))
} # end stock loop





# So everything will need to get wrapped up in a simulation loop
res.ts <- NULL
ts.unpack <- NULL
Pool.realized <- NULL
# define the trophic levels we have.


for(j in 1:n.sims)
{
  if(j == 1) count <- 0 
  Ks <- NULL
  K.tmp <- NULL
  results <- NULL
  st.time <- Sys.time()
  
  for(t in 1:n.yrs.proj)
  {
    K.tmp1 <- NULL
    # Now we need to run this through each trophic level, for bottom up we go 3, then 4, then 5....
    for(tl in troph.levels)
    {
      if(tl =="Low") tl.pool <- bm.sim.3[[j]][t]
      # This gets the Ks
      if(tl != "Low") tl.pool.space.prop <- Ks[[j]]$pool.space[Ks[[j]]$tl==tl & Ks[[j]]$Years==t]
      
      bm.stash <- NULL
      tl.stocks <- unique(sim.pool.stocks$Stock[sim.pool.stocks$troph.cat ==tl])
      
      # Now get the stock biomass from last year.
      if(t ==1)  tl.bm.last <- init.tl.bm[init.tl.bm$troph.cat==tl,]
      
      # Then we'll need to get these from the model simulations.
      #if(t==2) 
      if(t > 1)
      {
        # Use the handy av.wgt data.frame I made above
        bm.stocks <- data.frame(bm = NA,meta.dat)
        for(s in stock.eco) bm.stocks$bm[bm.stocks$Stock == s] <-  results[[s]]$bm[results[[s]]$Years == t-1]
        bm.stocks$bm <- bm.stocks$bm
        stock.bm.last <- bm.stocks
        com.bm.last <- sum(bm.stocks$bm)
        tl.bm.last <- bm.stocks |> collapse::fsubset(troph.cat ==tl) |> collapse::fsummarise(bm.tl = sum(bm))
      }  
      # Now we need to figure out what K space is available for each stock within the trophic level.
      # I'm going to base the K space for t+1 on the biomass available in the higher trophic level
      # in year t along with the 'transfer efficiency', which is currently the proportion of the total 
      # ecosystem biomass that the trophic level gets.
      # We can then adjust the stock K's by the available K space in each stock
      # So this is the K space available in a given trophic level in a year
      if(tl =="Low") 
      {
        #browser()
        base.com.pool.tmp <- sim.com.pool |> collapse::fsubset(sim == j & Years ==t)
        # Get the trophic level 3 only since we are going bottom up
        base.tl.pool.tmp <- sim.troph.pool |> collapse::fsubset(sim == j & Years ==t & troph.cat==tl)
        base.stock.pool.tmp <- sim.pool.stocks |> collapse::fsubset(sim == j & Years ==t & troph.cat==tl)
        #base.stock.pool.tmp$prop.bm.stock <- base.stock.pool.tmp
        base.stock.pool.tmp$tl.pool <- tl.pool 
        base.tl.pool.tmp$prop.pool.space <- base.tl.pool.tmp$bm.tl/tl.bm.last$bm.tl
        # We can then adjust the stock K's by the available K space in each stock
        base.stock.pool.tmp$pool.space <- NA
        base.stock.pool.tmp$pool.space <- base.stock.pool.tmp$bm.stock[base.stock.pool.tmp$troph.cat ==tl] * 
          (base.tl.pool.tmp$prop.pool.space-1)
        base.stock.pool.tmp$adj.pool <- base.stock.pool.tmp$bm.stock + base.stock.pool.tmp$pool.space
      } # end if tl =="Low"
      
      if(tl !="Low") 
      {
        #
        base.stock.pool.tmp <- sim.pool.stocks |> collapse::fsubset(sim == j & Years ==t & troph.cat == tl)
        base.tl.pool.tmp <- sim.troph.pool |> collapse::fsubset(sim == j & Years ==t & troph.cat == tl)
        base.stock.pool.tmp$tl.pool <- base.tl.pool.tmp$bm.tl 
        # Get the new trophic level values right
        base.tl.pool.tmp$prop.pool.space <- tl.pool.space.prop
        base.tl.pool.tmp$adj.pool <- base.tl.pool.tmp$prop.pool.space * base.tl.pool.tmp$bm.tl
        base.tl.pool.tmp$pool.space <- base.tl.pool.tmp$adj.pool - base.tl.pool.tmp$bm.tl
        
        # Now get the stock right, all we have to do is multiply the bm.stock
        # by the prop.pool.space (i.e. what the proportinonal chance in the avilable K-space is)
        base.stock.pool.tmp$adj.pool <- base.tl.pool.tmp$prop.pool.space * base.stock.pool.tmp$bm.stock
        base.stock.pool.tmp$pool.space <- base.stock.pool.tmp$adj.pool - base.stock.pool.tmp$bm.stock
        
      }
      # So now I have Carrying Capacities that take up (or lose) any available K space.
      # Now we can convert these to numbers using the historic 'average weight' of the stocks, to avoid complication
      # I'm just using the average of the average weight for each stock...
      #
      
      base.stock.pool.tmp <- left_join(base.stock.pool.tmp,av.wgt,by=c("Stock","troph.cat"))
      # And now we can get a K in numbers....
      base.stock.pool.tmp$adj.pool.num <- base.stock.pool.tmp$adj.pool/base.stock.pool.tmp$mn.wgt
      
      for(s in tl.stocks)
      {
        count <- count + 1
        #stock.lambdas <- bm.best[bm.best$Stock == s,]
        bm.ts.stock <- bm.best[bm.best$Stock == s,]  
        pred.mod <- mod.pred[[s]]
      
      if(t == 1) 
      { 
        bm.start <- bm.ts.stock$bm.stock[bm.ts.stock$Year == last.year]
        results[[s]] <- data.frame(bm = c(bm.start,rep(NA,n.yrs.proj)),
                                   removals = NA,
                                   ex.rate = NA,
                                   Stock = s,
                                   sim= j,
                                   lambda = NA,
                                   Years=0:n.yrs.proj,
                                   troph.cat = tl,
                                   K.bm = NA)
        
      } else{ bm.start <- results[[s]]$bm[results[[s]]$Years == t-1]}
      
      
     
      # Since I have Years and sim recorded, I should just be able to recursivly rbind this...
      
      # This should do the trick for creating this for the very first stop and then updating thereafter.
      # if(t ==1 & j ==1 & count ==1) 
      # {
      #   base.stock.pool <- base.stock.pool.tmp
      #   base.tl.pool <- base.tl.pool.tmp
      #   base.com.pool <- base.com.pool.tmp
      # } else {
      #   base.stock.pool <- rbind(base.stock.pool,base.stock.pool.tmp)
      #   base.tl.pool <- rbind(base.tl.pool,base.tl.pool.tmp)
      #   base.com.pool <- rbind(base.com.pool,base.com.pool.tmp)
        
      #} # end the else...
      #
      l.v.h <- 0.4
      #if(s == "ICES HAWG_NS 4,3a,7d_Clupea_harengus")  l.v.h <- 0.6 # DK Note, using 0.6 for herring stock didn't decline below 50% in this time period.
      #if(s == "ICES WGNSSK_NS 4_Scopthalmus_maximus")  l.v.h <- 0.9 # # DK Note, using 0.9 for this stock because it only declined to 66% of max in time period.
      #if(t==25) 
      cur.pool <- base.stock.pool.tmp$adj.pool[base.stock.pool.tmp$Stock ==s]
      init.pool <- base.stock.pool.tmp$bm.stock[base.stock.pool.tmp$Stock ==s]
      
      # Here is where we insert the management to set the catch.
      # Note this is using all the data, not just the data in recent years to calculate the RPs and RR...
      #If you don't do anything, the default is to fish it using the 
      # expoitation from the time series as the RR
      # the LRP is 40% of median
      # and the USR is the median
      # No fishing below the LRP, and fishing turns down linearily between the USR and LRP.
      #browser()
      er <- proj.catch.eqn(dat =  ex.dat[[s]],bm = bm.start)
      
      res <- pop.dam(stock.dat = bm.ts.stock,stock.pool = cur.pool,bm.start = bm.start, catch = er, stock = s,
                     low.vs.high = l.v.h,method=method,mod.preds = pred.mod)
      
 
      results[[s]]$bm[results[[s]]$Years==t] <- res$tst.res
      results[[s]]$removals[results[[s]]$Years==t-1] <- res$removals
      results[[s]]$ex.rate[results[[s]]$Years==t-1] <- res$ex.rate
      results[[s]]$lambda[results[[s]]$Years==t-1] <- res$lambda
      results[[s]]$K.bm[results[[s]]$Years==t-1] <- cur.pool
      #
      bm.stash[[s]] <- data.frame(pool.stock = res$pool,bm.stock = res$tst.res,pool.space = res$pool - res$tst.res,pool.init =init.pool)
      } # end stock loop by trophic level
      # How this is set up we have the biomass after fishing impacting the K in that year
      # We might want to have this impact the following years K.
      # K for the next trophic level
      if(tl =="Low") 
      {
        #
        bms <- do.call('rbind', bm.stash)
        # This way takes the proportion of k-space available (or missing) and removes it from the trophic levels below
        # i.e. if the lower trophic levels isn't a K, then the higher trophic levels are missing food
        # so 90% of their initially allocated K (I have been debating using the K-space, but I think is problematic 
        # in ways that don't make biological sense to propogate up the chain, I think
        # as we move up the chain means everyone below is hurting too.
        # Doing it this way to ensure that the K's remain positive (by adding actual space you can get negatives....)
        next.tl.pool.space.prop <- (sum(bms$bm.stock)/sum(bms$pool.init))#/mean(sim.tl.4.to.5[[j]]) # Old way this was used to give them the tl efficicnecy, worked badly...
        #next.tl.pool <- sum(do.call('rbind',bm.stash))*sim.tl.3.to.4[[j]][t]
      }
      if(tl =="Medium") 
      {
        #if(t%in% c(1,10,20,50)) 
        bms <- do.call('rbind', bm.stash)
        # This way takes the proportion of k-space available (or missing) and removes it from the trophic levels below
        # i.e. if the lower trophic levels isn't a K, then the higher trophic levels are missing food
        # so 90% of their initially allocated K (I have been debating using the reallocated K-space, but I think is problematic as the collapse of the lower levels makes the K space collapse
        # in ways that don't make biological sense to propogate up the chain, I think
        # as we move up the chain means everyone below is hurting too.
        # Doing it this way to ensure that the K's remain positive (by adding actual space you can get negatives....)
        next.tl.pool.space.prop <- (sum(bms$bm.stock)/sum(bms$pool.init))#/mean(sim.tl.4.to.5[[j]]) # Old way this was used to give them the tl efficicnecy, worked badly...
        #next.tl.pool <- sum(do.call('rbind',bm.stash))*sim.tl.3.to.4[[j]][t]
      }
      
      if(tl != "High") 
      {
        mv.to.next <- which(troph.levels == tl)+1
        Ks[[j]] <- rbind(Ks[[j]],data.frame(pool.space.prop=next.tl.pool.space.prop,tl = troph.levels[mv.to.next],Years=t,sim=j))
      }
      K.tmp1[[tl]] <- base.stock.pool.tmp
    }# End trophic level loop.
    #browser()
    K.tmp[[t]] <- do.call("rbind",K.tmp1)
  #
    
  } # end the t looping through each year.

res.ts[[j]] <- do.call("rbind",results)
Pool.realized[[j]] <- do.call("rbind",K.tmp)
#ggplot(base.stock.pool) + geom_line(aes(x= Years,y=bm.stock,group=Stock,color=Stock)) + facet_wrap(~troph.cat) + scale_y_log10()
  #
  # Unpack the results
  #st.time <- Sys.time()
  #ts.unpack[[j]] <- do.call('rbind',res.ts)
  
#ggplot(ts.unpack[[j]]) + geom_line(aes(x= Years,y=abund,group=Stock,color=Stock)) + facet_wrap(~troph.cat) + scale_y_log10()
  
  
  # Pop a note when done each simulation
  timer <- Sys.time() - st.time
  print(paste("Simulation ", j))
  print(signif(timer,digits=2))
  
} # end n.sims

# Unpack all the results.
ts.final <- do.call("rbind",res.ts)
Pool.real <- do.call("rbind",Pool.realized)
#browser()

#ggplot(ts.final) + geom_line(aes(x= Years,y=abund,group=sim,color=sim)) + facet_wrap(~Stock) + scale_y_log10()


#ts.final$fm <- ts.final$removals/ts.final$abund
#av.wgt$troph.cat <- as.numeric(av.wgt$troph.cat)
ts.final <- left_join(ts.final,av.wgt,by=c("Stock","troph.cat"))
ts.final$troph.cat <- factor(ts.final$troph.cat, levels = c("Low", "Medium" ,"High"),labels = c("Low", "Medium", "High"))

#ts.final$biomass <- ts.final$abund*ts.final$mn.wgt
#r.final <- do.call("rbind",r.unpack)
quants <- ts.final |>  collapse::fgroup_by(Years,Stock,troph.cat) |> collapse::fsummarize(L.50 = quantile(bm,probs=c(0.25),na.rm=T),
                                                                                          med = median(bm,na.rm=T),
                                                                                          U.50 = quantile(bm,probs=c(0.75),na.rm=T),
                                                                                          er.L.50 = quantile(ex.rate,probs=c(0.25),na.rm=T),
                                                                                          er.med = median(ex.rate,na.rm=T),
                                                                                          er.U.50 = quantile(ex.rate,probs=c(0.75),na.rm=T),
                                                                                          rem.L.50 = quantile(removals,probs=c(0.25),na.rm=T),
                                                                                          rem.med = median(removals,na.rm=T),
                                                                                          rem.U.50 = quantile(removals,probs=c(0.75),na.rm=T))

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



# Two simple plots. 
p.sims <- ggplot(ts.final) + geom_line(aes(x=Years+max(bm.best$Year),y=bm/1e3,group = sim,color=sim),alpha=0.8) +
                              facet_wrap(~interaction(Stock.short, troph.cat , sep = " TC: ",drop=T),scales='free_y') +
                              geom_line(data=bm.best,aes(x=Year,y=bm.stock/1e3)) +
                              scale_x_continuous(name='',breaks = few.breaks) +
                              scale_y_continuous(name = "Biomass (thousands of tonnes)") + 
                              theme(legend.position = 'none') 

save_plot(paste0(repo.loc,"/Figures/BU/raw_biomass_trends.png"),p.sims,base_height = 12,base_width = 24)

# How about a removals and exploitation rate time series.
#browser()
p.removals <- ggplot(ts.final) + geom_line(aes(x=Years+max(bm.best$Year),y=removals/1000,group = sim,color=sim),alpha=0.8) +
                                  facet_wrap(~interaction(Stock.short, troph.cat , sep = " TC: ",drop=T),scales='free') +
                                  geom_line(data=bm.best,aes(x=Year,y=catch/1000)) +
                                  scale_x_continuous(name='',breaks = few.breaks) +
                                  scale_y_continuous(name = "Removals (thousands of tonnes)") + 
                                  theme(legend.position = 'none') 

save_plot(paste0(repo.loc,"/Figures/BU/raw_removals_trends.png"),p.removals,base_height = 12,base_width = 24)

p.er <- ggplot(ts.final) + geom_line(aes(x=Years+max(bm.best$Year),y=ex.rate,group = sim,color=sim),alpha=0.8) +
                              facet_wrap(~interaction(Stock.short, troph.cat , sep = " TC: ",drop=T),scales='free') +
                              geom_line(data=bm.best,aes(x=Year,y=er)) +
                              scale_x_continuous(name='',breaks = few.breaks) +
                              scale_y_continuous(name = "Exploitation rate (proportional)") + 
                              theme(legend.position = 'none') 

save_plot(paste0(repo.loc,"/Figures/BU/raw_exploitation_rates.png"),p.er,base_height = 12,base_width = 24)


# 
p.sims.quants <- ggplot(quants) + geom_line(aes(x=Years+max(bm.best$Year),y=med/1e3,group=Stock,color=spec.tl),linewidth=2) +  
                                  #geom_line(data=bm.best,aes(x=Year,y=bm.stock,group=Stock,color=spec.tl)) +
                                  facet_wrap(~troph.cat) +  scale_y_log10(name="Biomass (thousands of tonnes)") +   theme(legend.position = 'top') +
                                  guides(colour = guide_legend(nrow = 5)) + scale_color_manual(values=pal) +
                                  scale_x_continuous(name='',breaks = few.breaks)
#geom_ribbon(data=quants, aes(x=Years,ymax=U.50,ymin = L.50),alpha=0.5,fill='blue',color='blue') 
save_plot(paste0(repo.loc,"/Figures/BU/Quantile_biomass_trends.png"),p.sims.quants,base_height = 8,base_width = 16)

p.sims.quants.by.stock <- ggplot(quants) + geom_line(aes(x=Years+max(bm.best$Year),y=med/1e3),linewidth=1.25,color='grey50') + 
                                            geom_ribbon(aes(x=Years+max(bm.best$Year),ymax=U.50/1e3,ymin=L.50/1e3),fill='firebrick2',alpha=0.2)+
                                            geom_line(data=bm.best,aes(x=Year,y=bm.stock/1e3),color='black')+#,linetype = 'dashed') +
                                            facet_wrap(~interaction(Stock.short, troph.cat , sep = " TC: ",drop=T),scales='free_y') +
                                            scale_y_continuous(name="Biomass (thousands of tonnes)") +   
                                            theme(legend.position = 'top') +
                                            guides(colour = guide_legend(nrow = 5)) +
                                            scale_x_continuous(name='',breaks = few.breaks)
#geom_ribbon(data=quants, aes(x=Years,ymax=U.50,ymin = L.50),alpha=0.5,fill='blue',color='blue') 
save_plot(paste0(repo.loc,"/Figures/BU/Quantile_biomass_trends_by_stock.png"),p.sims.quants.by.stock,base_height = 12,base_width = 24)

p.sims.rems.by.stock <- ggplot(quants) + geom_line(aes(x=Years+max(bm.best$Year),y=rem.med/1e3),linewidth=1.25,color='grey50') + 
                                          geom_ribbon(aes(x=Years+max(bm.best$Year),ymax=rem.U.50/1e3,ymin=rem.L.50/1e3),fill='firebrick2',alpha=0.2)+
                                          geom_line(data=bm.best,aes(x=Year,y=catch/1e3)) +
                                          facet_wrap(~interaction(Stock.short, troph.cat , sep = " TC: ",drop=T),scales='free_y') +
                                          scale_y_continuous(name="Removals (thousands of tonnes)") +   
                                          theme(legend.position = 'top') +
                                          guides(colour = guide_legend(nrow = 5)) +
                                          scale_x_continuous(name='',breaks = few.breaks)
#geom_ribbon(data=quants, aes(x=Years,ymax=U.50,ymin = L.50),alpha=0.5,fill='blue',color='blue') 
save_plot(paste0(repo.loc,"/Figures/BU/Quantile_removals_by_stock.png"),p.sims.rems.by.stock,base_height = 12,base_width = 24)
#
p.sims.er.by.stock <- ggplot(quants) + geom_line(aes(x=Years+max(bm.best$Year),y=er.med),linewidth=1.25,color='grey50') + 
                                        geom_ribbon(aes(x=Years+max(bm.best$Year),ymax=er.U.50,ymin=er.L.50),fill='firebrick2',alpha=0.2)+
                                        geom_line(data=bm.best,aes(x=Year,y=er)) +
                                        facet_wrap(~interaction(Stock.short, troph.cat , sep = " TC: ",drop=T),scales='free_y') +
                                        scale_y_continuous(name="Exploitation rate (proportional)") +   
                                        theme(legend.position = 'top') +
                                        guides(colour = guide_legend(nrow = 5)) +
                                        scale_x_continuous(name='',breaks = few.breaks)
#geom_ribbon(data=quants, aes(x=Years,ymax=U.50,ymin = L.50),alpha=0.5,fill='blue',color='blue') 
save_plot(paste0(repo.loc,"/Figures/BU/Quantile_er_by_stock.png"),p.sims.er.by.stock,base_height = 12,base_width = 24)


# Clean up the names of the biomass pools to distinguish them from the actual biomasses
names(Pool.real) <- c("Years","sim","Stock","troph.cat","prop.tl.pool", "cor.prop.tl.pool", "pool.init.stock", "trophic","Species","color",   
                   "spec.tl", "Stock.short","common","tl.pool","pool.space","pool.real","mn.wgt", "pool.real.num")
names(sim.pool.stocks) <- c("Years","sim","Stock","troph.cat","prop.tl.pool", "cor.prop.tl.pool", "pool.init.stock", "trophic","Species","color",   
                         "spec.tl", "Stock.short","common")
names(sim.troph.pool) <- c("Years","sim","pool.tl","troph.cat","pool.com","prop.com.pool")
names(sim.com.pool) <- c("com.pool","Years","sim")

return(list(sim.quantiles = quants,
            sim.ts = ts.final,
            past.bm = bm.best,
            sim.pool.stocks = sim.pool.stocks,
            sim.troph.pool = sim.troph.pool,
            sim.com.pool = sim.com.pool,
            Pool.realized = Pool.real
            ))
} # end function