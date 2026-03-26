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
          "https://raw.githubusercontent.com/dave-keith/Multispecies_model/main/Scripts/NS_catch_function.R"
          
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

# Community and trophic biomasses
com.tot.bm <- bm.best |> collapse::fgroup_by(Year) |> 
  collapse::fsummarize(num.eco = sum(num.stock),bm.eco = sum(bm.stock)+sum(catch))

trophic.bm <- bm.best |> collapse::fgroup_by(Year,troph.cat) |> 
  collapse::fsummarize(num.tl = sum(num.stock),bm.tl = sum(bm.stock)+sum(catch))
#trophic.bm$troph.cat <- factor(trophic.bm$troph.cat,levels = c("Low","Medium","High"))

tl.com.bm <- left_join(trophic.bm,com.tot.bm,by="Year")
tl.com.bm$prop.bm.tl <- tl.com.bm$bm.tl/tl.com.bm$bm.eco
tl.com.bm$prop.num.tl <- tl.com.bm$num.tl/tl.com.bm$num.eco


# The 'transfer efficiency' between our trophic levels
tl.3.to.4 <- bm.best$prop.bm.tl[bm.best$troph.cat=="Medium"][1:n.years]/bm.best$prop.bm.tl[bm.best$troph.cat=="Low"][1:n.years]
tl.4.to.5 <- bm.best$prop.bm.tl[bm.best$troph.cat=="High"][1:n.years]/bm.best$prop.bm.tl[bm.best$troph.cat=="Medium"][1:n.years]
tl.3.to.5 <- bm.best$prop.bm.tl[bm.best$troph.cat=="High"][1:n.years]/bm.best$prop.bm.tl[bm.best$troph.cat=="Low"][1:n.years]


# So now we want to look at stock level within a trophic level
# add some colors...

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
meta.dat$troph.cat <- meta.dat$troph.cat

colors <- distinct(bm.best, spec.tl, color)
pal <- colors$color
names(pal) <- colors$spec.tl



stock.prop.bm.plt <- ggplot(bm.best) + geom_line(aes(x=Year,y=prop.bm.stock.tl,group = Stock,color=spec.tl),linewidth=2) + 
  facet_wrap(~troph.cat) + guides(colour = guide_legend(nrow = 5)) + theme(legend.position = 'top',legend.title = element_blank()) +
  scale_y_log10(name= "Proportion of biomass Pool",n.breaks=10) + scale_x_continuous(name="",breaks=seq(1970,2200,by=5))+
  scale_color_manual(values=pal)
save_plot(paste0(repo.loc,"/Figures/NI/Historic_Prop_Biomass_ns_by_stock.png"),stock.prop.bm.plt,base_height = 8,base_width = 15)

stock.bm.plt <- ggplot(bm.best) + geom_line(aes(x=Year,y=bm.stock/1e3,group = Stock,color=spec.tl),linewidth=2) + 
  facet_wrap(~troph.cat) + scale_x_continuous(name="",breaks=seq(1970,2200,by=5))+
  scale_y_log10(name = "Biomass Pool (thousands of tonnes)",n.breaks=7) + theme(legend.position = 'top',legend.title = element_blank()) +
  guides(colour = guide_legend(nrow = 5)) + scale_color_manual(values=pal)
save_plot(paste0(repo.loc,"/Figures/NI/Historic_Biomass_ns_by_stock.png"),stock.bm.plt,base_height = 8,base_width = 15)

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
#tl.3.4.cor <- ccf(trophic.bm.best$bm.tl[trophic.bm.best$troph.cat =="Low"],trophic.bm.best$bm.tl[trophic.bm.best$troph.cat =="Medium"],plot = F)
#tl.3.5.cor <- ccf(trophic.bm.best$bm.tl[trophic.bm.best$troph.cat =="Low"],trophic.bm.best$bm.tl[trophic.bm.best$troph.cat =="High"],plot = F)
#tl.4.5.cor <- ccf(trophic.bm.best$bm.tl[trophic.bm.best$troph.cat =="Medium"],trophic.bm.best$bm.tl[trophic.bm.best$troph.cat =="High"],plot = F)
# Looking at proportions, need to stew a bit on this because there is necessarily some 
# correlation built into proportions, but what is interesting is that
# the correlation strength is really really high between TL 3 and TL 4, it is weaker (more diffuse really) at TL 5
# and there is no correlation between 4 and 5
tl.3.4.prop.cor <- ccf(bm.best$prop.bm.tl[bm.best$troph.cat =="Low"][1:n.years],bm.best$prop.bm.tl[bm.best$troph.cat =="Medium"][1:n.years],plot = F)
tl.3.5.prop.cor <- ccf(bm.best$prop.bm.tl[bm.best$troph.cat =="Low"][1:n.years],bm.best$prop.bm.tl[bm.best$troph.cat =="High"][1:n.years],plot = F)
tl.4.5.prop.cor <- ccf(bm.best$prop.bm.tl[bm.best$troph.cat =="Medium"][1:n.years],bm.best$prop.bm.tl[bm.best$troph.cat =="High"][1:n.years],plot = F)


# So now really what I need to do first is make a quick simulation that gets me ecosystem K, trophic level K, and stock K
# once I have those then we just run the models :-)
#mn.com.bm <- mean(com.tot.bm.best$bm.eco)
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
tl.4.5.prop.bm <- bm.best$bm.tl[bm.best$troph.cat =="High"][1:n.years]/(bm.best$bm.tl[bm.best$troph.cat =="Medium"][1:n.years]+
                                                                          bm.best$bm.tl[bm.best$troph.cat =="High"][1:n.years])
tl.4.5.prop.4.5.bm <- pacf(tl.4.5.prop.bm,plot = F)
#NEW: proportion of 5 in the ecosystem
tl.5.prop.bm.ts <- bm.best$prop.bm.tl[bm.best$troph.cat=="High"][1:n.years]
tl.5.prop.pacf <- pacf(tl.5.prop.bm.ts, plot=F)
tl.5.prop.bm.lag.1 <- tl.5.prop.pacf$acf[1]
tl.5.prop.bm.lag.2 <-tl.5.prop.pacf$acf[2]
#NEW: transfer efficiencies
#didn't make a time series variable, because the transfer efficiencies are already isolated
tl.3.to.4.pacf <- pacf(tl.3.to.4, plot=F)
tl.3.to.4.lag.1 <- tl.3.to.4.pacf$acf[1]
tl.3.to.4.lag.2 <- tl.3.to.4.pacf$acf[2]
tl.4.to.5.pacf <- pacf(tl.4.to.5, plot=F)
tl.4.to.5.lag.1 <- tl.4.to.5.pacf$acf[1]
tl.4.to.5.lag.2 <- tl.4.to.5.pacf$acf[2]


troph.levels <- sort(unique(bm.best$troph.cat))

# Get necessary data on logit scale
tl.3.logit <- logit(tl.3.prop.bm.ts)
tl.4.5.logit <- logit(tl.4.5.prop.bm)
#NEW
tl.5.logit <- logit(tl.5.prop.bm.ts)
tl.3.to.4.logit <- logit(tl.3.to.4)
tl.4.to.5.logit <- logit(tl.4.to.5)

# Starting values for the ecosystem and the proportions, logit needed for arima models with the proportions
start.com.sim <- max.com.bm
start.tl.3.prop.bm <- tl.3.prop.bm.ts[length(tl.3.prop.bm.ts)]
start.tl.3.logit <- tl.3.logit[length(tl.3.logit)]
start.tl.4.5.prop.bm <- tl.4.5.prop.bm[length(tl.4.5.prop.bm)]
start.tl.4.5.logit <- tl.4.5.logit[length(tl.4.5.logit)]
#NEW
start.tl.5.prop.bm <- tl.5.prop.bm.ts[length(tl.5.prop.bm.ts)]
start.tl.5.logit <- tl.5.logit[length(tl.5.logit)]
start.tl.3.to.4 <- tl.3.to.4[length(tl.3.to.4)]
start.tl.3.to.4.logit <- tl.3.to.4.logit[length(tl.3.to.4.logit)]
start.tl.4.to.5 <- tl.4.to.5[length(tl.4.to.5)]
start.tl.4.to.5.logit <- tl.4.to.5.logit[length(tl.4.to.5.logit)]

# Mean values for the trophic levels
mn.tl.3.prop.bm <- mean(tl.3.prop.bm.ts)
mn.tl.3.logit <- mean(tl.3.logit)
mn.tl.4.5.prop.bm <- mean(tl.4.5.prop.bm)
mn.tl.4.5.logit <- mean(tl.4.5.logit)
#NEW
mn.tl.5.prop.bm <- mean(tl.5.prop.bm.ts)
mn.tl.5.logit <- mean(tl.5.logit)
mn.tl.3.to.4 <- mean(tl.3.to.4)
mn.tl.3.to.4.logit <- mean(tl.3.to.4.logit)
m.tl.4.to.5 <- mean(tl.4.to.5)
mn.tl.4.to.5.logit <- mean(tl.4.to.5.logit)

# We could instead use the most recent year as the mean...
#mn.tl.3.prop.bm <- start.tl.3.prop.bm
#mn.tl.3.logit <- start.tl.3.logit
#mn.tl.4.5.prop.bm <- start.tl.4.5.prop.bm
#mn.tl.4.5.logit <- start.tl.4.5.logit
#NEW
#mn.tl.5.prop.bm <- start.tl.5.prop.bm
#mn.tl.5.logit <- start.tl.5.logit
#mn.tl.3.to.4 <- start.tl.3.to.4
#mn.tl.3.to.4.logit <- start.tl.3.to.4.logit
#mn.tl.4.to.5 <- start.tl.4.to.5
#mn.tl.4.to.5.logit <- start.tl.4.to.5.logit

# For the community I am starting it at the maximum and running from there
start.com.diff = start.com.sim- max.com.bm
# Tor the rest I am starting it at where they ended.
start.tl.3.diff <- start.tl.3.logit - mn.tl.3.logit
start.tl.4.5.diff <- start.tl.4.5.logit - mn.tl.4.5.logit
#NEW
start.tl.5.diff <- start.tl.5.logit - mn.tl.5.logit
start.tl.3.to.4.diff <- start.tl.3.to.4.logit - mn.tl.3.to.4.logit
start.tl.4.to.5.diff <- start.tl.4.to.5.logit - mn.tl.4.to.5.logit

# the standard deviations
sd.tl.3.logit <- sd(tl.3.logit)
sd.tl.4.5.logit <- sd(tl.4.5.logit)
#NEW
sd.tl.5.logit <- sd(tl.5.logit)
sd.tl.3.to.4.logit <- sd(tl.3.to.4.logit)
sd.tl.4.to.5.logit <- sd(tl.4.to.5.logit)

# Lag for the Arima model
tl.4.5.prop.bm.lag.1 <- tl.4.5.prop.4.5.bm$acf[1]

sim.pool.stock.lst <- NULL
sim.pools <- NULL
sim.com.bm <- NULL
bm.trophic.pools <- NULL
bm.sim.3 <- NULL
bm.sim.4 <- NULL
bm.sim.5 <- NULL

for(i in 1:n.sims) 
{
  #browser()
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
  sim.tl.3.prop.bm <-inv.logit(mn.tl.3.logit + 
                                 arima.sim(model =list(ar = c(tl.3.prop.bm.lag.1)),n = n.yrs.proj,
                                           n.start =1, start.innov = c(start.tl.3.diff/tl.3.prop.bm.lag.1), 
                                           innov = c(0,rnorm(n.yrs.proj-1,0,sd.tl.3.logit))))
  
  bm.sim.3[[i]] <- sim.tl.3.prop.bm * sim.com.bm[[i]]$bm
  # So this is what is left for 3 and 4
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
  #browser()
  bm.trophic.pools[[i]] <- data.frame(Years = rep(1:n.yrs.proj,3), sim =i,
                                      bm.tl = c(bm.sim.3[[i]],bm.sim.4[[i]],bm.sim.5[[i]]),troph.cat = sort(rep(troph.levels,n.yrs.proj)),
                                      bm.eco = rep(sim.com.bm[[i]]$bm,3))
  bm.trophic.pools[[i]]$prop.bm.tl <- bm.trophic.pools[[i]]$bm.tl/bm.trophic.pools[[i]]$bm.eco
  
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
      #count = count+1
      #browser()
      # Now get the time series for each stock...
        tmp.dat <- bm.best[bm.best$Stock ==s,]
        tmp.cor <- pacf(tmp.dat$prop.bm.stock.tl,plot=F) # Get the correlation, use AR1 and AR2 but no more.
        tmp.cor.lag.1 <- tmp.cor$acf[1]
        #tmp.cor.lag.2 <- tmp.cor$acf[2]
        #tmp.beta <- estBetaParams(mean(tmp.dat$prop.bm.tl),sd(tmp.dat$prop.bm.tl)^2)
        # Logit tranform the proportions and do the ARIMA on the logits
        bm.logit <- logit(tmp.dat$prop.bm.stock.tl)
        start.bm.logit <- bm.logit[length(bm.logit)]
        mn.bm.logit <- mean(bm.logit)
        # DK Note, I thought about using the most recent bm on logit scale as the 'mean' value for the simulation to 
        # start where we finished, but I don't like that behaviour in TL3, so going to use the mean which means the stocks
        # will want to go back to an average value of the 'sharing' of biomass. Shit-canning this will make
        # some of the below unnecessaril complicated.
        #mn.bm.logit <- start.bm.logit
        # And the standard deviation
        sd.bm.logit <- sd(bm.logit)
        diff.bm.logit <- start.bm.logit - mn.bm.logit
        #browser()
        # Then backtransform and everything will stay positive! Just using the AR1 term for these
        # FIX: SEE above comment for where I'm using the AR2, here using the AR2 would give some poor starting values
        # So I'm not comfy doing that (it works by luck in the above for the NS IMHO.)
        
        tmp.prop.bm <- c(inv.logit(arima.sim(model =list(ar = c(tmp.cor.lag.1)),
                                             n.start = 1, start.innov = c(diff.bm.logit/tmp.cor.lag.1),
                                             n = n.yrs.proj+1,innov = c(0,rnorm(n.yrs.proj,0,sd.bm.logit))) + mn.bm.logit))[2:(n.yrs.proj+1)]
        
        sim.pools[[s]] <- data.frame(Years = 1:n.yrs.proj, sim = i,
                                     Stock = s, troph.cat = tl,
                                     prop.bm.stock = tmp.prop.bm,
                                     cor.prop.bm = NA,
                                     bm.stock = bm.trophic.pools[[i]]$bm.tl[bm.trophic.pools[[i]]$troph.cat==tl])

    } # end the stocks loop
    
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
  } # end the trophic level loop
  sim.pool.stock.lst[[i]] <- do.call("rbind",sim.pools)
  
} # end the simulation loop

sim.pool.stocks.tmp <- do.call("rbind",sim.pool.stock.lst)
sim.troph.pool <- do.call("rbind",bm.trophic.pools)
sim.com.pool <- do.call("rbind",sim.com.bm)

# Get the meta data into the K stuff
sim.pool.stocks <- left_join(sim.pool.stocks.tmp,meta.dat,by=c("Stock",'troph.cat'))

# Wrap up the K time series for each simulation
#sim.K.stocks$Species <- substr(sim.K.stocks$Stock,14,100)
quants.pool.stocks <- sim.pool.stocks |> dplyr::group_by(Stock,common,Species,Stock.short,Years,troph.cat) |> dplyr::summarize(mn = mean(bm.stock,na.rm=T),
                                                                                                                         log.mn = mean(log(bm.stock),na.rm=T),
                                                                                                                         med = median(bm.stock,na.rm=T),
                                                                                                                         sd = sd(log(bm.stock),na.rm=T))
quants.pool.stocks$UCI <- exp(quants.pool.stocks$log.mn + quants.pool.stocks$sd)
quants.pool.stocks$LCI <- exp(quants.pool.stocks$log.mn - quants.pool.stocks$sd)
# If this CI goes negative than set it to 1.  probably should do this on the log scale...
quants.pool.stocks$LCI[quants.pool.stocks$LCI < 0] <- 1

quants.pool.tl <- sim.troph.pool  |> collapse::fgroup_by(troph.cat,Years) |> collapse::fsummarize(mn = mean(bm.tl,na.rm=T),
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
#sim.K.stocks$Species <- substr(sim.K.stocks$Stock,14,100)
sim.stock.pool.plt <- ggplot(sim.pool.stocks) + geom_line(aes(x=Years+max(bm.best$Year)-1,y=bm.stock/1e3,group=sim),linewidth=2,alpha=0.2,color='grey50') + 
  geom_line(data=bm.best,aes(x=Year,y=bm.stock/1e3),color='black') +
  facet_wrap(~troph.cat+Stock.short,scales='free_y') +  scale_x_continuous(name='',breaks = few.breaks)+ 
  scale_y_continuous(name="Biomass Pool (thousands of tonnes)",limits=c(0,NA)) + theme(legend.position = 'none') #+

save_plot(filename = paste0(repo.loc,"/Figures/NI/Simulation_stock_K.png"),sim.stock.pool.plt,base_height = 10,base_width = 20)


sim.tl.pool.plt <- ggplot(sim.troph.pool) + geom_line(aes(x=Years+max(bm.best$Year)-1,y=bm.tl/1e6,group=as.factor(sim),color=as.factor(sim))) + 
  geom_line(data=bm.best,aes(x=Year,y=bm.tl/1e6)) +
  facet_wrap(~troph.cat,scales='free_y') + theme(legend.position='none') + scale_x_continuous(name='',breaks = few.breaks) + 
  scale_y_continuous(name="Biomass Pool (millions of tonnes)",limits=c(0,NA))  
save_plot(filename = paste0(repo.loc,"/Figures/NI/Simulation_trophic_K.png"),sim.tl.pool.plt,base_height = 10,base_width = 20)

sim.com.pool.plt <- ggplot(sim.com.pool) + geom_line(aes(x=Years+max(bm.best$Year)-1,y=bm/1e6,group=as.factor(sim),color=as.factor(sim))) +
  geom_line(data=bm.best,aes(x=Year,y=bm.eco/1e6)) + scale_x_continuous(name='',breaks = lotsa.breaks)+ 
  scale_y_continuous(name="Biomass Pool (millions of tonnes)",limits=c(0,NA)) + theme(legend.position = 'none')
save_plot(filename = paste0(repo.loc,"/Figures/NI/Simulation_eco_K.png"),sim.com.pool.plt,base_height = 8,base_width = 11)

# Now get the quantile plots with some uncertainty
#
sim.stock.pool.quant.plt <-  ggplot(quants.pool.stocks) + geom_line(aes(x=Years+max(bm.best$Year)-1,y=mn/1e3),linewidth=1,alpha=0.2) + 
  geom_ribbon(aes(x=Years+max(bm.best$Year)-1,ymax=UCI/1e3,ymin=LCI/1e3),linewidth=1,alpha=0.2,fill='blue') + 
  geom_line(data=bm.best,aes(x=Year,y=bm.stock/1e3)) +
  facet_wrap(~Stock.short,scales='free_y') + scale_x_continuous(name='',breaks = few.breaks)+ 
  scale_y_continuous(name="Biomass Pool (thousands of tonnes)",limits=c(0,NA)) + theme(legend.position = 'none') #+
#guides(colour = guide_legend(nrow = 7))

save_plot(filename = paste0(repo.loc,"/Figures/NI/Simulation_quantile_stock_K.png"),sim.stock.pool.quant.plt,base_height = 10,base_width = 20)

# Trophic level quantile plots

sim.TL.pool.quant.plt <- ggplot(quants.pool.tl) + geom_line(aes(x=Years+max(trophic.bm.best$Year)-1,y=mn/1e6),linewidth=1.5,alpha=0.2) + 
  geom_ribbon(aes(x=Years+max(bm.best$Year)-1,ymax=UCI/1e6,ymin=LCI/1e6),linewidth=2,alpha=0.2,fill='blue') + 
  geom_line(data=trophic.bm.best,aes(x=Year,y=bm.tl/1e6)) +
  facet_wrap(~troph.cat,scales='free_y') +
  scale_x_continuous(name='',breaks = few.breaks)+ 
  scale_y_continuous(name="Biomass Pool (millions of tonnes)",limits=c(0,NA)) + theme(legend.position = 'none') #
#guides(colour = guide_legend(nrow = 7))

save_plot(filename = paste0(repo.loc,"/Figures/NI/Simulation_quantile_TL_K.png"),sim.TL.pool.quant.plt,base_height = 10,base_width = 20)


# Ecosystem quantile plots
sim.com.pool.quant.plt <- ggplot(quants.pool.eco) +  geom_line(aes(x=Years+max(com.tot.bm.best$Year)-1,y=mn/1e6),linewidth=1.5,alpha=0.2) + 
  geom_ribbon(aes(x=Years+max(com.tot.bm.best$Year)-1,ymax=UCI/1e6,ymin=LCI/1e6),linewidth=2,alpha=0.2,fill='blue') + 
  geom_line(data=com.tot.bm.best,aes(x=Year,y=bm.eco/1e6)) +
  scale_x_continuous(name='',breaks = lotsa.breaks)+ 
  scale_y_continuous(name="Biomass Pool (millions of tonnes)",limits=c(0,NA)) + theme(legend.position = 'none') #
#guides(colour = guide_legend(nrow = 7))

save_plot(filename = paste0(repo.loc,"/Figures/NI/Simulation_quantile_eco_K.png"),sim.com.pool.quant.plt,base_height = 10,base_width = 20)



# Comparing TL and ecosystem K going stock by stock with the trophic level and ecosystem K's that I originally made up
# And it's not perfect, but I think for a first pass this work, they keep the characteristics we want in terms of
# correlation and the K's are quite similar to the original ones. For TL3 it is perfect, for 4 and 5 it can be slightly off
# because I aimed to keep the trophic level having the correlation over focusing on getting the K exactly right
# There could be ways to do both I haven't thought of, but think this is ok for now.
 #tst <- sim.K.stocks |> collapse::fgroup_by(troph.cat,Years,sim) |> collapse::fsummarise(tot.bm = sum(bm.stock))
# tst2 <- left_join(tst,sim.troph.K,by=c("Years","troph.cat","sim"))
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
#rowser()
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
      ex.dat[[ss]]$usr <- 0
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
    bm.n.er.hist <- bm.best[bm.best$Stock == ss,]
    #stk.ex.dat <- exploit[exploit$Stock == s,] 
    #bm.n.er.hist <- bm.best[bm.best$Stock == s,]
    ex.dat[[ss]] <- data.frame(Stock == ss,lrp = NA, urp = NA, rr= NA,er.mn = NA,er.below.lrp = NA)
    # If no data, we'll make the lrp be 40% of the median historic biomass
    ex.dat[[ss]]$lrp <- 0.4*median(bm.n.er.hist$bm.stock,na.rm=T)
    # If no data, we'll make the urp be 100% of the median historic biomass
    ex.dat[[ss]]$urp <- median(bm.n.er.hist$bm.stock,na.rm=T)
    # If no data, we'll make the rr be the median historic fishing mortality, ths is set as instantaneous!
    ex.dat[[ss]]$rr <- 1-exp(-median(bm.n.er.hist$er,na.rm=T))
    # If we didn't set the rr, we'll make it be the rr
    ex.dat[[ss]]$er.mn <- -log(1-ex.dat[[ss]]$rr)
    ex.dat[[ss]]$fm.mn <- ex.dat[[ss]]$rr
    # If we didn't set the exploitation rate below the lrp, we'll make it 0
    ex.dat[[ss]]$er.below.lrp <- 0
    ex.dat[[ss]]$er.sd <- 0.1
  } # end if(is.null(exploit))
} # end stock loop


# So everything will need to get wrapped up in a simulation loop
results <- NULL
res.ts <- NULL
#ts.unpack <- NULL
for(j in 1:n.sims)
{
  st.time <- Sys.time()
  
  for(t in 1:n.yrs.proj)
  {
    # Get some starting points. These are for the current year
    base.com.pool.tmp <- sim.com.pool |> collapse::fsubset(sim == j & Years ==t)
    base.tl.pool.tmp <- sim.troph.pool |> collapse::fsubset(sim == j & Years ==t)
    base.stock.pool.tmp <- sim.pool.stocks |> collapse::fsubset(sim == j & Years ==t)
    # Now get the stock biomass from last year.
    if(t ==1)
    {
      stock.bm.last <- init.stock.bm
      stock.bm.last <- stock.bm.last[order(stock.bm.last$troph.cat),]
      com.bm.last <- init.com.bm
      tl.bm.last <- init.tl.bm
    }
    # Then we'll need to get these from the model simulations.
    if(t > 1)
    {
      # Use the handy av.wgt data.frame I made above
      bm.stocks <- data.frame(bm = NA,meta.dat)
      for(s in stock.eco) bm.stocks$bm[bm.stocks$Stock == s] <- results[[s]]$bm[results[[s]]$Years == t-1]
      bm.stocks$bm <- bm.stocks$bm
      stock.bm.last <- bm.stocks
      com.bm.last <- sum(bm.stocks$bm)
      tl.bm.last <- bm.stocks |> collapse::fgroup_by(troph.cat) |> collapse::fsummarise(bm.tl = sum(bm))
    }  
   # Now we need to figure out what K space is available for each stock within the trophic level.
   # First is our Trophic level above the K we have for it.
    # So this is the K space available in a given trophic level in a year
    base.tl.pool.tmp$prop.pool.space <- base.tl.pool.tmp$bm.tl/tl.bm.last$bm.tl
    # We can then adjust the stock K's by the available K space in each stock
    
    base.stock.pool.tmp$K.space <- NA
    base.stock.pool.tmp$K.space[base.stock.pool.tmp$troph.cat =="Low"] <- base.stock.pool.tmp$bm.stock[base.stock.pool.tmp$troph.cat =="Low"] * 
                                                        (base.tl.pool.tmp$prop.pool.space[base.tl.pool.tmp$troph.cat =="Low"]-1)
    base.stock.pool.tmp$K.space[base.stock.pool.tmp$troph.cat =="Medium"] <- base.stock.pool.tmp$bm.stock[base.stock.pool.tmp$troph.cat =="Medium"] * 
      (base.tl.pool.tmp$prop.pool.space[base.tl.pool.tmp$troph.cat =="Medium"]-1)
    base.stock.pool.tmp$K.space[base.stock.pool.tmp$troph.cat =="High"] <- base.stock.pool.tmp$bm.stock[base.stock.pool.tmp$troph.cat =="High"] * 
      (base.tl.pool.tmp$prop.pool.space[base.tl.pool.tmp$troph.cat =="High"]-1)
    
    base.stock.pool.tmp$adj.pool <- base.stock.pool.tmp$bm.stock + base.stock.pool.tmp$K.space
    
    # So now I have Carrying Capacities that take up (or lose) any available K space.
    # Now we can convert these to numbers using the historic 'average weight' of the stocks, to avoid complication
    # I'm just using the average of the average weight for each stock...
    base.stock.pool.tmp <- left_join(base.stock.pool.tmp,av.wgt,by=c("Stock","troph.cat"))
    # And now we can get a K in numbers....
    base.stock.pool.tmp$adj.pool.num <- base.stock.pool.tmp$adj.pool/base.stock.pool.tmp$mn.wgt
    # Since I have Years and sim recorded, I should just be able to recursivly rbind this...
    # if(t ==1 & j == 1) 
    # {
    #   base.stock.pool <- base.stock.pool.tmp
    #   base.tl.pool <- base.tl.pool.tmp
    #   base.com.pool <- base.com.pool.tmp
    # } else {
    #         base.stock.pool <- rbind(base.stock.pool,base.stock.pool.tmp)
    #         base.tl.pool <- rbind(base.tl.pool,base.tl.pool.tmp)
    #         base.com.pool <- rbind(base.com.pool,base.com.pool.tmp)
    #         } # end the else...
  
    #browser()
  for(s in stock.eco)
  {
      # Reset samples
      pred.mod <- mod.pred[[s]]
      bm.ts.stock <- bm.best[bm.best$Stock == s,]
      tmp.bm.last <- stock.bm.last |> collapse::fsubset(Stock == s)
      tmp.stock.pool <- base.stock.pool.tmp |> collapse::fsubset(Stock == s)
      #bm.ts.stock <- bm.final[bm.final$Stock == s & bm.final$Year %in% first.year:last.year,]  
      tl <- unique(bm.ts.stock$troph.cat)
      #a=s
      #browser()
      # Now get the final year bm
      if(t == 1) 
      { 
        #browser()
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
      
      # Sort out which of the years are low or high bm
      # I'm using 0.5 as the cut off, other options are valid (0.4 is my fav...)
       l.v.h <- 0.4
      #if(s == "ICES-HAWG_ NS-IV 3a,7d_Clupea_harengus")   l.v.h <- 0.6 # DK Note, using 0.6 for herring stock didn't decline below 50% in this time period.
      #if(s == "ICES-WGNSSK_NS4 _Scopthalmus_maximus")   l.v.h <- 0.9 # # DK Note, using 0.9 for this stock because it only declined to 66% of max in time period.
      #if(s == "ICES-HAWG_NS_Ammodytes_tobianus")  l.v.h <- 0.6 # DK Note, trying to make dynamics more realistic
       cur.pool <- tmp.stock.pool$adj.pool
       
       er <- proj.catch.eqn(dat = ex.dat[[s]],bm = bm.start)
       
       
       res <- pop.dam(stock.dat = bm.ts.stock,stock.pool = cur.pool,bm.start = bm.start, catch = er, stock = s,
                      low.vs.high = l.v.h,method=method,mod.preds = pred.mod)
          
      
      results[[s]]$bm[results[[s]]$Years==t] <- res$tst.res
      results[[s]]$removals[results[[s]]$Years==t-1] <- res$removals
      results[[s]]$ex.rate[results[[s]]$Years==t-1] <- res$ex.rate
      results[[s]]$lambda[results[[s]]$Years==t-1] <- res$lambda
      results[[s]]$K.bm[results[[s]]$Years==t-1] <- cur.pool
       

#if(tst$r$r[1] > 14) stop("WTF")
#abund.new[[s]] <- data.frame(abund = tst$Pop$abund[2])
#res.r[[s]] <- data.frame(tst$r[1,-2],stock=s,sim=j)

  } # end stock loop
    #res.ts[[t]] <- do.call("rbind",results)
    
  } # end the t looping through each year.
  res.ts[[j]] <- do.call("rbind",results)
#ggplot(base.stock.pool) + geom_line(aes(x= Years,y=bm.stock,group=Stock,color=Stock)) + facet_wrap(~troph.cat) + scale_y_log10()
  
  # Unpack the results
  #ts.unpack[[j]] <- do.call('rbind',res.ts)
  
#ggplot(ts.unpack[[j]]) + geom_line(aes(x= Years,y=abund,group=Stock,color=Stock)) + facet_wrap(~troph.cat) + scale_y_log10()
  
  
  # Pop a note when done each simulation
  timer <- Sys.time() - st.time
  print(paste("Simulation ", j))
  print(signif(timer,digits=2))
  
} # end n.sims
#browser()
# Unpack all the results.
ts.final <- do.call("rbind",res.ts)
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
saveRDS(object = ts.final,file = paste0(repo.loc,"/Results/NI/NS_projections_",n.sims,"_sims_",min(years),"_to_",max(years),
                                        "_time_series_projections.Rds"))

saveRDS(object = quants,file = paste0(repo.loc,"/Results/NI/NS_projections_",n.sims,"_sims_",min(years),"_to_",max(years),
                                      "_time_series_quantiles.Rds"))
# 
# saveRDS(object = r.final,file = paste0(repo.loc,"/Results/NI/NS_projections_",n.sims,"_sims_",min(years),"_to_",max(years),
#                                        "_r_projections.Rds"))



# Two simple plots. 
p.sims <- ggplot(ts.final) + geom_line(aes(x=Years+max(bm.best$Year),y=bm/1e3,group = sim,color=sim),alpha=0.8) +
  facet_wrap(~interaction(Stock.short, troph.cat , sep = " TC: ",drop=T),scales='free_y') +
  geom_line(data=bm.best,aes(x=Year,y=bm.stock/1e3)) +
  scale_x_continuous(name='',breaks = few.breaks) +
  scale_y_continuous(name = "Biomass (thousands of tonnes)") + 
  theme(legend.position = 'none') 

save_plot(paste0(repo.loc,"/Figures/NI/raw_biomass_trends.png"),p.sims,base_height = 12,base_width = 24)

# How about a removals and exploitation rate time series.
#browser()
p.removals <- ggplot(ts.final) + geom_line(aes(x=Years+max(bm.best$Year),y=removals/1000,group = sim,color=sim),alpha=0.8) +
  facet_wrap(~interaction(Stock.short, troph.cat , sep = " TC: ",drop=T),scales='free') +
  geom_line(data=bm.best,aes(x=Year,y=catch/1000)) +
  scale_x_continuous(name='',breaks = few.breaks) +
  scale_y_continuous(name = "Removals (thousands of tonnes)") + 
  theme(legend.position = 'none') 

save_plot(paste0(repo.loc,"/Figures/NI/raw_removals_trends.png"),p.removals,base_height = 12,base_width = 24)

p.er <- ggplot(ts.final) + geom_line(aes(x=Years+max(bm.best$Year),y=ex.rate,group = sim,color=sim),alpha=0.8) +
  facet_wrap(~interaction(Stock.short, troph.cat , sep = " TC: ",drop=T),scales='free') +
  geom_line(data=bm.best,aes(x=Year,y=er)) +
  scale_x_continuous(name='',breaks = few.breaks) +
  scale_y_continuous(name = "Exploitation rate (proportional)") + 
  theme(legend.position = 'none') 

save_plot(paste0(repo.loc,"/Figures/NI/raw_exploitation_rates.png"),p.er,base_height = 12,base_width = 24)


# 
p.sims.quants <- ggplot(quants) + geom_line(aes(x=Years+max(bm.best$Year),y=med/1e3,group=Stock,color=spec.tl),linewidth=2) +  
  #geom_line(data=bm.best,aes(x=Year,y=bm.stock,group=Stock,color=spec.tl)) +
  facet_wrap(~troph.cat) +  scale_y_log10(name="Biomass (thousands of tonnes)") +   theme(legend.position = 'top') +
  guides(colour = guide_legend(nrow = 5)) + scale_color_manual(values=pal) +
  scale_x_continuous(name='',breaks = few.breaks)
#geom_ribbon(data=quants, aes(x=Years,ymax=U.50,ymin = L.50),alpha=0.5,fill='blue',color='blue') 
save_plot(paste0(repo.loc,"/Figures/NI/Quantile_biomass_trends.png"),p.sims.quants,base_height = 8,base_width = 16)

p.sims.quants.by.stock <- ggplot(quants) + geom_line(aes(x=Years+max(bm.best$Year),y=med/1e3),linewidth=1.25,color='grey50') + 
  geom_ribbon(aes(x=Years+max(bm.best$Year),ymax=U.50/1e3,ymin=L.50/1e3),fill='firebrick2',alpha=0.2)+
  geom_line(data=bm.best,aes(x=Year,y=bm.stock/1e3),color='black')+#,linetype = 'dashed') +
  facet_wrap(~interaction(Stock.short, troph.cat , sep = " TC: ",drop=T),scales='free_y') +
  scale_y_continuous(name="Biomass (thousands of tonnes)") +   
  theme(legend.position = 'top') +
  guides(colour = guide_legend(nrow = 5)) +
  scale_x_continuous(name='',breaks = few.breaks)
#geom_ribbon(data=quants, aes(x=Years,ymax=U.50,ymin = L.50),alpha=0.5,fill='blue',color='blue') 
save_plot(paste0(repo.loc,"/Figures/NI/Quantile_biomass_trends_by_stock.png"),p.sims.quants.by.stock,base_height = 12,base_width = 24)

p.sims.rems.by.stock <- ggplot(quants) + geom_line(aes(x=Years+max(bm.best$Year),y=rem.med/1e3),linewidth=1.25,color='grey50') + 
  geom_ribbon(aes(x=Years+max(bm.best$Year),ymax=rem.U.50/1e3,ymin=rem.L.50/1e3),fill='firebrick2',alpha=0.2)+
  geom_line(data=bm.best,aes(x=Year,y=catch/1e3)) +
  facet_wrap(~interaction(Stock.short, troph.cat , sep = " TC: ",drop=T),scales='free_y') +
  scale_y_continuous(name="Removals (thousands of tonnes)") +   
  theme(legend.position = 'top') +
  guides(colour = guide_legend(nrow = 5)) +
  scale_x_continuous(name='',breaks = few.breaks)
#geom_ribbon(data=quants, aes(x=Years,ymax=U.50,ymin = L.50),alpha=0.5,fill='blue',color='blue') 
save_plot(paste0(repo.loc,"/Figures/NI/Quantile_removals_by_stock.png"),p.sims.rems.by.stock,base_height = 12,base_width = 24)
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
save_plot(paste0(repo.loc,"/Figures/NI/Quantile_er_by_stock.png"),p.sims.er.by.stock,base_height = 12,base_width = 24)


# Clean up the names of the biomass pools to distinguish them from the actual biomasses

names(sim.pool.stocks) <- c("Years","sim","Stock","troph.cat","prop.tl.pool", "cor.prop.tl.pool", "pool.init.stock", "trophic","Species","color",   
                            "spec.tl", "Stock.short","common")
names(sim.troph.pool) <- c("Years","sim","pool.tl","troph.cat","pool.com","prop.com.pool")
names(sim.com.pool) <- c("com.pool","Years","sim")

# ggplot(bm.best) + geom_line(aes(x=Year,y=prop.bm.stock.tl,group = Stock,color=spec.tl),linewidth=2) + 
#   facet_wrap(~troph.cat) + guides(colour = guide_legend(nrow = 5)) + theme(legend.position = 'top') +
#   scale_y_log10(name= "Proportion of biomass",n.breaks=10) + scale_x_continuous(name="",labels = c(1990,2000,2010),breaks=c(1990,2000,2010))+
#   scale_color_manual(values=pal)

return(list(sim.quantiles = quants,
            sim.ts = ts.final,
            past.bm = bm.best,
            sim.pool.stocks = sim.pool.stocks,
            sim.troph.pool = sim.troph.pool,
            sim.com.pool = sim.com.pool))
} # end function