# Here we develop a multi-species model for the North Sea.

#################  Section 1 Loading #################  Section 1 Loading #################  Section 1 Loading  ###############################################
library(tidyverse)
library(GGally)
library(cowplot)
library(ggthemes)
library(rfishbase)
library(boot)
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

# Want the local version for now...
source("D:/Github/ICM/Scripts/functions/simple_forward_sim.r")


#load(file = "D:/Github/ICM/Results/model_inputs.Rdata")
dat.loc <- 'D:/GitHub/ICM'
repo.loc <- "D:/GitHub/Multispecies_model/"
#loc <- "C:/Users/Owner/Documents/Github/ICM"
load(file = paste0(dat.loc,"/Results/all_cleaned_forward_tune_summaries_no_age_corection_fec_nm.Rdata"))
load(file = paste0(dat.loc,"/Results/model_inputs_no_age_correction.Rdata"))

########################### End Section 1 Loading  ###################################### End Section 1 Loading ###############################################

########################## Section 2 Parameters ########################## Section 2 Parameters ########################## Section 2 Parameters

n.yrs.proj <- 100 # How many years into the future we are going to project the stocks
n.sims <- 500 # The numbers of simulations to run, keeping low for testing...

# We've spelt the name of Turbot wrong and the name of arrow-tooth flounder. Tidy up of other names because of spaces and capitialization.
ASR.long$Genus[ASR.long$Genus == 'Scopthalmus'] <- "Scophthalmus" 
ASR.long$Genus[ASR.long$Genus == 'clupea'] <- "Clupea" 
ASR.long$Genus[ASR.long$Genus == ' Hippoglossoides'] <- "Hippoglossoides" 
ASR.long$Genus[ASR.long$Genus == 'Hippoglossoides '] <- "Hippoglossoides" 
ASR.long$Genus[ASR.long$Genus == 'Scomber '] <- "Scomber" 
ASR.long$Genus[ASR.long$Genus == 'Dicentrarchus '] <- "Dicentrarchus" 
ASR.long$Genus[ASR.long$Genus == 'Pollachius '] <- "Pollachius" 
ASR.long$Genus[ASR.long$Genus == 'Sardina '] <- "Sardina" 
ASR.long$Species[ASR.long$Species == 'Aeglefinus'] <- "aeglefinus" 
ASR.long$Species[ASR.long$Species == 'Chrysops'] <- "chrysops" 
ASR.long$Species[ASR.long$Species == ' harengus'] <- "harengus" 
ASR.long$Species[ASR.long$Species == 'stomais'] <- "stomias" 
ASR.long$Species[ASR.long$Species == 'Solea'] <- "solea" 
ASR.long$Gen.Spec <- paste(ASR.long$Genus,ASR.long$Species,sep=' ')


# Get the right stocks
Stocks <- names(for.tune.all)
ns.stocks <- Stocks[grep("NS",Stocks)]
ns.stocks <- ns.stocks[ns.stocks != "ICES-WGHANSA_SP8abd_Sardina _pilchardus"]

# Get the trophic levels from Fishbase

gen.spec <- unique(paste(ASR.long$Genus,ASR.long$Species,sep=" "))


# In fishbase "Limanda ferruginea" is Myzopsetta ferruginea, so add that, then change the name back to what we've been using.  Seems both are in use but Limanda feels more common.
# here's the trophic levels from Fishbase, looks like we have 2 trophic levels with 5 species in each
all.troph <- ecology(c(gen.spec,"Myzopsetta ferruginea"), fields=c("Species","SpecCode", "FoodTroph", "FoodSeTroph", "DietTroph", "DietSeTroph"))
all.troph$Species[all.troph$Species == "Myzopsetta ferruginea"] <- "Limanda ferruginea"

# A few species pop up multple times for some reason
mutli.entries <- names(which(table(all.troph$Species) > 1))
# Clean up the multiple entries
if(length(mutli.entries) > 0)
{
  drop <- NULL
  for(i in 1:length(mutli.entries))
    drop[[i]] <- which(all.troph$Species == mutli.entries[i])[-1] # Drop all but the first one
}
drop <- do.call('c',drop)

all.troph <- all.troph[-drop,]
# ns.troph <- all.troph |> collapse::fsubset(Species %in% gen.spec)
# Add the trophic level information to ASR.long
ASR.long <- left_join(ASR.long,all.troph,by = join_by(Gen.Spec==Species))

# Trophic levels from Simon Jennings Paper in 2002 for NS...or fishbase if not in Jennings (e.g. lesser sand eel and Turbot)
ns.troph <- data.frame(Stock = ns.stocks,
                       Common = c("Herring","Lesser Sand eel","Sole","Atlantic cod",
                                  "Haddock","European plaice","Norway pout","Saithe",
                                  "Atlantic cod", "Whiting", "Sole", "European plaice","Turbot","Sole"),
                       TL = c(3.8,3.08,5.0,5.2,
                              4.7,4.5,4.2,4.6,
                              5.2,5.3,5.0,4.5,4.4,5.0))

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
for(s in  ns.stocks)
{
  years.ns[[s]] <- years.tmp[[s]]
  #vpa.ns[[s]] <- vpa.tmp[[s]]
  ages.ns[[s]] <- ages.tmp[[s]]

  num.ns[[s]] <- ASR.long |> collapse::fsubset(Stock == s & type == "Num")
  num.ns[[s]] <- num.ns[[s]] |> collapse::fsubset(age != "tot")
  waa.ns[[s]] <- ASR.long|> collapse::fsubset(Stock == s & type == "WA")
  tl <- ns.troph |> collapse::fsubset(Stock == s)
  if(s == "ICES-HAWG_NS_Ammodytes_tobianus") waa.ns[[s]]$value <- waa.ns[[s]]$value/1000
  bm.ns[[s]] <- data.frame(Year = num.ns[[s]]$Year,Stock = num.ns[[s]]$Stock,age = num.ns[[s]]$age,
                           bm = num.ns[[s]]$value*waa.ns[[s]]$value,
                           num = num.ns[[s]]$value,
                           trophic = tl$TL,
                           troph.cat = as.character(floor(tl$TL)),
                           Species = num.ns[[s]]$Gen.Spec)
  pnm.ns[[s]] <- 1-exp(-for.tune.all[[s]]$nm.opt)
  mx.ns[[s]] <- for.tune.all[[s]]$fecund.opt
  vpa.ns[[s]] <- for.tune.all[[s]]$res$est.abund
  rem.ns[[s]] <- rem.tmp[[s]]
  rem.ns[[s]]$Stock <- s
  am.ns[[s]] <- am.tmp[[s]]
} # end for(s in  ns.stocks)
# Combine the biomass and abundance data into a dataframe
bm.tst <- do.call("rbind",bm.ns)


# Look at the biomass and abundance in the ecosystem
bm.tot <- bm.tst |> collapse::fgroup_by(Stock,Year,trophic,Species,troph.cat) |> 
                    collapse::fsummarize(bm = sum(bm,na.rm=T),
                                         num = sum(num,na.rm=T))

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

# Biomass by trophic level over time
ggplot(bm.best) + geom_line(aes(x=Year,y=bm.tl,group=troph.cat,color=troph.cat)) + 
  scale_color_manual(values = c("blue","red","darkgrey","lightgreen")) + scale_y_log10(name="Biomass")

# This is real good now...
ggplot(bm.best) + geom_line(aes(x=Year,y=prop.bm.tl,group=troph.cat,color=troph.cat)) + 
  scale_color_manual(values = c("blue","red","darkgrey","lightgreen")) + scale_y_continuous(name="Proportion of Biomass")
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

colors <- distinct(bm.best, spec.tl, color)
pal <- colors$color
names(pal) <- colors$spec.tl

ggplot(bm.best) + geom_line(aes(x=Year,y=prop.bm.stock.tl,group = Stock,color=spec.tl),linewidth=2) + 
                  facet_wrap(~troph.cat) +
                  scale_y_log10(n.breaks=10) + scale_color_manual(values=pal)

ggplot(bm.best) + geom_line(aes(x=Year,y=bm.stock,group = Stock,color=species)) + 
                  facet_wrap(~troph.cat) +
                  scale_y_log10(n.breaks=7)


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
K.cor <- pacf(eco.tot.bm.best$bm.eco)
# The cross correlation between the ecosystem biomass trend and the trophic level biomasses
# All correlated, but strongest is unsurprisingly the link between the the ecosystem and the biomass in the
# lowest TL. I suspect this may structurally come out even without explicity building in a lot of
# correlation structure to the models.
K.tl.3.cor <- ccf(eco.tot.bm.best$bm.eco,trophic.bm.best$bm.tl[trophic.bm.best$troph.cat==3])
K.tl.4.cor <- ccf(eco.tot.bm.best$bm.eco,trophic.bm.best$bm.tl[trophic.bm.best$troph.cat==4])
K.tl.5.cor <- ccf(eco.tot.bm.best$bm.eco,trophic.bm.best$bm.tl[trophic.bm.best$troph.cat==5])
# Within trophic levels...
# So these 3 mostly say if the biomass is up one TL, it is up in all TLs, tho there might be some negative between 3 and 4
# at Lag -1 (though that's not quite significant)
tl.3.4.cor <- ccf(trophic.bm.best$bm.tl[trophic.bm.best$troph.cat == 3],trophic.bm.best$bm.tl[trophic.bm.best$troph.cat == 4])
tl.3.5.cor <- ccf(trophic.bm.best$bm.tl[trophic.bm.best$troph.cat == 3],trophic.bm.best$bm.tl[trophic.bm.best$troph.cat == 5])
tl.4.5.cor <- ccf(trophic.bm.best$bm.tl[trophic.bm.best$troph.cat == 4],trophic.bm.best$bm.tl[trophic.bm.best$troph.cat == 5])
# Looking at proportions, need to stew a bit on this because there is necessarily some 
# correlation built into proportions, but what is interesting is that
# the correlation strength is really really high between TL 3 and TL 4, it is weaker (more diffuse really) at TL 5
# and there is no correlation between 4 and 5
tl.3.4.prop.cor <- ccf(bm.best$prop.bm.tl[bm.best$troph.cat == 3][1:n.years],bm.best$prop.bm.tl[bm.best$troph.cat == 4][1:n.years])
tl.3.5.prop.cor <- ccf(bm.best$prop.bm.tl[bm.best$troph.cat == 3][1:n.years],bm.best$prop.bm.tl[bm.best$troph.cat == 5][1:n.years])
tl.4.5.prop.cor <- ccf(bm.best$prop.bm.tl[bm.best$troph.cat == 4][1:n.years],bm.best$prop.bm.tl[bm.best$troph.cat == 5][1:n.years])


# So now really what I need to do first is make a quick simulation that gets me ecosystem K, trophic level K, and stock K
# once I have those then we just run the models :-)
mn.eco.bm <- mean(eco.tot.bm.best$bm.eco)
start.eco.sim <- eco.tot.bm.best$bm.eco[length(eco.tot.bm.best$bm.eco)]
sd.eco.bm <- sd(eco.tot.bm.best$bm.eco)
# trophic level biomass and proportions... for the proportion will probably wanna sample from a beta distro
# So not sure how to do that nicely...
mn.tl3.bm <- mean(trophic.bm.best$bm.tl[trophic.bm.best$troph.cat ==3])
sd.tl3.bm <- sd(trophic.bm.best$bm.tl[trophic.bm.best$troph.cat ==3])
mn.tl3.prop.bm <- mean(bm.best$prop.bm.tl[bm.best$troph.cat == 3][1:n.years])
sd.tl3.prop.bm <- sd(bm.best$prop.bm.tl[bm.best$troph.cat == 3][1:n.years])
mn.tl4.bm <- mean(trophic.bm.best$bm.tl[trophic.bm.best$troph.cat ==4])
sd.tl4.bm <- sd(trophic.bm.best$bm.tl[trophic.bm.best$troph.cat ==4])
mn.tl4.prop.bm <- mean(bm.best$prop.bm.tl[bm.best$troph.cat == 4][1:n.years])
sd.tl4.prop.bm <- sd(bm.best$prop.bm.tl[bm.best$troph.cat == 4][1:n.years])
mn.tl5.bm <- mean(trophic.bm.best$bm.tl[trophic.bm.best$troph.cat ==5])
sd.tl5.bm <- sd(trophic.bm.best$bm.tl[trophic.bm.best$troph.cat ==5])
mn.tl5.prop.bm <- mean(bm.best$prop.bm.tl[bm.best$troph.cat == 5][1:n.years])
sd.tl5.prop.bm <- sd(bm.best$prop.bm.tl[bm.best$troph.cat == 5][1:n.years])

# This little puppy will take the mean/variance parameters from our data into the beta parameters
# Extract the Beta parameters from the mean and variance of your data
estBetaParams <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}
tl3.prop.bm.params <- estBetaParams(mn.tl3.prop.bm,sd.tl3.prop.bm^2)
tl4.prop.bm.params <- estBetaParams(mn.tl4.prop.bm,sd.tl4.prop.bm^2)
tl5.prop.bm.params <- estBetaParams(mn.tl5.prop.bm,sd.tl5.prop.bm^2)
# Doesn't do badly, gets the mean/spread about right, not enough data in our distro to say much else...
hist(bm.best$prop.bm.tl[bm.best$troph.cat == 3][1:n.years])
hist(rbeta(10000,tl3.prop.bm.params$alpha,tl3.prop.bm.params$beta))
hist(bm.best$prop.bm.tl[bm.best$troph.cat == 4][1:n.years])
hist(rbeta(10000,tl4.prop.bm.params$alpha,tl4.prop.bm.params$beta))
hist(bm.best$prop.bm.tl[bm.best$troph.cat == 5][1:n.years])
hist(rbeta(10000,tl5.prop.bm.params$alpha,tl5.prop.bm.params$beta))


# First get the ecosystem biomass in a correlated time series, there are a whole lot of ways one could do this, this
# is one of many different ideas. I think we could get the 4 and 5 correlations better another way, but
# For a first pass I'm ok with this.
# Ok, duh, use the mean of the time series then the arima gives us the deviations from that mean and we get a nice time series.

#ggplot(bm.trophic.Ks[[i]]) + geom_line(aes(x=Years,y=prop.bm.tl,group = troph.cat,color=troph.cat))

# Used for simulations to get good time series for the K for TL3,4, and 5 
tl.3.prop.bm.ts <- bm.best$prop.bm.tl[bm.best$troph.cat==3][1:n.years]
# Extract the frist two components from the pacf to get the two AR components from the model.
tl.3.prop.pacf <- pacf(tl.3.prop.bm.ts,plot = F)
tl.3.prop.bm.lag.1 <- tl.3.prop.pacf$acf[1]
tl.3.prop.bm.lag.2 <- tl.3.prop.pacf$acf[2]
# TL 4 and 5 splits historically
tl.4.5.prop.bm <- bm.best$bm.tl[bm.best$troph.cat == 5][1:n.years]/(bm.best$bm.tl[bm.best$troph.cat == 4][1:n.years]+bm.best$bm.tl[bm.best$troph.cat == 5][1:n.years])
# This is the correlation between 4 and 5
tl.4.5.prop.4.5.bm <- pacf(tl.4.5.prop.bm)


troph.levels <- sort(unique(bm.best$troph.cat))

# 
sim.K.stock <- NULL
sim.Ks <- NULL
sim.eco.bm <- NULL
bm.trophic.Ks <- NULL
# Starting values for the ecosystem and the proportions
start.eco.sim <- eco.tot.bm.best$bm.eco[nrow(eco.tot.bm.best)]
start.eco.diff = start.eco.sim- mn.eco.bm
# Starting and mean values for the trophic levels
mn.tl.3.prop.bm <- mean(tl.3.prop.bm.ts)
start.tl.3.prop.bm <- tl.3.prop.bm.ts[length(tl.3.prop.bm.ts)]
# convert to logit scale for the arima models
tl.3.logit <- logit(tl.3.prop.bm.ts)
start.tl.3.logit <- tl.3.logit[length(tl.3.logit)]
mn.tl.3.logit <- mean(tl.3.logit)
strat.tl.3.diff <- start.tl.3.logit - mn.tl.3.logit
sd.tl.3.logit <- sd(tl.3.logit)
# Now the proprotion between 4 and 5
mn.tl.4.5.prop.bm <- mean(tl.4.5.prop.bm)
start.tl.4.5.prop.bm <- tl.4.5.prop.bm[length(tl.4.5.prop.bm)]
tl.4.5.prop.bm.lag.1 <- tl.4.5.prop.4.5.bm$acf[1]
# convert to logit scale for the arima models
tl.4.5.logit <- logit(tl.4.5.prop.bm)
start.tl.4.5.logit <- tl.4.5.logit[length(tl.4.5.logit)]
mn.tl.4.5.logit <- mean(tl.4.5.logit)
strat.tl.4.5.diff <- start.tl.4.5.logit - mn.tl.4.5.logit
sd.tl.4.5.logit <- sd(tl.4.5.logit)



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
                                           n.start =2, start.innov = c(strat.tl.3.diff/tl.3.prop.bm.lag.1,strat.tl.3.diff/tl.3.prop.bm.lag.1), 
                                           innov = c(0,rnorm(n.yrs.proj-1,0,sd.tl.3.logit))))
  
  bm.sim.3 <- sim.tl.3.prop.bm * sim.eco.bm[[i]]$bm
  # So this is what is left for 3 and 4
  bm.left.4.5<- sim.eco.bm[[i]]$bm - bm.sim.3
  # So then we use the historical split between 4 and 5 can see 5 gets about 1/3-1-5 of 3
   # so then simulate this split
  sim.tl.5.4.prop.bm <- inv.logit(mn.tl.4.5.logit + 
                                    arima.sim(model =list(ar = tl.4.5.prop.bm.lag.1),n = n.yrs.proj,
                                              n.start =1, start.innov = c(strat.tl.4.5.diff/tl.4.5.prop.bm.lag.1), 
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
                                   Stock = s, troph.cat = tl,
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
#ggplot(sim.K.stocks) + geom_line(aes(x=Years,y=bm.stock,group=Stock,color=as.factor(sim))) + facet_wrap(~troph.cat)
#ggplot(sim.troph.K) + geom_line(aes(x=Years,y=bm.tl,group=as.factor(sim),color=as.factor(sim))) + facet_wrap(~troph.cat)
#ggplot(sim.eco.K) + geom_line(aes(x=Years,y=bm,group=as.factor(sim),color=as.factor(sim))) 
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
# ggplot(sim.K.stocks |> collapse::fsubset(sim == 1)) + geom_line(aes(x=Years,y=bm.stock,group=Stock,color=troph.cat)) + scale_y_log10()
# ggplot(bm.best) + geom_line(aes(x=Year,y=bm.stock,group=Stock,color=troph.cat))+ scale_y_log10()



# Fix: This is not perfect way to get the past exploitation rates as the removals we have here are in numbers
# Give we have the database with the age specific removals and age specific weights, this
# should be tweaked to use that data. That said, for the moment this should be 'good enough' 
rem.tst <- do.call("rbind",rem.ns)
fm.dat <- left_join(bm.best,rem.tst,by=c("Stock",'Year'))
# This is where we go from numbers to a biomass and get an exploitation rate in biomass.
fm.dat$exploit <- (fm.dat$rem*fm.dat$avg.weight)/fm.dat$bm.stock

ns.mod.fit <- for.tunes |> collapse::fsubset(Stock %in% ns.stocks) 
ns.mod.fit <- ns.mod.fit |> collapse::fgroup_by(Stock) |> collapse::fmutate(max.num = max(est.abund,na.rm=T)) |> as.data.frame()
ns.mod.fit$prop.max <- ns.mod.fit$est.abund/ns.mod.fit$max.num
#ggplot(ns.mod.fit,aes(x=prop.max,y=log(mean.fec),group=Stock,color=Stock)) + geom_point() + geom_smooth(method = 'lm')
#ggplot(ns.mod.fit,aes(x=prop.max,y=log(mean.nm),group=Stock,color=Stock)) + geom_point() + geom_smooth(method = 'lm')
#ggplot(ns.mod.fit,aes(x=prop.max,y=log(lambda),group=Stock,color=Stock)) + geom_point() + geom_smooth(method = 'lm')


############### Section 4 Multi-species model of North Sea ############### Section 4 Multi-species model of North Sea ############### Section 4 Multi-species model of North Sea

# Now we have our carrying capacity for each stock and we can get to business and running a model.
# Initialize some things, or maybe no things
res.ts <- NULL
ts.unpack <- NULL
# Get the year range, going from the 'last' year to n.yrs.proj in the future, note this will go 1 year less than your intuition because
# we want n.yrs of data, i.e., 20 years is 2000 to 2019, not 2020... )
years <- (last.year+1):(last.year+n.yrs.proj)
# Take the biomass data for the north sea and subset it to the years we have data
bm.mod.yrs <- bm.tst |> collapse::fsubset(Year %in% first.year:last.year)
bm.2014 <- bm.mod.yrs |> collapse::fsubset(Year == last.year) |> 
                         collapse::fgroup_by(Stock,trophic,troph.cat,Species) |> 
                         collapse::fsummarise(bm.tot = sum(bm,na.rm=T))
# Get the initial ecosystem biomass..
init.eco.bm <- sum(bm.2014$bm.tot)
init.tl.bm <- bm.2014 |> collapse::fgroup_by(troph.cat) |> collapse::fsummarise(bm.tl = sum(bm.tot))
init.stock.bm <- bm.2014
# Get the average weight of the fish in the stocks so we can go from biomass to abundance for the model
# FIX: This could definitely be done more sophisisticatedly!
av.wgt <- fm.dat |> collapse::fgroup_by(Stock,troph.cat) |> collapse::fsummarise(mn.wgt = mean(avg.weight,na.rm=T))


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
      bm.stocks <- data.frame(abund = NA, bm = NA,av.wgt)
      for(s in ns.stocks) bm.stocks$abund[bm.stocks$Stock == s] <- res.ts[[s]]$abund[res.ts[[s]]$Years == t-1]
      bm.stocks$bm <- bm.stocks$abund*bm.stocks$mn.wgt
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

  for(s in ns.stocks)
  {
      # Reset samples
      mx.samp <- NA
      nm.samp <- NA
      stock.fit <- ns.mod.fit |> collapse::fsubset(Stock == s)
      tmp.bm.last <- stock.bm.last |> collapse::fsubset(Stock == s)
      tmp.stock.K <- base.stock.K.tmp |> collapse::fsubset(Stock == s)
      fm.stock <- fm.dat |> collapse::fsubset(Stock ==s) |> collapse::fsummarize(mn = median(exploit,na.rm=T),
                                                                                      sd = sd(log(exploit[exploit > 0]),na.rm=T))
      
      
      # Now get the final year abundance
      if(t == 1) 
      {
        # The number of individuals (using our retrospective simulations)
        vpa.ns  <- bm.final$num.stock[bm.final$Stock == s]
        N.start <- vpa.ns[length(vpa.ns)]
        res.ts[[s]] <- data.frame(abund = N.start,removals = NA,Stock = s,sim= j,r = NA,Years=t-1,
                                  troph.cat = floor(ns.troph$TL[ns.troph$Stock ==s]),
                                  K.num = NA)
        # Get the ages here too...
        ages <- ages.ns[[s]]
      } else{  N.start <- res.ts[[s]]$abund[res.ts[[s]]$Years == t-1]}
      
      # Sort out which of the years are low or high abundance
      # I'm using 0.5 as the cut off, other options are valid (0.4 is my fav...)
      low.vs.high <- 0.5
      low.abund.years <- which(stock.fit$est.abund < low.vs.high*stock.fit$max.num[1])
      if(length(low.abund.years) == 0) low.years <- F else low.years <- T
      high.abund.years <- which(stock.fit$est.abund >= low.vs.high*stock.fit$max.num[1])
        
      
      # FIX: THERE ARE MANY MANY WAYS WE COULD DO THIS AND IT IS SUPER IMPORTANT TO THE SIMULATIONS, SO EPLOXRE OTHER OPTIONS.
      # FIX, my solution here is to adjust the mx and lx vectors that are sampled using the density dependence relationship
      # and then the distance to the stock is from the carrying capacity. 
      #Obtain a sample for the fecundity and natural mortality, we are doing this by splitting the time series based on the
      # current stock status and the available K space.  This piece is a most interesting part of the simulations and is
      # really important. Based on our LTR paper I'm going to split this into 3 components
      # If the stock is 
      # 1) above the available K, we take the minimum fecundity observed, put some variability on that, and get a 
      # fecundity estimate for the stock
      # 2) when biomass was > 50% of K, we sample the fecundity observed when the stock was > 50% of the maximum observed fecundity
      # 3) below 50% of K we sample the fecundity when the stock was < 50% of max biomass
      # For the moment I ignore the potential for DD in natural mortality, but following similar logic we could do the same
      # thing for the natural mortality term. One reason I hesitate on that is that the natural mortalities often
      # are static both temporally and across age classes, another reason is I'm lazy.
      # So first, get a sample from the natural moratlites
      prop.nat.mort <- pnm.ns[[s]] 
      method <- "not_sample"
      if(method == 'sample')
      {
        # Pick one of these to sample if that's how we want to roll
        samp <- sample(nrow(-(log(1-prop.nat.mort))),1)
        # The simple way to do it is just to sample from the natural mortality distribution
        nm.samp <- -(log(1-prop.nat.mort))[samp,] # Get the sample years.  
      }
      
      # Or do it the fun way...
      if(method != "sample")
      {
        # The fun way to do it is to do something multivariate! Note these are instantaneous now!!
        nm.mns <- colMeans(-(log(1-prop.nat.mort)))
        sd.mns <- apply(-(log(1-prop.nat.mort)),2,sd)
        nm.covar <- cov(-(log(1-prop.nat.mort)))
        # I may need this to be the correlatoin not covariance!! Also need to deal with the potential for negative values 
        # This should also be done on the instantaneous scale, not the proportional, thus what we have here is instantaneous NM!!
        nm.samp <- mvtnorm::rmvnorm(n = 1, mean = nm.mns, sigma = nm.covar)
      }

      
      # Now we can do the same thing with the fecundity, but accounting for the stock status
      # FIX: For the moment I haven't accounted for the age of the recruits, not too complicated to do, but
      # It should be done in theory!
      ratioed <- N.start/tmp.stock.K$adj.K.num
      mx <- mx.ns[[s]] 
      # What rows are NA's?
      na.rows <- which(is.na(mx[,1]))
      if(length(na.rows) > 0)
      {
        mx <- mx[-na.rows,] # removal any rows that are all NAs
        # I also have to remove these from the high/low data above...
        low.abund.years <- low.abund.years[low.abund.years != na.rows]
        high.abund.years <- high.abund.years[high.abund.years != na.rows]  
      }
      # Also drop the immature ages for a sec
      immature <- which(colMeans(mx) == 0)
      if(length(immature) > 0) mx <- mx[,-immature]
      # Now get the year with the lowest mx...
      mx.mns.year <- rowMeans(mx)
      min.mx.year <- which(mx.mns.year == min(mx.mns.year))
      fec.covar <- cov(log(mx))
      min.fec <- c(unlist(mx[min.mx.year,]))
      # And get the year with the high natural mortality
      nm.years <- rowMeans(-(log(1-prop.nat.mort)))
      max.nm.year <- which(nm.years == max(nm.years))
      max.nm <- c(unlist(-(log(1-prop.nat.mort))[max.nm.year,]))
      
      # If there are still some 0 years mixed in here, we'll do this, just for the covariance matrix to work.
      if(any(mx == 0)) 
      {
        the.zeros <- data.frame(which(mx==0, arr.ind=TRUE))
        mx.replace <- 0.5*min(unlist(mx)[which(unlist(mx) > 0)])
        mx.cov.tmp <- mx
        for(z in 1:nrow(the.zeros)) mx.cov.tmp[the.zeros$row[z],the.zeros$col[z]] <- mx.replace
        fec.covar <- cov(log(mx.cov.tmp))
        min.fec <- c(unlist(mx.cov.tmp[min.mx.year,]))
      }
      
      # This works, note that I've made the covariance matrix values be 10% of the whole thing, so that keeps the fecundity low there.
      if(ratioed >= 1) 
      {
        mx.samp <-  exp(mvtnorm::rmvnorm(n = 1, mean = log(min.fec), sigma = 0.1*fec.covar,rnorm=stats::rnorm))
        nm.samp <-  exp(mvtnorm::rmvnorm(n = 1, mean = log(max.nm), sigma = 0.1*nm.covar,rnorm=stats::rnorm))
      }
      
      if(ratioed >= low.vs.high & ratioed < 1)
      {
        mn.fec <- colMeans(mx[high.abund.years,])
        if(any(mx == 0))  fec.covar <- cov(log(mx.cov.tmp[high.abund.years,])) else { fec.covar <- cov(log(mx[high.abund.years,]))}
        mx.samp <-  exp(mvtnorm::rmvnorm(n = 1, mean = log(mn.fec), sigma = fec.covar,rnorm=stats::rnorm))
      }
      
      # FIX: This isn't able to get the really high recruitment events from what I'm seeing
      # so worth thinking about how to account for those events
      
      
      # PROBLEM!  SOME YEARS HAVE FECUNDITY OF 0 in first AGE CLASS FOR TOBIAS!!
      
      if(ratioed < low.vs.high & low.years ==T)
      {
        mn.fec <- colMeans(mx[low.abund.years,])
        # For the covariance matrix, if you have some 0's, make them 50% of the minimum observed for that row, only doing so 
        # for the covariance calcs since they are needed on log scale...
        if(any(mx == 0))  fec.covar <- cov(log(mx.cov.tmp[low.abund.years,])) else { fec.covar <- cov(log(mx[low.abund.years,]))}
        mx.samp <-  exp(mvtnorm::rmvnorm(n = 1, mean = log(mn.fec), sigma = fec.covar,rnorm=stats::rnorm))
      }
      
      # If we don't have an mx sample yet, then we get one using all the data... This should only be needed when we don't have biomasses observed below the low.vs.high value.
      if(low.years == F)
      {
        mn.fec <- colMeans(mx)
        mx.samp <-  exp(mvtnorm::rmvnorm(n = 1, mean = log(mn.fec), sigma = fec.covar,rnorm=stats::rnorm))
      }
      if(length(immature) > 0) mx.samp <- c(rep(0,length(immature)),mx.samp)
      # That in this set up this picks the fecundity and natural mortality from the 'same' year, which is 'fine', but lots of options here.
      # FIX: We could do this sample as a proper statistical sample rather than picking exactly what has been observed in the past.
    
    # Run the simulations
    tst <- simp.for.sim(years= 1,
                        nm = nm.samp, # Converted to instantaneous already
                        ages = ages,
                        rems =  list(0.01,fm.stock$sd), #fm,
                        fecund = mx.samp,
                        N.start = N.start,
                        pop.model = 'exponential', 
                        sim= "project",
                        n.sims = 1,
                        #K = Ks,
                        repo='preload')
    
# Because I'm only doing this one year at a time, there's something in here I need to mess around with to get the output tidy...
res.ts[[s]] <- rbind(res.ts[[s]] ,data.frame(abund = tst$Pop$abund[2],removals = tst$Pop$removals[1],
                                             Stock = s,sim= j,r = tst$r$r[1],Years=t,
                                             troph.cat = floor(ns.troph$TL[ns.troph$Stock ==s]),
                                             K.num = tmp.stock.K$adj.K.num))
if(tst$r$r[1] > 14) stop("WTF")
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

ggplot(ts.final) + geom_line(aes(x= Years,y=abund,group=sim,color=sim)) + facet_wrap(~Stock) + scale_y_log10()


ts.final$fm <- ts.final$removals/ts.final$abund
r.final <- do.call("rbind",r.unpack)

quants <- ts.final |>  collapse::fgroup_by(Years,Stock,troph.cat) |> collapse::fsummarize(L.50 = quantile(abund,probs=c(0.25),na.rm=T),
                                                                          med = median(abund,na.rm=T),
                                                                          U.50 = quantile(abund,probs=c(0.75),na.rm=T))#,
                                                                          #fml.50 = quantile(fm,probs=c(0.25),na.rm=T),
                                                                          #fm = median(fm,na.rm=T),
                                                                          #fmu.50 = quantile(fm,probs=c(0.75),na.rm=T))
# If happy save the 2 objects
saveRDS(object = ts.final,file = paste0(repo.loc,"/Results/NS_projections_",n.sims,"_sims_",min(years),"_to_",max(years),
                                         "_time_series_projections.Rds"))

saveRDS(object = quants,file = paste0(repo.loc,"/Results/NS_projections_",n.sims,"_sims_",min(years),"_to_",max(years),
                                       "_time_series_quantiles.Rds"))
# 
# saveRDS(object = r.final,file = paste0(repo.loc,"/Results/NS_projections_",n.sims,"_sims_",min(years),"_to_",max(years),
#                                        "_r_projections.Rds"))


# Two simple plots. 
p.sims <- ggplot(ts.final ) + geom_line(aes(x=Years,y=abund,group = sim,color=sim),alpha=0.8) +
  facet_wrap(~Stock,scales = 'free_y') + ylim(c(0,NA)) + scale_x_continuous(breaks = seq(2010,2200,by=10)) 
#save_plot(paste0("D:/Github/ICM/Figures/NS_sims/NS_all_realizations_climate_starts_at_year_",c.effect,
#                 "_mx_decade_effect_",climate.mx.effect, "_nm_decade_effect_",climate.nm.effect,".png"),p.sims,base_height = 12,base_width = 20)

p.sims.quants <- ggplot(quants) + geom_line(aes(x=Years,y=med,group=Stock,color=Stock)) + facet_wrap(~troph.cat,scales = 'free_y') + ylim(c(0,NA)) #+
  #geom_ribbon(data=quants, aes(x=Years,ymax=U.50,ymin = L.50),alpha=0.5,fill='blue',color='blue') 
#save_plot(paste0("D:/Github/ICM/Figures/NS_sims/NS_quantiles_climate_starts_at_year_",c.effect,
#                 "_mx_decade_effect_",climate.mx.effect, "_nm_decade_effect_",climate.nm.effect,".png"),p.sims.quants,base_height = 12,base_width = 20)


