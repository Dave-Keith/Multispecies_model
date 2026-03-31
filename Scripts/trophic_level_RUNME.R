# FIX For the MAESTRO work, lets simulate the ecosystem a fixed number of times, then run the population dynamics on each of these ecosystem
# simulations, instead of simulating the ecosystem for each of the population dynamics scenarios, for example, we can simulate 1,000 
# ecosystems.  Then, for the fishery dynamics, if we are testing 10 scenarios and simulating the population dynamics 100 times, that 
# all happens within a set of 1000 ecosystems scenarios (so each ecosystems would have 1000 population dynamics simulations run on it in this scenario)


# Three scenarios
# 1) NO fishing, top down vs bottom up just to show what happens in our model world.
# 2) Fishing with the following harvest control rules... our DFO ones
# Should get a metric of landings for each stock out of this too.
# For another day
# 3) Some sort of degreation of primary productiving hitting TL and and implaications of this for TD and bottom up worlds with fishing
#    as has been done historically (i.e., using the 'old' HCRs).
# 4) Using the model to recreate what was observed in the ecosystem and see if top down or bottom up works better. If this worked then we'd be talking 
#    about this as an assessment tool rather than a new way to explore your world.


library(readxl)
library(tidyverse)
library(ggplot2)
library(cowplot)
library(collapse)
n.yrs.proj <- 50 # How many years into the future we are going to project the stocks
n.sims <- 20 # The numbers of simulations to run, keeping low for testing...


#dat.loc <- 'C:/Users/keithd/Documents/GitHub/ICM'
#repo.loc <- "C:/Users/keithd/Documents/GitHub/Multispecies_model"
dat.loc <- repo.loc <- "D:/GitHub/Multispecies_model"

#loc <- "C:/Users/keithd/Documents/Github/ICM"
for.tune.summary <- readRDS(file = paste0(dat.loc,"/Results/forwards_fully_tuned_summary.Rds"))
res.lambda.final <- readRDS(file = paste0(dat.loc,"/Results/lambdas.Rds"))
load(file = paste0(dat.loc,"/Results/model_inputs.Rdata"))
# load in the trophic mocdel function

# We've spelt the name of Turbot wrong and the name of arrow-tooth flounder. Tidy up of other names because of spaces and capitialization.

ASR.long$Gen.Spec <- paste(ASR.long$Genus,ASR.long$Species,sep=' ')


# Get the right stocks, this is from our model runs
eco.stocks <- names(for.tune.summary)


# Get the trophic levels
# Trophic levels from Simon Jennings Paper in 2002 for NS...or fishbase if not in Jennings (e.g. Turbot)
# The TL for Sand eel and witch flounder is 4, but I am lumping 3.8 and 4 together as one, doing so in the code by making the TL be the floor()
# so these as 3.99 makes that work without doing anything greasy.
eco.troph <- read_xlsx( paste0(dat.loc,"/Data/species_TL.xlsx"))
# This gets some of the stock data we need.
stocks <- ASR.long |> collapse::fsubset(Stock %in% eco.stocks)
# Get the trophic levels for the stocks
stocks <- merge(stocks,eco.troph,by=c("Stock"))
# Now make a short stock ID
stocks$Stock.short <- paste(stocks$S.short,stocks$area)

bm <- NULL
num <- NULL
waa <- NULL
fm <- NULL
lambs <- NULL
rem.age <- NULL
for(s in  eco.stocks)
{
  tmp <- stocks |> collapse::fsubset(Stock == s)
  #years[[s]] <- unique(tmp$year)
  #vpa[[s]] <- vpa.tmp[[s]]
  #lambs[[s]] <- res.lambda.final[res.lambda.final$stock==s,]
  num[[s]] <- tmp |> collapse::fsubset(type == "Num")
  num[[s]] <- num[[s]] |> collapse::fsubset(age != "tot")
  waa[[s]] <- tmp |> collapse::fsubset(type == "WA")
  fm[[s]] <- tmp |>collapse::fsubset(type == "FM")
  #rem.age[[s]] <- stocks[[s]]|> collapse::fsubset(type == "catch")
  tl <- rep(unique(tmp$TL),nrow(waa[[s]]))
  tc <- rep(unique(tmp$troph.cat),nrow(waa[[s]]))
  
  #if(s == "ICES-HAWG_NS_Ammodytes_tobianus") waa[[s]]$value <- waa[[s]]$value/1000
  
  bm[[s]] <- data.frame(Year = num[[s]]$Year,Stock = num[[s]]$Stock,age = num[[s]]$age,
                        bm = num[[s]]$value*waa[[s]]$value,
                        #lam.no.fish.ltr = lambs[[s]]$lam.no.fish,
                        #lam.fish.ltr = lambs[[s]]$lam.fish,
                        catch.bm = (1-exp(-fm[[s]]$value))*waa[[s]]$value*num[[s]]$value,
                        #er = 1-exp(-fm[[s]]$value,
                        num = num[[s]]$value,
                        trophic = tl,
                        troph.cat = tc,
                        Species = num[[s]]$Gen.Spec,
                        Stock.short = num[[s]]$Stock.short,
                        common = num[[s]]$common)
  #Need to clip out the years we don't have biomass data for...
  #bm[[s]] <- bm[[s]] |> collapse::fsubset(Year %in% years[[s]])
  #pnm[[s]] <- 1-exp(-lambdas[[s]]$nm.opt)
  #mx[[s]] <- lambdas[[s]]$fecund.opt
  #vpa[[s]] <- lambdas[[s]]$res$est.abund
} # end for(s in  stock.eco)
# Combine the biomass and abundance data into a dataframe
bm.tst <- do.call("rbind",bm)
bm.tst$troph.cat <- factor(bm.tst$troph.cat,levels = c("≤ 4.0","4.1-4.9","≥ 5.0"),labels = c("Low", "Medium", "High"))

# Look at the biomass and abundance in the ecosystem
# FIX, about 1% of the catch biomasses are larger than the actual biomass observed, take a look
# and make sure that there isn't something mis-aligned for one of the stocks.
bm.tmp <- bm.tst |> collapse::fgroup_by(Stock,Year,trophic,Species,troph.cat,Stock.short,common) |> 
                    collapse::fsummarize(bm = sum(bm,na.rm=T), #+ sum(catch.bm,na.rm=T),
                       num = sum(num,na.rm=T),
                       catch = sum(catch.bm,na.rm=T))# sum(catch.num,na.rm=T),

bm.tst <- NULL
for(s in eco.stocks)
{
  tmp <- bm.tmp[bm.tmp$Stock ==s,]
  tmp$lam.no.fish <- c((tmp$bm[2:nrow(tmp)] + tmp$catch[1:(nrow(tmp)-1)]) / tmp$bm[1:(nrow(tmp)-1)],NA)
  tmp$lam.fish <- c(tmp$bm[2:nrow(tmp)] / tmp$bm[1:(nrow(tmp)-1)],NA)
  ltr.years <- res.lambda.final$year[res.lambda.final$stock==s]
  tmp.years <- tmp$Year
  lam.locs <- which(tmp.years %in% ltr.years)
  tmp$lam.fish.ltr <- NA
  tmp$lam.no.fish.ltr <- NA
  tmp$lam.no.fish.ltr[lam.locs] <- res.lambda.final$lam.no.fish[res.lambda.final$stock==s]
  tmp$lam.fish.ltr[lam.locs] <- res.lambda.final$lam.fish[res.lambda.final$stock==s]
  bm.tst[[s]] <- tmp
}

bm.tot <- do.call('rbind',bm.tst)
#catch = sum(catch.bm,na.rm=T))

#bm.tot$er <- bm.tot$catch/bm.tot$bm # I could add catch back into the denominator of this equation...

# The 'ecosystem' biomass and numbers
com.tot.bm <- bm.tot |> collapse::fgroup_by(Year) |> 
  collapse::fsummarize(num.eco = sum(num),bm.eco = sum(bm))

com.tot.bm$prop.max <- com.tot.bm$bm.eco/max(com.tot.bm$bm.eco)
# Trophic level biomass and numbers.
trophic.bm <- bm.tot |> collapse::fgroup_by(Year,troph.cat) |> 
  collapse::fsummarize(num.tl = sum(num),bm.tl = sum(bm))

tl.com.bm <- left_join(trophic.bm,com.tot.bm,by="Year")
tl.com.bm$prop.bm.tl <- tl.com.bm$bm.tl/tl.com.bm$bm.eco
tl.com.bm$prop.num.tl <- tl.com.bm$num.tl/tl.com.bm$num.eco
# So now take the bm.tot and merge that with the total biomass and the trophic level biomass so we can
# look at what the stock does within it's TL.


# Now we combine the ecosystem results with the stock biomass's

bm.final <- left_join(bm.tot,tl.com.bm,by=c("Year","troph.cat"))
names(bm.final) <- c("Stock","Year","trophic","Species","troph.cat","Stock.short",'common',"bm.stock","num.stock","catch",'lam.no.fish',
                     'lam.fish',"lam.fish.ltr","lam.no.fish.ltr","num.tl","bm.tl",'num.eco','bm.eco','prop.max','prop.bm.tl','prop.num.tl')
# Get the proportion of the total biomass each stock accounts for
bm.final <- bm.final |> collapse::fmutate(prop.bm.stock.eco = (bm.stock+catch)/bm.eco,
                                          prop.num.stock.eco = num.stock/num.eco,
                                          prop.bm.stock.tl = (bm.stock+catch)/bm.tl,
                                          prop.num.stock.tl = num.stock/num.tl)
# Remove 0s from the data
bm.final <- bm.final[bm.final$bm.stock > 0,]
bm.final <- as.data.frame(bm.final)
#bm.final$troph.cat <- factor(bm.final$troph.cat,levels = c("Low","Medium","High"))
# This gets the average weight of individuals in each stock, we'll need this later to get an approximate exploitation rate
bm.final$avg.weight <- bm.final$bm.stock/bm.final$num.stock
# Exploitation Rate
bm.final$er <-  bm.final$catch/(bm.final$catch + bm.final$bm.stock)
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
bm.tl.plt <- ggplot(bm.best) + geom_line(aes(x=Year,y=bm.tl/1e6,group=troph.cat,color=troph.cat)) + 
  scale_x_continuous(name="",breaks=seq(1970,2200,by=5))+
  scale_color_manual(values = c("blue","red","darkgrey","lightgreen")) + 
  scale_y_log10(name="Biomass Pool (millons of tonnes)") + theme(legend.title = element_blank()) 
save_plot(paste0(repo.loc,"/Figures/Historic_Biomass_by_trophic_level.png"),bm.tl.plt,base_height = 8,base_width = 11)
# This is real good now...
prop.bm.tl.plt <-  ggplot(bm.best) + geom_line(aes(x=Year,y=prop.bm.tl,group=troph.cat,color=troph.cat)) + 
  scale_color_manual(values = c("blue","red","darkgrey","lightgreen")) + scale_x_continuous(name="",breaks=seq(1970,2200,by=5))+
  scale_y_continuous(name="Proportion of Biomass Pool") + theme(legend.title = element_blank()) 
save_plot(paste0(repo.loc,"/Figures/Historic_Prop_biomass_by_trophic_level.png"),prop.bm.tl.plt,base_height = 8,base_width = 11)

# The biomass for the ecosystem
bm.com.plt <- ggplot(bm.best) + geom_line(aes(x=Year,y=bm.eco/1e6)) + scale_x_continuous(name="",breaks=seq(1970,2200,by=5))+
  scale_y_continuous(name="Biomass Pool (millons of tonnes)",limits = c(0,NA))
save_plot(paste0(repo.loc,"/Figures/Historic_Biomass_ns_ecosystem.png"),bm.com.plt,base_height = 8,base_width = 11)

# Now get linear model predictions for the lambda-biomass for each stock...
lm.pred.list <- NULL
for(s in eco.stocks)
{
  tmp <- bm.final[bm.final$Stock == s,]
  tmp$bm.prop <- tmp$bm.stock/max(tmp$bm.stock,na.rm=T)
  #mod.tmp <- glm(lam.no.fish.ltr ~ bm.prop,data = tmp,family = gaussian(link = "log"))
  mod.tmp <- lm(log(lam.no.fish.ltr) ~ bm.prop,data = tmp)
  p.val <- summary(mod.tmp)$coefficients[2,4]
  pred.tmp <- data.frame(lambda = NA,bm.prop = seq(0,1,by=0.01),Stock = s)
  if(p.val < 0.2)
  {
    print(summary(mod.tmp)$sigma/sqrt(nrow(tmp)))
    pred.tmp$lambda <- exp(predict(mod.tmp,newdata = pred.tmp))
    pred.tmp$sd <- summary(mod.tmp)$sigma#/sqrt(nrow(tmp))
  } else {
    pred.tmp$lambda <- exp(mean(log(tmp$lam.no.fish.ltr),na.rm=T)) 
    pred.tmp$sd <- sd(log(tmp$lam.no.fish.ltr),na.rm=T)#/sqrt(nrow(tmp))
    } #end else and if statement
  # Shouldn't be necessary, but seems it is to make sure I get a match in the loops...
  pred.tmp$bm.prop <- round(pred.tmp$bm.prop,digits=2)
  pred.tmp$inv.mn.lam <- 1/mean(tmp$lam.no.fish.ltr,na.rm=T)
  lm.pred.list[[s]] <- pred.tmp
}


bm.final |> collapse::fgroup_by(Stock) |> collapse::fsummarise(mn = mean(lam.no.fish.ltr,na.rm=T),
                                                               median = median(lam.no.fish.ltr,na.rm=T),
                                                               inv.mn = 1/mean(lam.no.fish.ltr,na.rm=T),
                                                               inv.med =1/median(lam.no.fish.ltr,na.rm=T))

ggplot(bm.best,aes(x=bm.stock,y=lam.fish.ltr)) + geom_point() + facet_wrap(~Stock,scale='free') + geom_smooth(method='lm') + scale_y_log10()

# get the 

# # You could also get the trophic levels from Fishbase... 
# gen.spec <- unique(paste(ASR.long$Genus,ASR.long$Species,sep=" "))
# library(rfishbase)
# 
# # In fishbase "Limanda ferruginea" is Myzopsetta ferruginea, so add that, then change the name back to what we've been using.  Seems both are in use but Limanda feels more common.
# # here's the trophic levels from Fishbase, looks like we have 2 trophic levels with 5 species in each
# all.troph <- ecology(c(gen.spec,"Myzopsetta ferruginea"), fields=c("Species","SpecCode", "FoodTroph", "FoodSeTroph", "DietTroph", "DietSeTroph"))
# all.troph$Species[all.troph$Species == "Myzopsetta ferruginea"] <- "Limanda ferruginea"
# 
# # A few species pop up multple times for some reason
# mutli.entries <- names(which(table(all.troph$Species) > 1))
# # Clean up the multiple entries
# if(length(mutli.entries) > 0)
# {
#   drop <- NULL
#   for(i in 1:length(mutli.entries))
#     drop[[i]] <- which(all.troph$Species == mutli.entries[i])[-1] # Drop all but the first one
# }
# drop <- do.call('c',drop)
# 
# all.troph <- all.troph[-drop,]
# # ns.troph <- all.troph |> collapse::fsubset(Species %in% gen.spec)
# # Add the trophic level information to ASR.long
# ASR.long <- left_join(ASR.long,all.troph,by = join_by(Gen.Spec==Species))

# Try the function

# Let's try to lay the groundwork for some HDR's given the historical patterns...
manage.strat <- read_xlsx(paste0(repo.loc,"/Data/management_strategy.xlsx"),sheet="strategy")
manage.strat$use.hcr <- F
#manage.strat$relative.er[manage.strat$troph.cat == "≥ 5.0"] <- 3
manage.strat$relative.er[manage.strat$troph.cat == "≤ 4.0"] <- 0.5
#manage.strat$relative.er[manage.strat$troph.cat == "≥ 5.0"] <- 0.5

#manage.strat$er <- 0
#manage.strat$rp.fun <- "min"
# here are the models we have

#source(paste0(repo.loc,"/Scripts/trophic_model_function.R"))





source(paste0(repo.loc,"/Scripts/Population_dynamics_function.R"))
source(paste0(repo.loc,"/Scripts/NS_catch_function.R"))

source(paste0(repo.loc,"/Scripts/trophic_model_function_no_trophic_interactions.R")) # working ok.
#source(paste0(repo.loc,"/Scripts/trophic_model_function_bottom_up.R")) # 
#source(paste0(repo.loc,"/Scripts/trophic_model_function_top_down.R"))


# There are things wrong with the top down function at the moment, fix it. I think the real
# issue is related to indexing....
#exploit.mn$ex.mn <- 0.01
#Scopthalmus_maximus isn't making sense yet with the simple no trophic interaction model.
# First guess the K for the trophic level is being used instead of the K for the stock


#1 NEED TO FIX THE NO INTERACTION, I BROKE IT MESSNING AROUND WITH REMOVING THE TL 3 species 2 case scenario, since that isn't a thing
#possible I broke it standardizing the proportions to =1.

#"bp_log_normal"
#"log_linear
#"bp_resample"
result <- trophic.mod(dat = bm.best,
                                  n.yrs.proj = 50,
                                  n.sims = 500,
                                  manage = manage.strat, 
                                  repo.loc = repo.loc,
                                  method = 'bp_log_normal',
                                  mod.pred = lm.pred.list)


tst <- bm.best[bm.best$Stock=="ICES HAWG_NS 4,3a,7d_Clupea_harengus",c("bm.stock","lam.no.fish","lam.fish","Year","catch","er")]

bm.best |> collapse::fgroup_by(Stock) |> collapse::fsummarise(mn = mean(lam.fish,na.rm=T),
                                                              med = median(lam.fish,na.rm=T),
                                                              min = min(lam.fish,na.rm=T),
                                                              max = max(lam.fish,na.rm=T))
            
bester <- bm.best |> collapse::fgroup_by(Stock) |> collapse::fsummarise(mn = mean(lam.no.fish.ltr,na.rm=T),
                                                                        med = median(lam.no.fish.ltr,na.rm=T),
                                                                        min = min(lam.no.fish.ltr,na.rm=T),
                                                                        max = max(lam.no.fish.ltr,na.rm=T))
              
simmer <- result$sim.ts |> collapse::fgroup_by(Stock) |> collapse::fsummarise(mn.lam = mean(lambda,na.rm=T),
                                                                              med.lam = median(lambda,na.rm=T),
                                                                              mn.ex = mean(ex.rate,na.rm=T))
simmer$mn.lam/bester$mn
simmer$med.lam/bester$med
summary(simmer$med.lam/bester$med)

ggplot(bm.best) + geom_point(aes(x=bm.stock,y=lam.fish)) + facet_wrap(~Stock,scale='free')
dd <- result$sim.ts



dd |> collapse::fgroup_by(Stock) |> collapse::fsummarise(mn = mean(lambda,na.rm=T),
                                                         med = median(lambda,na.rm=T),
                                                         sd = sd(lambda,na.rm=T))

ks <- result$sim.pool.stocks
ks <- result$Pool.realized
ks |> collapse::fgroup_by(Stock) |> collapse::fsummarise(mn = mean(pool.real/pool.init.stock))

ggplot(ks) + geom_line(aes(x=Years,y=pool.real/pool.init.stock,group = sim,color=sim),alpha=0.8) +
  facet_wrap(~interaction(Stock.short, troph.cat , sep = " TC: ",drop=T),scales='free_y') +
  #geom_line(data=bm.best,aes(x=Year,y=bm.stock/1e3)) +
  #scale_x_continuous(name='',breaks = few.breaks) +
  scale_y_continuous(name = "Biomass Pool (thousands of tonnes)") + 
  theme(legend.position = 'none') 

ggplot(ks) + geom_line(aes(x=Years,y=,group = sim,color=sim),alpha=0.8) +
  facet_wrap(~interaction(Stock.short, troph.cat , sep = " TC: ",drop=T),scales='free_y') +
  #geom_line(data=bm.best,aes(x=Year,y=bm.stock/1e3)) +
  #scale_x_continuous(name='',breaks = few.breaks) +
  scale_y_continuous(name = "Biomass Pool (thousands of tonnes)") + 
  theme(legend.position = 'none') 


ggplot(dd) + geom_point(aes(x=bm,y=lambda),alpha=0.8,color='grey') + 
  facet_wrap(~interaction(Stock.short, troph.cat , sep = " TC: ",drop=T),scales='free') +
  geom_point(data=bm.best,aes(x=bm.stock,y=lam.fish.ltr),color='black') +
  #scale_x_continuous(name='',breaks = few.breaks) +
  scale_y_log10(name = "Lambda") + scale_x_log10(name = "biomass") +
  theme(legend.position = 'none') 

ggplot(old.all) + geom_point(aes(x=bm.stock,y=lam.no.fish),alpha=0.8) +
  facet_wrap(~interaction(Stock, troph.cat , sep = " TC: ",drop=T),scales='free') +
  #geom_line(data=bm.best,aes(x=Year,y=bm.stock/1e3)) +
  #scale_x_continuous(name='',breaks = few.breaks) +
  scale_y_continuous(name = "Lambda") + 
  theme(legend.position = 'none') 


# I am getting some impossible results with exploitation rates > 1.... so I am going to have to do the 
# relative increase on the F instead of the exploitation rate!

summary(result$sim.ts$ex.rate[result$sim.ts$troph.cat == "≥ 5.0"])

tst <- result$sim.ts[result$sim.ts$species == "Scopthalmus maximus",]
tst <- result$sim.ts[result$sim.ts$species == "Ammodytes marinus",]
tst <- result$sim.ts[result$sim.ts$species == "Clupea harengus",]

ggplot(tst) + geom_line(aes(x=Years,y=bm,group=as.character(sim),color=as.character(sim))) + scale_y_log10()

summary(tst$K.bm)

tst <- result$sim.K.stocks[result$sim.K.stocks$Species == "-IV 3a,7d_Clupea_harengus",]
summary(tst$cor.prop.bm)

ggplot(result$sim.ts) + geom_line(aes(x=Years,y=ex.rate,group=as.character(sim),color=as.character(sim))) + facet_wrap(~Stock)

tmp <- result$sim.ts |>collapse::fgroup_by(Stock) |> collapse::fsummarise(mn = mean(ex.rate))
mean(tmp$mn)
ggplot(result$sim.troph.K) + geom_line(aes(x=Years,y=bm.tl,group=as.character(sim),color=as.character(sim))) + scale_y_log10() + facet_wrap(~troph.cat,scales='free_y')

hist(eco.lambdas$`ICES-HAWG_ NS-IV 3a,7d_Clupea_harengus`$lam.no.fish)


dd.1 <- dd |> collapse::fsubset(sim == 1)
dd.mer <- dd.1[dd.1$Stock ==s,]

plot(dd.mer$bm/dd.mer$K.bm, dd.mer$lambda

     
dat <-  result$sim.ts |> collapse::fsubset(sim == 1 & Stock == s)

ggplot(dat) + geom_point(aes(x=Years,y=lambda))
ggplot(dat) + geom_point(aes(x=bm,y=lambda))
ggplot(dat) + geom_line(aes(x=Years,y=bm))
ggplot(dat) + geom_point(aes(x=K.bm,y=lambda))
ggplot(dat) + geom_point(aes(x=K.bm,y=bm)) + geom_abline(intercept = 0,slope=1)
ggplot(dat) + geom_point(aes(x=bm/K.bm,y=lambda))
ggplot(dat) + geom_point(aes(x=bm/K.bm,y=K.bm))

one <- result$sim.ts |> collapse::fgroup_by(Stock) |> collapse::fsummarise(mn.k = mean(K.bm,na.rm=T))
two <- bm.best |> collapse::fgroup_by(Stock) |> collapse::fsummarise(mn.bm = mean(bm.stock,na.rm=T))
three <- result$Pool.realized |> collapse::fgroup_by(Stock) |> collapse::fsummarise(mn.k = mean(pool.init.stock,na.rm=T))
two$mn.bm/three$mn.k
two$mn.bm/one$mn.k
one$mn.k/three$mn.k

