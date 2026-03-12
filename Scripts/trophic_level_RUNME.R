# FIX For the MAESTRO work, lets simulate the ecosystem a fixed number of times, then run the population dynamics on each of these ecosystem
# simulations, instead of simulating the ecosystem for each of the population dynamics scenarios, for example, we can simulate 1,000 
# ecosystems.  Then, for the fishery dynamics, if we are testing 10 scenarios and simulating the population dynamics 100 times, that 
# all happens within a set of 1000 ecosystems scenarios (so each ecosystems would have 1000 population dynamics simulations run on it in this scenario)


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
for.tune.all <- for.tune.summary
eco.stocks <- names(for.tune.all)


# Get the trophic levels
# Trophic levels from Simon Jennings Paper in 2002 for NS...or fishbase if not in Jennings (e.g. Turbot)
# The TL for Sand eel and witch flounder is 4, but I am lumping 3.8 and 4 together as one, doing so in the code by making the TL be the floor()
# so these as 3.99 makes that work without doing anything greasy.
eco.troph <- read_xlsx( paste0(dat.loc,"/Data/species_TL.xlsx"))
# This gets some of the stock data we need.
stocks <- ASR.long |> collapse::fsubset(Stock %in% eco.stocks)
# Get the trophic levels for the stocks
stocks <- merge(stocks,eco.troph,by="Stock")
# Get the ages for the stocks



# Make it a list
stock.lst <- NULL
exploit <- NULL
eco.lambdas <- NULL
for(s in eco.stocks) 
{
  stock.lst[[s]] <- stocks |> collapse::fsubset(Stock == s)
  exploit[[s]] <- data.frame(stock =s,
                            removals = for.tune.summary[[s]]$removals,
                            er = for.tune.summary[[s]]$removals/(for.tune.summary[[s]]$est.abund+for.tune.summary[[s]]$removals))
  eco.lambdas[[s]] <- res.lambda.final[res.lambda.final$stock ==s,]
}



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

exploit.dat <- data.frame(er.mn = NA,
                          er.sd = 0.1,
                          stock = names(stock.lst),
                          lrp = NA,
                          urp = NA,
                          rr = NA,
                          er.below.lrp = NA)


# here are the models we have
source(paste0(repo.loc,"/Scripts/Population_dynamics_function.R"))
source(paste0(repo.loc,"/Scripts/NS_caatch_function.R"))

#source(paste0(repo.loc,"/Scripts/trophic_model_function.R"))




source(paste0(repo.loc,"/Scripts/trophic_model_function_no_trophic_interactions.R")) # working ok.
source(paste0(repo.loc,"/Scripts/trophic_model_function_bottom_up.R")) # 
source(paste0(repo.loc,"/Scripts/trophic_model_function_top_down.R"))


# There are things wrong with the top down function at the moment, fix it. I think the real
# issue is related to indexing....
#exploit.mn$ex.mn <- 0.01
#Scopthalmus_maximus isn't making sense yet with the simple no trophic interaction model.
# First guess the K for the trophic level is being used instead of the K for the stock

result <- trophic.mod(stocks = stock.lst,
                                  lambdas = eco.lambdas,
                                  n.yrs.proj = 50,
                                  n.sims = 2,
                                  exploit = exploit,
                                  repo.loc = repo.loc)

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


tst <- data.frame(result$sim.K.stocks)

tst |> group_by(Years,sim) |> mutate(sum(bm.stock))