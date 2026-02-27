# FIX For the MAESTRO work, lets simulate the ecosystem a fixed number of times, then run the population dynamics on each of these ecosystem
# simulations, instead of simulating the ecosystem for each of the population dynamics scenarios, for example, we can simulate 1,000 
# ecosystems.  Then, for the fishery dynamics, if we are testing 10 scenarios and simulating the population dynamics 100 times, that 
# all happens within a set of 1000 ecosystems scenarios (so each ecosystems would have 1000 population dynamics simulations run on it in this scenario)


n.yrs.proj <- 50 # How many years into the future we are going to project the stocks
n.sims <- 50 # The numbers of simulations to run, keeping low for testing...


#dat.loc <- 'C:/Users/keithd/Documents/GitHub/ICM'
#repo.loc <- "C:/Users/keithd/Documents/GitHub/Multispecies_model"
dat.loc <- 'D:/GitHub/ICM'
repo.loc <- "D:/GitHub/Multispecies_model"

#loc <- "C:/Users/keithd/Documents/Github/ICM"
load(file = paste0(dat.loc,"/Results/all_cleaned_forward_tune_summaries_no_age_corection_fec_nm.Rdata"))
load(file = paste0(dat.loc,"/Results/model_inputs_no_age_correction.Rdata"))
# load in the trophic mocdel function
source(paste0(repo.loc,"/Scripts/trophic_model_function.R"))

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


# Get the right stocks, this is from our model runs
Stocks <- names(for.tune.all)
eco.stocks <- Stocks[grep("NS",Stocks)]
# This one doesn't have the data we need
eco.stocks <- eco.stocks[eco.stocks != "ICES-WGHANSA_SP8abd_Sardina _pilchardus"]

# Get the trophic levels
# Trophic levels from Simon Jennings Paper in 2002 for NS...or fishbase if not in Jennings (e.g. lesser sand eel and Turbot)
eco.troph <- data.frame(Stock = eco.stocks,
                        Common = c("Herring","Lesser Sand eel","Sole","Atlantic cod",
                                   "Haddock","European plaice","Norway pout","Saithe",
                                   "Atlantic cod", "Whiting", "Sole", "European plaice","Turbot","Sole"),
                        TL = c(3.8,3.08,5.0,5.2,
                               4.7,4.5,4.2,4.6,
                               5.2,5.3,5.0,4.5,4.4,5.0))
# This gets some of the stock data we need.
stocks <- ASR.long |> collapse::fsubset(Stock %in% eco.stocks)
# Get the trophic levels for the stocks
stocks <- merge(stocks,eco.troph,by="Stock")
# Get the ages for the stocks



# Make it a list
stock.lst <- NULL
for(s in eco.stocks) stock.lst[[s]] <- stocks |> collapse::fsubset(Stock == s)

eco.lambdas <- NULL
for(s in names(for.tune.all)) if(s %in% eco.stocks) eco.lambdas[[s]] <- res.lambda.final[res.lambda.final$Stock ==s,]

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
exploit.mn <- data.frame(ex.mn = c(0.01, 0.01, 0.30,  0.30, 0.01, 0.01, 0.01, 0.01, 0.30, 0.30, 0.30, 0.01, 0.01, 0.30),Stock = names(stock.lst)) # Fish TL 5 hard
exploit.mn <- data.frame(ex.mn = c(0.01, 0.01, 0.01, 0.01, 0.20, 0.20, 0.20, 0.20, 0.01, 0.01, 0.01, 0.20 ,0.20, 0.01),Stock = names(stock.lst)) # Fish TL 4 hard
exploit.mn <- data.frame(ex.mn = c(0.20, 0.20, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01),Stock = names(stock.lst)) # Fish TL 3 hard
exploit.mn <- data.frame(ex.mn = rep(0.01,length(stock.lst)),Stock = names(stock.lst)) # minimal fishing
exploit.sd <- data.frame(ex.sd = c(rep(0.1,length(stock.lst))),Stock = names(stock.lst))
st.time <- Sys.time()
# test <- trophic.mod(stocks = stock.lst,lambdas= eco.lambdas,n.sims=n.sims,
#                     catch = list(catch =NULL,er.mn = exploit.mn,er.sd = exploit.sd),
#                     n.yrs.proj= n.yrs.proj,repo.loc=repo.loc)
# Sys.time() - st.time
# 
# 
# catch <- data.frame(catch = c(1e5,1e5,100,100,1e4,1000,1e4,1e4,1e4,1e4,100,100,100,100),Stock = names(stock.lst))
# test <- trophic.mod(stocks = stock.lst,lambdas= eco.lambdas,n.sims=n.sims,
#                     catch = list(catch =catch,er.mn = NULL,er.sd = NULL),
#                     n.yrs.proj= n.yrs.proj,repo.loc=repo.loc)

# Look at this...
test$sim.ts


saveRDS(object = test,file = paste0(repo.loc,"/Results/NS_projections_",n.sims,"_sims.Rds"))

haddock.K <- test$sim.K.stocks[test$sim.K.stocks$Stock=="ICES-WGNSSK_NS  4-6a-20_Melanogrammus_aeglefinus",]
haddock.bm <- test$sim.ts[test$sim.ts$Stock=="ICES-WGNSSK_NS  4-6a-20_Melanogrammus_aeglefinus",]

plot(haddock.bm$bm[haddock.bm$sim == 1],type='l')
plot(haddock.K$bm.stock[haddock.K$sim ==1],type='l')

plot(haddock.bm$bm[haddock.bm$sim == 7]/ haddock.K$bm.stock[haddock.K$sim ==7],type='l')                       

# here are the models we have
source(paste0(repo.loc,"/Scripts/Population_dynamics_function.R"))
#source(paste0(repo.loc,"/Scripts/trophic_model_function.R"))

source(paste0(repo.loc,"/Scripts/trophic_model_function_full_cascade.R"))


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
                                  n.sims = 10,
                                  catch = list(catch =NULL,er.mn = exploit.mn,er.sd = exploit.sd),
                                  repo.loc = repo.loc)

tst <- result$sim.ts[result$sim.ts$species == "Scophthalmus maximus",]
tst <- result$sim.ts[result$sim.ts$species == "Ammodytes tobianus",]
tst <- result$sim.ts[result$sim.ts$species == "Clupea harengus",]

ggplot(tst) + geom_line(aes(x=Years,y=bm,group=as.character(sim),color=as.character(sim))) + scale_y_log10()

summary(tst$K.bm)

tst <- result$sim.K.stocks[result$sim.K.stocks$Species == "-IV 3a,7d_Clupea_harengus",]
summary(tst$cor.prop.bm)


ggplot(result$sim.troph.K) + geom_line(aes(x=Years,y=bm.tl,group=as.character(sim),color=as.character(sim))) + scale_y_log10() + facet_wrap(~troph.cat,scales='free_y')

hist(eco.lambdas$`ICES-HAWG_ NS-IV 3a,7d_Clupea_harengus`$lam.no.fish)


tst <- data.frame(result$sim.K.stocks)

tst |> group_by(Years,sim) |> mutate(sum(bm.stock))