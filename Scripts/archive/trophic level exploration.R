# Playing around to see how to get the 'trophic levels' and bin the biomass

library(tidyverse)
library(GGally)
library(cowplot)
library(ggthemes)
library(FishLife)
library(rfishbase)
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
load(file = paste0(dat.loc,"/Results/all_cleaned_forward_tune_summaries_no_age_corection_fec_nm.Rdata"))
load(file = paste0(dat.loc,"/Results/model_inputs_no_age_correction.Rdata"))


n.yrs.proj <- 25 # How many years into the future we are going to project the stocks
n.sims <- 10 # The numbers of simulations to run, keeping low for testing...

# Get the right stocks
Stocks <- names(for.tune.all)
ns.stocks <- Stocks[grep("NS",Stocks)]
ns.stocks <- ns.stocks[ns.stocks != "ICES-WGHANSA_SP8abd_Sardina _pilchardus"]

gen.spec <- unique(paste(ASR.long$Genus,ASR.long$Species,sep=" "))

# In fishbase "Limanda ferruginea" is Myzopsetta ferruginea, so add that, then change the name back to what we've been using.  Seems both are in use but Limanda feels more common.
# here's the trophic levels from Fishbase, looks like we have 2 trophic levels with 5 species in each
all.troph <- ecology(c(gen.spec,"Myzopsetta ferruginea"), fields=c("Species","SpecCode", "FoodTroph", "FoodSeTroph", "DietTroph", "DietSeTroph"))
all.troph$Species[all.troph$Species == "Myzopsetta ferruginea"] <- "Limanda ferruginea"

mutli.entries <- names(which(table(all.troph$Species) > 1))

# Herring shows up twice for some reason, drop the second one... make it an if just in case this stops happening...
if(length(mutli.entries) > 0)
{
drop <- NULL
for(i in 1:length(mutli.entries))
drop[[i]] <- which(all.troph$Species == mutli.entries[i])[-1] # Drop all but the first one
}
drop <- do.call('c',drop)

all.troph <- all.troph[-drop,]

# Now stick the trophic level on the ASR data...
tst <- left_join(all.troph,ASR.long,by = join_by(Species==Gen.Spec))

ns.species <- ASR.long$Gen.Spec[ASR.long$Stock %in% ns.stocks]
ns.troph <- all.troph |> collapse::fsubset(Species %in% ns.species)

ggplot(ns.troph) + geom_point(aes(x=FoodTroph,y=DietTroph))

# Going to use DietTroph with trophic levels of 3 and 4.


