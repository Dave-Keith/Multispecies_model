### FISHING FUNCTION ###

#This function is called in the stock loop of trophic_model_function simulation,
#and determines the catch and net biomass of a given stock in a given projection year

#DK Note
#1 DK has taken Danielle's function and messed around with it.


# DB Notes
#1: mgmt.stock:     # stock-specific management plan details (thresholds, assessment interval, stock ID #, 
#and slope/intercept of equation between u and biomass)
#2: repo.loc:       # Location of the Github repo, defaults to "D:/GitHub/Multispecies_model/"

#NOTE: writing this as if we are fishing population first, then growing population (discrete-time model)
#NOTE 2: each updated exploitation rate is applied not to the current projection year, but to the next one 
#(this delay models lags in mgmt implementation)
#



proj.catch.eqn <- function(dat, bm = NULL,repo.loc = "D:/GitHub/Multispecies_model/")
{

#browser()
# If we are using harvest control rules  
  if(dat$use.hcr == T)
  {
  if (bm <= dat$lrp) 
    {
      if(dat$fm.below.lrp == 0) ex.next <- 0 
      if(dat$fm.below.lrp > 0) ex.next <- 1-exp(-rlnorm(1, mean = log(dat$fm.below.lrp), sd = dat$er.sd))
    } # end the below lrp if
    if(bm > dat$lrp  & bm < dat$urp)  ex.next <- 1-exp(-rlnorm(1, mean = log(dat$fm.mn*(bm-dat$lrp)/dat$lrp + dat$fm.below.lrp) , sd = dat$er.sd))
    if(bm >= dat$urp) ex.next <- 1-exp(-rlnorm(1, mean = log(dat$fm.mn), sd = dat$er.sd))
  }
# If we are not using harvest control rules
  if(dat$use.hcr == F) ex.next <- 1-exp(-rlnorm(1, mean = log(dat$fm.mn), sd =dat$er.sd))
  

  return(list(ex.next = ex.next))
  
}#end function