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
  #isolate management criteria from stock-specific management plan
  lrp <- dat$lrp
  urp <- dat$urp # Was called high threshold
  er.sd <- dat$er.sd
  er.below.lrp <- dat$er.below.lrp
  rr <- dat$rr
  

  #lognormal sample exploitation rates from distribution of mean around desired value; set sd to 0 to get exact value
  #use rlnorm()

  
  #start "assessment year" loop
  #if(t %% a.interval == 0){
  if (bm <= lrp) 
  {
    if(er.below.lrp == 0) ex.next <- 0 
    if(er.below.lrp > 0) ex.next <- rlnorm(1, mean = log(er.below.lrp), sd = er.sd)
  } # end the below lrp if
  if(lrp <= bm  & bm < urp)  ex.next <- rlnorm(1, mean = log(rr*(bm-lrp)/lrp + er.below.lrp) , sd = er.sd)
  if(bm >= urp) ex.next <- rlnorm(1, mean = log(rr), sd = er.sd)
      

  return(list(ex.next = ex.next))
  
}#end function