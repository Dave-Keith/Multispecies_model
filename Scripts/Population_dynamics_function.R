pop.dam <- function(lambda = stock.lambdas,stock.K = tl.K,
                    bm.stock = bm.ts.stock,low.vs.high = 0.5,method="not_sample")
{

  #browser()
  # Sort out which of the years are low or high bm
  # I'm using 0.5 as the cut off, other options are valid (0.4 is my fav...)
  #low.vs.high <- low.vs.high
  low.vs.high.bm <- low.vs.high * max(bm.stock$bm.stock)
  # Have to drop the final year because we don't have a lambda estimate for the final year
  low.bm.years <- which(bm.stock$bm.stock[-nrow(bm.stock)] < low.vs.high.bm)
  if(length(low.bm.years) == 0) low.years <- F else low.years <- T
  high.bm.years <- which(bm.stock$bm.stock[-nrow(bm.stock)] >= low.vs.high.bm)
  # ANd get what our carrying capacity is at this moment
  cur.K <- stock.K
  
  # So first, get a sample 
  method <- "not_sample"
  if(method == 'sample')
  {
    # Pick one of these to sample if that's how we want to roll, if we have low biomass years (as
    # defined by our cut off low.vs.high)
    if(low.years == T)
    {
      #browser()
      if(bm.start < low.vs.high.bm) samp <- sample(low.bm.years,1)
      if(bm.start >= low.vs.high.bm) samp <- sample(high.bm.years,1)
    } # end If we have low years
    if(low.years == F) samp <- sample(high.bm.years,1)
    # The simple way to do it is just to sample from the natural mortality distribution
    # Now using the right lambda, go look at trends from the stocks that are declining to see what's up there.
    lam.samp <- lambda$lam.no.fish[samp] # Get the sample years.  
  } # end the sample method.
  
  # Or do it the fun way...
  if(method != "sample")
  {
    
    # The fun way to do it is to do something multivariate! Note these are instantaneous now!!
    #browser()
    if(bm.start < low.vs.high.bm) 
    {
      if(length(low.bm.years) >0)
      {
        lam.mn <- mean(lambda$lam.no.fish[low.bm.years],na.rm=T)
        lam.sd <- sd(log(lambda$lam.no.fish[low.bm.years]),na.rm=T)
        if(length(low.bm.years) == 1) lam.sd <- 0.2 # In case there is just one low biomass year
      }
      if(length(low.bm.years) ==0)
      {
        lam.mn <- mean(lambda$lam.no.fish,na.rm=T)
        lam.sd <- sd(log(lambda$lam.no.fish),na.rm=T)
      }
      lam.samp <- rlnorm(1,log(lam.mn),lam.sd)
      while(is.na(lam.samp)) lam.samp <- rlnorm(1,log(lam.mn),lam.sd)
    } # end if(bm.start < low.vs.high.bm) 
   
    if(bm.start >= low.vs.high.bm & bm.start < cur.K) 
    {
      lam.mn <- mean(lambda$lam.no.fish[high.bm.years],na.rm=T)
      lam.sd <- sd(log(lambda$lam.no.fish[high.bm.years]),na.rm=T)
      lam.samp <- rlnorm(1,log(lam.mn),lam.sd)
      
      while(is.na(lam.samp)) lam.samp <- rlnorm(1,log(lam.mn),lam.sd)
      } # end if(bm.start < low.vs.high.bm) 
    
   # end if(method != "sample")
  #if(is.na(lam.samp)) browser()
  #while(is.na(lam.samp)) lam.samp <- rlnorm(1,log(lam.mn),lam.sd)
  # Final one, if we are above the K, we are just doing the high biomass scenario for now.
  # Solution, sample from the lambdas at high biomass, but only take lambdas that are <= 1
  if(bm.start >= cur.K) 
  {
    lam.mn <- mean(lambda$lam.no.fish[high.bm.years],na.rm=T)
    lam.sd <- sd(log(lambda$lam.no.fish[high.bm.years]),na.rm=T)
    lam.samp <- rlnorm(1,log(lam.mn),lam.sd)
    #if(is.na(lam.samp)) browser()
    while(lam.samp >1) lam.samp <- rlnorm(1,log(lam.mn),lam.sd)
  }
  }
  #while(is.na(lam.samp)) lam.samp <- rlnorm(1,log(lam.mn),lam.sd)
  
  # Simple way to include removals
  
  if(is.null(catch$er.mn)) {ex.rate = 0; ex.sd = 0}
  if(!is.null(catch$er.mn))
  {
    er.mn <- catch$er.mn
    if(!is.null(catch$er.sd)) er.sd <- catch$er.sd
    if(is.null(catch$er.sd)) er.sd <- data.frame(er.sd=0,Stock =s)
    ex.rate = er.mn$ex.mn[er.mn$Stock == s]
    ex.sd = er.sd$ex.sd[er.sd$Stock == s]
  } # end if(!is.null(catch$er.mn))
  if(is.null(catch$catch)) 
  {
    #print(ex.rate)
    # Convert proportion to F
    ex.rate <- -log(1-ex.rate)
    er <- rlnorm(1,log(ex.rate),ex.sd)
    # Go from F back to proportion
    er <- 1-exp(-er)
    removals <- bm.start*er
  } # end if(is.null(catch$catch)) 
  
  # If you have a catch estimate for the stock
  if(!is.null(catch$catch))
  {
    limit.er <- 0.4
    removals.tmp <- catch$catch$catch[catch$catch$Stock ==s]
    er <- removals.tmp/(bm.start+removals.tmp)
    if(er > limit.er) 
    {
      removals.tmp <- limit.er*bm.start
      er <- limit.er
    }
    removals <- removals.tmp
  } # end if(!is.null(catch$catch))
  
  
  #print(er)
  
  #lam.samp <- rlnorm(1,log(1),0.1)
  #tst.res <- (lam.samp)*bm.start - removals
  # We want to grow after removals because otherwise we can get negative values given exploitation was
  # calculated using the initial biomass
  tst.res <- lam.samp*(bm.start - removals)
  
 

return(list(tst.res = tst.res,removals = removals,ex.rate = er,lambda = lam.samp,K = stock.K))
} #end function