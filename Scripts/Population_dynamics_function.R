# now have 3 methods, "bp_resample", ""bp_log_normal"

pop.dam <- function(stock.dat = bm.ts.stock,stock.pool = stock.pool,bm.start=bm.start, catch = catch,stock = s,
                    low.vs.high = 0.4,method="not_sample",mod.preds = NULL)
{

  #browser()
  # Sort out which of the years are low or high bm
  # I'm using 0.5 as the cut off, other options are valid (0.4 is my fav...)
  #low.vs.high <- low.vs.high
  #browser()
  lambda.bump <- mean(stock.dat$er,na.rm=T)
  lambda.bump <- 0
  min.lam <- min(stock.dat$lam.no.fish.ltr,na.rm=T)
  max.lam <- max(stock.dat$lam.no.fish.ltr,na.rm=T)
  mn.lam <- mean(stock.dat$lam.no.fish.ltr,na.rm=T)
  #stock.pool <- stock.pool *(1+2*lambda.bump)
  low.vs.high.bm <- low.vs.high * max(stock.dat$bm.stock)
  # If the stock is in a low vs high 'phase' is determined by
  # the biomass relative to the stock pool.  In this scneario I am going to say low is 40% of the current pool allocation.
  # Don't use this I dont't think, seems to lead to a trap where stocks get stuck where they are.
  low.relative.to.bm.pool <- low.vs.high*stock.pool
  # Have to drop the final year because we don't have a lambda estimate for the final year
  low.bm.years <- which(stock.dat$bm.stock[-nrow(stock.dat)] < low.vs.high.bm)
  if(length(low.bm.years) == 0) low.years <- F else low.years <- T
  high.bm.years <- which(stock.dat$bm.stock[-nrow(stock.dat)] >= low.vs.high.bm)
  
  
  #if(stock == "ICES WGNSSK_NS 4_Scopthalmus_maximus") browser()
  # So first, get a sample 
  method <- method

  if(method == 'bp_resample')
  {
    # Pick one of these to sample if that's how we want to roll, if we have low biomass years (as
    # defined by our cut off low.vs.high)
    if(low.years == T)
    {
      #browser()
      if(bm.start < low.vs.high.bm) samp <- sample(low.bm.years,1)
      if(bm.start >= low.vs.high.bm & bm.start <= stock.pool) samp <- sample(high.bm.years,1)
    } # end If we have low years
    if(low.years == F) samp <- sample(high.bm.years,1)
    # The simple way to do it is just to sample from the natural mortality distribution
    # Now using the right lambda, go look at trends from the stocks that are declining to see what's up there.
    if(bm.start <= stock.pool) 
    {
      lam.samp <- stock.dat$lam.no.fish.ltr[samp] # Get the sample years.  
      # constraining the world to not increase by more than an order of magnitude in a year
      while(lam.samp >(1.5*max.lam) | lam.samp < (0.9*min.lam)) lam.samp <- stock.dat$lam.no.fish.ltr[samp] # Get the sample years.  
    }
    # Still to this the same way if above stock.pool
    if(bm.start > stock.pool) 
    {
      # lam.mn <- mean(stock.dat$lam.no.fish.ltr[high.bm.years],na.rm=T)
      # lam.sd <- sd(log(stock.dat$lam.no.fish.ltr[high.bm.years]),na.rm=T)
      # lam.samp <- rlnorm(1,log(lam.mn),lam.sd)
      lam.mn <- stock.pool/bm.start
      lam.sd <- mod.preds$sd[1]
      lam.samp <- rlnorm(1,log(lam.mn),lam.sd)
      #if(is.na(lam.samp)) browser()
      # If above K, you tend to decline
      if(lam.mn > 0.25) while(lam.samp > mn.lam | lam.samp < (0.5*min.lam)) lam.samp <- rlnorm(1,log(lam.mn),lam.sd)
      if(lam.mn < 0.25) while(lam.samp > mn.lam | lam.samp < 0.25) lam.samp <- rlnorm(1,0.5,lam.sd)
    } # end if(bm.start >= stock.pool) 
  } # end the sample method.
  
  # Or do it the fun way...
  if(method == "bp_log_normal")
  {
    
    # The fun way to do it is to do something multivariate! Note these are instantaneous now!!
    #browser()
    if(bm.start < low.vs.high.bm) 
    {
      if(length(low.bm.years) >0)
      {
        lam.mn <- mean(stock.dat$lam.no.fish.ltr[low.bm.years],na.rm=T)*(1+lambda.bump)
        lam.sd <- sd(log(stock.dat$lam.no.fish.ltr[low.bm.years]),na.rm=T)
        if(length(low.bm.years) == 1) lam.sd <- 0.2 # In case there is just one low biomass year
      }
      if(length(low.bm.years) ==0)
      {
        lam.mn <- mean(stock.dat$lam.no.fish.ltr,na.rm=T)*(1+lambda.bump)
        lam.sd <- sd(log(stock.dat$lam.no.fish.ltr),na.rm=T)
      }
      lam.samp <- rlnorm(1,log(lam.mn),lam.sd)
      while(is.na(lam.samp)) lam.samp <- rlnorm(1,log(lam.mn),lam.sd)
      while(lam.samp >(1.5*max.lam) | lam.samp < (0.9*min.lam)) lam.samp <- rlnorm(1,log(lam.mn),lam.sd)
    } # end if(bm.start < low.vs.high.bm) 
   
    if(bm.start >= low.vs.high.bm & bm.start <= stock.pool) 
    {
      lam.mn <- mean(stock.dat$lam.no.fish.ltr[high.bm.years],na.rm=T)*(1+lambda.bump)
      lam.sd <- sd(log(stock.dat$lam.no.fish.ltr[high.bm.years]),na.rm=T)
      lam.samp <- rlnorm(1,log(lam.mn),lam.sd)
      
      while(is.na(lam.samp)) lam.samp <- rlnorm(1,log(lam.mn),lam.sd)
      while(lam.samp >(1.5*max.lam) | lam.samp < (0.9*min.lam)) lam.samp <- rlnorm(1,log(lam.mn),lam.sd)
      
      } # end if(bm.start < low.vs.high.bm) 
    
   # end if(method != "sample")
  #if(is.na(lam.samp)) browser()
  #while(is.na(lam.samp)) lam.samp <- rlnorm(1,log(lam.mn),lam.sd)
  # Final one, if we are above the K, we are just doing the high biomass scenario for now.
  # Solution, sample from the lambdas at high biomass, but only take lambdas that are <= 1
    if(bm.start > stock.pool) 
    {
      # lam.mn <- mean(stock.dat$lam.no.fish.ltr[high.bm.years],na.rm=T)
      # lam.sd <- sd(log(stock.dat$lam.no.fish.ltr[high.bm.years]),na.rm=T)
      # lam.samp <- rlnorm(1,log(lam.mn),lam.sd)
      lam.mn <- stock.pool/bm.start
      lam.sd <- mod.preds$sd[1]
      lam.samp <- rlnorm(1,log(lam.mn),lam.sd)
      #if(is.na(lam.samp)) browser()
      if(lam.mn > 0.5) while(lam.samp > mn.lam | lam.samp < (0.5*min.lam)) lam.samp <- rlnorm(1,log(lam.mn),lam.sd)
      if(lam.mn < 0.5) while(lam.samp > mn.lam | lam.samp < 0.25) lam.samp <- rlnorm(1,log(0.5),lam.sd)
    } # end if(bm.start >= stock.pool) 
  } #  if(method == "bp_log_normal")
  #while(is.na(lam.samp)) lam.samp <- rlnorm(1,log(lam.mn),lam.sd)
  
  if(method == "log_linear")
  {
    # If no data, run the linear model and calculate it. This will slow things down a bit as you are running
    # a shit tonne of linear models needlessly.
    if(is.null(mod.preds))
    {
        stock.dat$bm.prop <- stock.dat$bm.stock/max(stock.dat$bm.stock,na.rm=T)
        #mod.tmp <- glm(lam.no.fish.ltr ~ bm.prop,data = tmp,family = gaussian(link = "log"))
        mod.tmp <- lm(log(lam.no.fish.ltr) ~ bm.prop,data = stock.dat)
        p.val <- summary(mod.tmp)$coefficients[2,4]
        mod.preds <- data.frame(lambda = NA,bm.prop = seq(0,1,by=0.01),Stock = s)
        browser()
        if(p.val < 0.2)
        {
          mod.preds$lambda <- exp(predict(mod.tmp,newdata = mod.preds))
          mod.preds$sd <- summary(mod.tmp)$sigma/sqrt(nrow(stock.dat))
        } else {
          mod.preds$lambda <- exp(median(log(stock.dat$lam.no.fish.ltr),na.rm=T)) 
          mod.preds$sd <- sd(log(stock.dat$lam.no.fish.ltr),na.rm=T)/sqrt(nrow(stock.dat))
        } #end else and if statement
        
      } # end the if you don't have the data...
    
      if(bm.start <=stock.pool)
      {
        bm.prop <- round(bm.start/stock.pool,digits=2)
        lam.mn <- mod.preds$lambda[mod.preds$bm.prop == bm.prop]*(1+lambda.bump)
        lam.sd <- mod.preds$sd[1]
        lam.samp <- rlnorm(1,log(lam.mn),lam.sd)
        while(lam.samp >(1.5*max.lam) | lam.samp < (0.9*min.lam)) lam.samp <- rlnorm(1,log(lam.mn),lam.sd)
        
      }
      if(bm.start >stock.pool)
      {
        # Note I am using lambda fish here, not lambda no fish, this ensures we have values < 1 for all stocks.
        lam.mn <- stock.pool/bm.start
        lam.sd <- mod.preds$sd[1]
        lam.samp <- rlnorm(1,log(lam.mn),lam.sd)
        
        if(lam.mn > 0.25) while(lam.samp > mn.lam | lam.samp < (0.5*min.lam)) lam.samp <- rlnorm(1,log(lam.mn),lam.sd)
        if(lam.mn < 0.25) while(lam.samp > mn.lam | lam.samp < 0.25) lam.samp <- rlnorm(1,0.5,lam.sd)
      } # end if above pool biomass
  } # end if(method == "log_linear")
  
  # Simple way to include removals
  
  # Removals come out before growth
  #removals <- bm.start*catch$ex.next
  #tst.res <- lam.samp*(bm.start) - removals)
  # Or removals come out after growth?
  
  if(bm.start <= stock.pool) 
  {
    #lam.samp <- lam.samp*(1+lambda.bump)
    bm.end <- lam.samp*(bm.start)
  }
  
  if(bm.start > stock.pool) bm.end <- (lam.samp)*(bm.start)
  removals <- bm.end*catch$ex.next
  tst.res <-  bm.end - removals

  
 

return(list(tst.res = tst.res,removals = removals,ex.rate = catch$ex.next,lambda = lam.samp,pool = stock.pool))
} #end function