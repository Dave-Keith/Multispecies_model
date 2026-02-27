# This function projects the populations 'forward' (a similar function can be used to project 'backwards' from a starting point).

# All the data you need to make it dance...
# years:      The years you are running the simulation for
# nm:         Natural mortality, this wants the instantaneous natural mortality, not the proportional. Expects a matrix (dataframe) with 
#             a column for each age class, and a row for each year.  Should have same number of rows as there are years, 
#             and same number of columns as in fecund.
# fecund:     The fecundity, this is the average number of offspring produced by a female in a given year.Expects a matrix (dataframe) with 
#             a column for each age class, and a row for each year.  Should have same number of rows as there are years, 
#             and same number of columns as in nm.
# ages:       The age classes, should be a vector that has the same length as the columns of nm and fecund
# rems        Removals.  For retro simultions you just have a vector of the 'realized' removals here and that is used. For the projections
#             this is set up so that you provide the mean and sd of fishing mortality for the forward projections. 
#             also added an option ("R_based") to make F a function of R from the euler-lotka model, this give general stability in the model.

# K:          If you have a logistic or 'bounded exponential' model then you'll need to give it a K (carrying capacity). Can be a vector the same length
#             as the number of years in simulation, or can be a single fixed value.
# N.start:    The abundance in the first year
# pop.model:  What method you going to use to get the population growth modeled. You can use exponential model, bounded exponential model, or a logistic model
# n.sims:     The number of simulations you are running
# sim:        Are you doing a retrospecitive "retro" or a projection 'project' simulations. Main difference is really that you don't know removals in the
#             project case, but do know projections in the retrospective.
# repo.loc    Where is your repo.  Deaults to pulling from online Github Repo using 'repo'.  If you put in the directory to point at
#             like D:/Github/ICM that'll work, or you can go with 'preload', which means you have the necessary functions already loaded 


simp.for.sim<-function(years,nm=NULL,fecund = NULL,ages =NULL,rems,K=NULL,N.start = NULL,
                       pop.model = "exponential", n.sims=1,sim = 'retro', repo.loc = 'repo')
{

  if(repo.loc == 'repo')
  {
  # Download the function to go from inla to sf
  funs <- c("https://raw.githubusercontent.com/dave-keith/ICM/main/Scripts/functions/simple_Lotka_r.r",
            "https://raw.githubusercontent.com/dave-keith/ICM/main/Scripts/functions/forward_project.r")
  
  # Now run through a quick loop to load each one, just be sure that your working directory is read/write!
  for(fun in funs) 
  {
    download.file(fun,destfile = basename(fun))
    source(paste0(getwd(),"/",basename(fun)))
    file.remove(paste0(getwd(),"/",basename(fun)))
  }} # close the load from repo....
  
  if(repo.loc != 'repo' & repo.loc != 'preload')
  {
  source(paste0(repo.loc,"/Scripts/functions/simple_Lotka_r.r"))
  source(paste0(repo.loc,"/Scripts/functions/forward_project.r"))
  } # close the load from a location
  
  #st.time <- Sys.time()
  # In case I try to be lazy and shorten names...
  if(pop.model == 'exp') pop.model <- 'exponential'
  if(pop.model == 'log') pop.model <- 'logistic'
  # If doing the projections we start from the year before the simulations period, so we add in the N.start year.
  #if(years > 1)  
  years <- c(min(years)-1,years)
  #Initialize a bunch of objects
  n.years<-length(years)
  Pop<-data.frame(abund = rep(NA,n.years*n.sims),sim = rep(1:n.sims,n.years),years = sort(rep(years,n.sims)))   
  Pop.vec<-rep(NA,n.years)
  rem<-rep(NA,n.years)
  r.vec<-NULL
  B1.vec<-rep(0,n.sims) ### backwards projections 
  # If you have just a single value of K, repeat that number for all the years.
  if(length(K) == 1) K <- rep(K,n.years)

  #Calculate r for your example
  for(i in 1:n.sims)
  {
    #browser()
    # Using the fecundity and mortality data get the lotka r's for the stock
    # DK Note: It might make more sense to put the fishing mortality right into here as we've now done with the backwards simulations.
    if(max(years) > 1) junk<-simple.lotka.r(yrs = years,mort = nm,ages=ages,fecund=fecund)   
    if(max(years) == 1) junk<-simple.lotka.r(yrs = years,mort = c(nm),ages=ages,fecund=c(fecund))
    tmp <-junk$res[,2] 
    if(length(tmp) == 1) tmp <- rep(tmp,n.years)
    r.vec[[i]] <- data.frame(r = c(tmp[-length(tmp)],NA),years = years[],n.sims=i) # How I have it set up the last entry should be an NA  as we don't use that
    #browser()
  } # for(i in 1:n.sims)
  #browser() 
  #print(st.time2 <- Sys.time())
  #unwrap your r vector
  r.vec <- do.call('rbind',r.vec)
  
  #temp.r.vec<-r.vec[r.vec>0 & r.vec<r.cutoff]
  #r.vec<-temp.r.vec[1:n.sims]
  ############################################################
  
  #then do for each population simulation
  
  for(i in 1:n.sims)
  {
    # Get the final year estimate of your population
    Pop.vec[1] <- pop.last <- N.start  # Assuming this is our known starting point, could add uncertainty to this as well if we want to.
    removals <- NA
    r.tmp <- r.vec[r.vec$n.sims == i,]
    #browser()
    #r.store <- r.vec[r.vec$n.sims == i,-1] # We don't use the first r value 
    #Now run your model backwards with the r from the Lotka function and
    # if you use the logistic growth model the K estimated for the population.
    
    # Get your fishing mortality, if rems[2] = 0 this means you exactly know your fishing mortality.
    if(sim == 'project') 
    {
      if(is.numeric(rems[[1]])) fm <- rlnorm(n.years,log(rems[[1]]),rems[[2]])
      
      # Now we can calculate removals.
      if(rems[[1]] == "R_based")
      {
        # So here we harvest some percentage of the long term R estimate for the stock
        # Going to rescale r to be instantaneous, won't have much different for most stocks
        # but will help keep F reasonable for the highly productive stocks
        # Note that I'm using this as proportional removals hereafter (rescaling back properly would lead to fm values > 1 for some stocks)
       
        r.insta <- 1-exp(median(-r.tmp$r,na.rm=T))
        if(r.insta < 0) r.insta <- 0.02 # If this negative, make mortality like 2% for this scenario.
        mn.fm <- rems[[2]]*r.insta
        sd.r <- sd(1-exp(-r.tmp$r),na.rm=T) # This is a bit weird because of the negatives, so will need to think on that.
        fm <- c(NA,rlnorm((n.years-1),log(mn.fm),sd.r))
      } # end if(rems[[1]] == "R_based")
    } # end if(sim == "project")

    
    for(y in 2:n.years)
    {
      # For retro we 'know' removals, for project we set removals based on a fishing mortality
      if(sim == 'retro') removals.next <- rems[y]
      if(sim == 'project') 
      {
        # FIX: I've added in harvest controls for the projections, basically if population is < 20% of K (B0), then harvesting rate declines by 90%
        # This is something to discuss and make customizable.
        if(!is.null(K))
        {
          if(pop.last < 0.2*K[y-1]) removals.next <- 0.1*fm[y-1]*pop.last
          if(pop.last >= 0.2*K[y-1]) removals.next <- fm[y-1]*pop.last
        } 
        # If you have not supplied a carrying capacity
        if(is.null(K)) removals.next <- fm[y-1]*pop.last
        
      } # end the project sims 
     
      r.up <- r.tmp$r[y-1] # So we grab the r associated with the lead in year so r aligns with the initial population numbers
      # The exponential model
      if(pop.model == 'exponential') 
      {
        exp.res <- for.proj(option = "exponential",pop.last = pop.last -removals.next ,r=r.up)
        if(exp.res$Pop.current < 0) exp.res$Pop.current =0 # don't let it drop below 0
        pop.last <- exp.res$Pop.current 
      } # end if(pop.model == 'exponential') 
      
      # So this one is the exponential, but when above whatever bound you have set the population growth rate becomes more negative
      # the more above K you are
      if(pop.model == 'bounded_exp') 
      {
        # If the populations was above the carrying capacity, we force r to suck based on 'how above' the population is
        if(pop.last > K[y-1] & r.up > 0) 
        {
          # FIX?: This is what we assume happens to the population growth rate when the population is above the carrying capacity
          # I think it has nice properties (r tends to be slightly negative just above K, and r quickly declines as K increases)
          # But also r doesn't get ridiculously negative (asymptotes around -0.63). But also, this is just DK freelancing
          ratio <- K[y-1]/pop.last-1
          r.up <- rlnorm(1,ratio,0.05)-1
        }
        exp.res <- for.proj(option = "exponential",pop.last = pop.last -removals.next,r=r.up)
        #browser()
        if(exp.res$Pop.current < 0) exp.res$Pop.current =0 # don't let it drop below 0
        pop.last <- exp.res$Pop.current 
      }
      
      # If you are running the logistic model it's all pretty straightforward
      # FIX: Note in this formulation the 'r' calculated by the euler-lotka model is always down weighted because of K
      # for this reason I don't like this logistic model given we have 'realized' values of r from the past and thus I prefer the bounded exponential
      if(pop.model == 'logistic')
      {
        #browser()
        exp.res <- for.proj(option = "logistic",pop.last = pop.last -removals.next,K=K[y-1],r=r.up)
        if(exp.res$Pop.current < 0) exp.res$Pop.current =0 # don't let it drop below 0
        pop.last <- min(exp.res$Pop.current)
      }
     
      r.vec$r[r.vec$n.sims == i][y-1] <- r.up # Replace it with the realized r.vec
      Pop.vec[y] <- pop.last
      removals[y-1] <- removals.next 

    } #Loop through all the years.
    #browser()
    # Extract the data
    Pop$abund[Pop$sim == i]<-Pop.vec
    Pop$removals[Pop$sim == i]<- c(removals,NA)

  } #end backwards projection
  #print(Sys.time() - st.time2)

    return(list(Pop=Pop,r = r.vec))
}
  
  