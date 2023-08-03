library(poweRlaw)
library(tidyverse)

setwd("~/Work/projects/covid/long-deletions/")

read_data = function(path) {
  df = read_tsv(path)
  df$mut_type = factor(df$mut_type,levels=c("synonymous","missense","nonsense","undoStop"))
  df <- df %>%
    mutate(daysPos = days_circulated + 1)
  return(df)
}

estimate_dist = function(distr,xmin) {
  distr$setXmin(xmin)
  est_distr = estimate_pars(distr)
  distr$setPars(est_distr)
  return(distr)
}

estimate_alpha = function(data){
  muts = list('synonymous','missense','nonsense')
  A = matrix(0,nrow=1,ncol=3)
  colnames(A) = muts
  for (mut in muts) {
    data %>%
      filter(mut_type == mut) %>%
      pull(leaf_count) -> vec
    pl = displ$new(vec)
    pl$setXmin(min(vec))
    alpha = estimate_pars(pl)
    A[1,mut] = alpha$pars
  }
  return(A)
}


compared = function(dist1,dist2,lab1,lab2) {
  comp = compare_distributions(dist1,dist2)
  cat(lab1, "vs.",lab2,"\n")
  print(comp$p_one_sided)
  print(comp$p_two_sided)
}

checkPowerLaw = function(data) {
  pl = displ$new(data)
  bs_pl = bootstrap_p(pl,no_of_sims=500,threads=8,xmins=seq(1,100,1))
  return(plot(bs_pl,trim=0.1))
}

checkOtherDistributions = function(data) {
  pl = displ$new(data)
  
  ## MLE of initial power law
  est = estimate_xmin(pl,xmins=seq(1,100,1))
  print(est)
  pl$setXmin(est)
  est_pl = estimate_pars(pl)
  pl$setPars(est_pl)
  
  # MLE of lognormal 
  lnorm = dislnorm$new(data)
  lnorm = estimate_dist(lnorm,est)
  
  ## MLE of poisson
  pois = dispois$new(data)
  pois = estimate_dist(pois,est)
  
  ## MLE of exponential
  exp = disexp$new(data)
  exp = estimate_dist(exp,est)
  
  ## Compare distributions
  compared(pl,lnorm,"power law", "log normal")
  compared(pl,pois,"power law", "poisson")
  compared(pl,exp,"powr law", "exponential")
  compared(lnorm,pois, "log normal", "poisson")
  compared(lnorm, exp, "log normal", "exponential")
  compared(pois, exp, "poisson", "exponential")
}

# Read in data
orf8 = read_data('usher/trimmed/clades/ORF8_clades.tsv')
spike = read_data('usher/trimmed/clades/S_clades.tsv')
S1 = read_data('usher/trimmed/clades/S1_clades.tsv')
orf1a = read_data('usher/trimmed/clades/ORF1a_clades.tsv')

orf8_nested = read_data('usher/trimmed/clades_nested/ORF8_clades.tsv')
spike_nested = read_data('usher/trimmed/clades_nested/S_clades.tsv')
S1_nested = read_data('usher/trimmed/clades_nested/S1_clades.tsv')
orf1a_nested = read_data('usher/trimmed/clades_nested/ORF1a_clades.tsv')

orf8_alpha = estimate_alpha(orf8_nested)
spike_alpha = estimate_alpha(spike_nested)
S1_alpha = estimate_alpha(S1_nested)
orf1a_alpha = estimate_alpha(orf1a_nested)






# Does any of this look like a power law?
checkPowerLaw(orf8$leaf_count) # NOPE
checkPowerLaw(orf8$daysPos) # NADA
checkPowerLaw(orf8_nested$leaf_count) # also NO
checkPowerLaw(orf8_nested$daysPos) # yah noo

## What distributions do these look like?
checkOtherDistributions(orf8$leaf_count) # Log normal is the best
checkOtherDistributions(orf8$daysPos) # Can't distinguish between lognormal & power
checkOtherDistributions(orf8_nested$leaf_count) # Log normal is the best
checkOtherDistributions(orf8_nested$daysPos) # Can't distinguish between lognormal & power

## Are these log normal?
shapiro.test(log(sample(orf8$leaf_count,size=3000),10)) # NO
shapiro.test(log(sample(orf8$daysPos,size=3000),10)) # NO
shapiro.test(log(sample(orf8_nested$leaf_count,size=3000),10)) # NO
shapiro.test(log(sample(orf8_nested$daysPos,size=3000),10)) # NO

## Okay, so these data remain skewed as hell & they don't fit any standard ass distributions.
## I guess I just do a Mann Whitney U