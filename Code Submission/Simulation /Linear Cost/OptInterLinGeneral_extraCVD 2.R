#### OptInterLin calculate the optimal intervention for a center with given z value ####

OptInterLinGeneral <- function(beta.vec, gamma.vec = 0, lin.cost.coef, x.min, x.max, p.bar = 0.9,  z = 0, intercept = T)
{
if(length(x.max)!=length(x.min)) stop("x.min and x.max are not at the same length")
if(length(z)!=length(gamma.vec)) stop("gamma.vec and z are not at the same length")
if (any(x.min<0)) stop("The intervention must have non-negative values only")
if (any(x.min>=x.max)) stop("x.max must be larger than x.min")
if(intercept==T)
{
  beta0 <- beta.vec[1]
  beta.vec <- beta.vec[-1]
} else
{
  beta0 <- 0
}
if(length(x.max)!=length(beta.vec)) stop("beta.vec and x.max are not at the same length")

n.coeffs <- length(x.max)
est.x.opt <- x.min
cost.effect <- beta.vec/lin.cost.coef
# If All betas are non-positive the optimal intervnetion is x.min
if (all(cost.effect<=0)) {
  warning("Estimated effects are all negative or zero")
  est.x.opt <- rep(0,length(x.min))
  p.opt <- gamma.vec*z
  return(list(est.x.opt =est.x.opt, p.bar = p.bar, p.opt = p.opt))
}
#est.x.opt <- x.min
# If some betas are negative set their value in the estimated optimal intervention to x.min, and add warning
neg.eff <- cost.effect<=0 # which components have negative effect?
if(sum(neg.eff) >0) {
warning("Some estimated effects are negative. Their value is set to x.min")}
pos.eff <- cost.effect>0 # which components have positive effect?

# Setting the baseline mean by summing the intercept, the effect of z, and the minimal values of the intervention
beta0.center <- beta0 + sum(gamma.vec*z) + sum(beta.vec[neg.eff]*x.min[neg.eff])
if (beta0.center + sum(beta.vec[pos.eff]*x.min[pos.eff]) >=p.bar) {
  warning("The minimal intervention achieves at least p.bar")
  return(list(est.x.opt =est.x.opt, 
              p.bar = beta0.center + sum(beta.vec[pos.eff]*x.min[pos.eff]), 
              p.opt = p.bar))
  }
order.effect <- order(cost.effect[pos.eff], decreasing = T)
# Can we get to p.bar if setting all positive effect to x.max?
# I subtract x.min for the pos.eff (cause I have those in the beta0.center)
if(beta0.center + sum(beta.vec[pos.eff]*x.max[pos.eff]) < p.bar) 
  {
  est.x.opt[pos.eff] <- x.max[pos.eff] # Best you can do is to set the intervention to the max
  p.opt <-   beta0.center + sum(beta.vec[pos.eff]*x.max[pos.eff])
  warning("Probability under most effective pacakge is less than p.bar /n Returnning x.max")
} else {
  p.opt <- p.bar
  for (i in 1:max(order.effect))
  {
    beta.temp <- beta.vec[pos.eff][order.effect][1:i] # this beta is the 1 to i--th most effective betas (for i=1 the most effective, for i=2 the two most effective)
    x.max.temp <- x.max[pos.eff][order.effect][1:i] # same for x.max 
    if (i < length(order.effect)) # set other components to x.min
    {
    beta.other.temp <- beta.vec[pos.eff][order.effect][(i+1):length(order.effect)]
    x.other.temp <- x.min[pos.eff][order.effect][(i+1):length(order.effect)]
    } else 
    {
      beta.other.temp <- x.other.temp <- 0 # in case we are in the last componenet
    }
    if(beta0.center + sum(beta.temp*x.max.temp) + sum(beta.other.temp*x.other.temp)<p.bar) {next}
    if(beta0.center + sum(beta.temp*x.max.temp) + sum(beta.other.temp*x.other.temp)==p.bar) {
      est.x.opt[pos.eff][order.effect][1:i] <- x.max.temp
      est.x.opt[pos.eff][order.effect][1:i] <- x.max.temp
      break
    }
    if (i==1)
    {
      est.x.opt[pos.eff][order.effect][i] <- ((p.bar) - beta0.center - sum(x.other.temp*beta.other.temp))/beta.temp[i]
      break 
    } else {
      est.x.opt[pos.eff][order.effect][1:(i-1)] <- x.max.temp[-i] # the top i-1 gets their maximal values
      est.x.opt[pos.eff][order.effect][i] <- ((p.bar)-sum(beta.temp[1:(i-1)]*x.max.temp[1:(i-1)]) - beta0.center - sum(x.other.temp*beta.other.temp))/beta.temp[i]
      break 
    }
  }}
return(list(est.x.opt =est.x.opt, p.bar = p.bar, p.opt = p.opt)) # maybe you want to change it such that p.opt is also returned
  }
