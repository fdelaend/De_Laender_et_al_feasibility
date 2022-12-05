
# Functions -------------------------------------------------------------
## To compute equilibrium densities -------------------------------------
# Of a 3-sp consumer (n_1) resource (n_2 and n_3) system: 
# b_i is the intrinsic growth rate of the resources and the mortality rate of the consumer
# a is the consumption rate of the consumer on the resources (same on both resources)
# efficiency is 1
get_eq_3CR <- function(b_1, b_2, b_3, a) {
  tibble(n_1 = (b_3 + a*(b_1 - b_2))/2/a,
       n_2 = (b_3 - a*b_1 + a*b_2)/2/a, 
       n_3 = (-b_3 + a*(b_1 + b_2))/2/a^2)
}

# Of a 3-sp food chain with one resource (n_1), one consumer (n_2), and one predator (n_3): 
# b_i are the intrinsic growth rates (positive for the resource; negative for the consumer and predator (mortality))
# a_2 and a_3 are the consumption rates of the consumer n_2 on the resource n_1,
# and of the predator n_3 on the consumer n_2
# efficiency is 1
get_eq_3chain <- function(b_1, b_2, b_3, a_2, a_3) {
  tibble(n_1 = b_1 - a_2*b_3/a_3, 
       n_2 = b_3/a_3, 
       n_3 = (a_2*a_3*b_1 - a_3*b_2 - a_2^2*b_3)/a_3^2)
}

## To make all sorts of matrices -------------------------------------
#Make block diagonal matrix by stacking p times A
# on the diagonal. The matrix A can be simply perfectly copied (vary=0),
# or varied to some extend among copies, as quantified by  
# "vary" = the width of the uniform (max-min) around the elements of A
make_block_diagonal <- function(A, p, vary = 0, ...) {
  A_temp <- array(0, dim=dim(A)*p) # all zeros
  for (i in 1:p) {
    A_vary <- array(data = runif(prod(dim(A)), 1 - vary, 1 - vary), dim = dim(A)) #local deviation of A in block i
    A_temp[(ncol(A)*i-ncol(A)+1):(ncol(A)*i), (ncol(A)*i-ncol(A)+1):(ncol(A)*i)] <- A * A_vary
    }
  return(A_temp)
}
#Make a dispersal (D) matrix
# all diagonal blocks are zero
# off diagonals blocks are diagonal matrices containing the sp-specific dispersal rates
# if two sites are connected, they exchange individuals in both directions (but allowing for different rates)
# d_max = max of uniform dist of dispersal rate 
# p = nr of locations; take at least 10 to avoid problems
# n = nr of sp
# connectivity = *fraction* of all possible combos that are connected
make_D <- function(p=3, n=2, connectivity=0.5, d_max=1e-3, ...){
  D    <- array(0, dim=c(n*p, n*p))#matrix w correct dimensions and all zeros
  pos_links<- combn(c(1:p), 2) # all possible combos
  nr_links <- round(connectivity * ncol(pos_links), 0) #nr of links between locations
  links    <- pos_links[,sample(c(1:ncol(pos_links)), nr_links)]#sample that nr of links
  for (i in 1:ncol(links)) {
    D_12 <- diag(runif(n, 0, d_max))#make random dispersal matrix from 1 to 2
    D_21 <- diag(runif(n, 0, d_max))#make random dispersal matrix from 2 to 1
    pos_1 <- (n*links[1,i]-n+1):(n*links[1,i]) #position of first location
    pos_2 <- (n*links[2,i]-n+1):(n*links[2,i]) #position of second location
    D[pos_1, pos_2]<-D_12
    D[pos_2, pos_1]<-D_21
  }
  return(D)
}

# Make consumer resource matrix: parameters are s_r=nr of resources, s_c=nr of consumers
# a_mean=mean of uniformly distributed attack rates; a_range=range of attack rates (max-min);
# a_0_mean=mean of uniformly distributed competition coefficient; a_0_range=range of comp.coeff (max-min);
# Random among-resource interaction matrices are made until a negative def. is found.
make_A_CR <- function(s_r=4, s_c=2, a_c_mean=1, a_c_range=0.2, 
                       a_0_mean=0, a_0_range=0, ...){
  alpha <- matrix(data = runif(s_r*s_c,a_c_mean-a_c_range/2,a_c_mean+a_c_range/2), 
                  ncol=s_c, nrow=s_r) #basic alpha (consumer-resource interactions)
  negdif <- FALSE
  while(!negdif)
  {
    comp  <- matrix(data = runif(s_r*s_r,a_0_mean-a_0_range/2,a_0_mean+a_0_range/2), #basic alpha
                    ncol=s_r, nrow=s_r)
    diag(comp) <- 1
    negdif<- max(eigen(-comp+t(-comp))$values)<0
  }
  zeros <- diag(s_c)*0
  return(rbind(cbind(-comp, -alpha), cbind(t(alpha), zeros)))
}

# Make predator-consumer submatrix
# This is a block matrix, with the upper block consisting of all zeros
# Needed to tell that predators don't eat resources
# and the lower block consisting of randomly (uniformally dist.d) drawn 
# consumption rates (mean = a_p_mean; range = a_p_range)
# s_r, s_c, s_p are nr of resources, consumers, predators
make_A_PC <- function(s_r=4, s_c=2, s_p=1, a_p_mean=1, a_p_range=0.5, ...) {
  zeros_rp <- matrix(data=0, nrow=s_r, ncol=s_p)
  return(-rbind(zeros_rp, matrix(data=runif(s_c*s_p, a_p_mean - a_p_range/2, 
                                     a_p_mean + a_p_range/2), nrow=s_c, ncol=s_p)))
}

# Make whole food chain(/web) matrix by putting together the consumer-resource and PC matrices
# in the right way
make_A_chain <- function(s_p=1, A_CR, A_PC, ...) {
  #Compute overall A for the chain
  zeros_p  <- diag(s_p)*0
  A        <- rbind(cbind(A_CR, A_PC), cbind(t(-A_PC), zeros_p))
  return(A)
}

# Get diagonal matrix containing the f^-1 for a consumer-resource system on its diagonal, 
# Pars are s_r=nr of resource sp; s_c=nr of consumer sp; 
# f_i_mean =mean uniform describing effects 
# always a range of 50% of the mean on these uniforms
get_Finv_CR <- function(s_r=4, s_c=2, f_c_mean=1.2, f_r_mean=0.8, ...) {
  f_c_range <- f_c_mean/2 
  f_r_range <- f_r_mean/2 
  f_r       <- runif(s_r, f_r_mean-f_r_range/2, f_r_mean+f_r_range/2)
  f_c       <- runif(s_c, f_c_mean-f_c_range/2, f_c_mean+f_c_range/2)
  list(Finv_cr=diag(c(f_r, f_c)^-1), 
       f_r_mean_true=mean(f_r),
       f_c_mean_true=mean(f_c))
}
# Same as get_Finv_CR but for predators
get_Finv_P <- function(s_p=1, f_p_mean=1.2, ...) {
  f_p_range <- f_p_mean/2 
  f_p       <- runif(s_p, f_p_mean-f_p_range/2, f_p_mean+f_p_range/2)
  list(Finv_p=diag(f_p^-1, nrow=s_p, ncol=s_p), 
       f_p_mean_true=mean(f_p))
}
#Put f^-1 together for resources, consumers, and predators
get_Finv_RCP <- function(Finv_cr, Finv_p, ...) {
  if(is.null(dim(Finv_p))) {Finv_p <- array(Finv_p, dim=c(1,1))}
  zeros_topright   <- array(0, dim=c(nrow(Finv_cr), ncol(Finv_p)))
  zeros_bottomleft <- array(0, dim=c(nrow(Finv_p), ncol(Finv_cr)))
  return(rbind(cbind(Finv_cr, zeros_topright),
               cbind(zeros_bottomleft, Finv_p)))
}

# Get diagonal matrix containing the effects g on the consumption rates in a consumer-resource system on its diagonal, 
# Pars are s_r=nr of resource sp; s_c=nr of consumer sp; 
# g_mean =mean uniform describing effects 
# always a range of 50% of the mean on these uniforms
get_G_C <- function(s_r=4, s_c=2, g_c_mean=1.2, ...) {
  g_c_range <- g_c_mean/2 
  g_c       <- runif(s_c, g_c_mean-g_c_range/2, g_c_mean+g_c_range/2)
  list(G_C=diag(c(rep(1,s_r),g_c)), g_c_mean_true = mean(g_c))
}

# Get diagonal matrix containing the effects g on the consumption rates of predators on its diagonal, 
# Pars are s_c=nr of consumer sp; s_p=nr of predator sp; 
# g_p_mean =mean uniform describing effects on predator's attack rate
# always a range of 50% of the mean on these uniforms
get_G_P <- function(s_c=1, s_p=1, g_p_mean=1.2, ...) {
  g_p_range <- g_p_mean/2 
  g_p       <- runif(s_p, g_p_mean-g_p_range/2, g_p_mean+g_p_range/2)
  list(G_P=diag(g_p, nrow=s_p, ncol=s_p), g_p_mean_true = mean(g_p))
}
  
#generates p intrinsic growth rate vectors, stacked underneath each other
#for use in spatial LV simulations
#negative_sign is a vector telling which species needs to have a negative intrinsic growth rate
#example: c(2,3) is the case where the intrinsic growth rate of sp 2 and 3 carry a negative sign
make_R_spatial <- function(p, negative_sign, ...){
  Rs     <- as_tibble(t(sample_sphere_surface(3, n = p, radius = 1, positive = TRUE)))%>%
    rename_with(~gsub("V","", .x, fixed = TRUE)) %>% #give correct sign to intrinsic growth rate of consumer
    mutate(location=c(1:p)) %>%
    pivot_longer(!location) %>%
    mutate(value=value*(-1)*name%in%negative_sign + value*(1-name%in%negative_sign))
  return(Rs$value)
}

## Calculation of feasibility (Xi) with equations (no environmental change effects yet) -------------------------------------

# compute Xi for 3 sp food chain
# a is the attack rate (efficiency is 1)
get_Xi_3chain <- function(a_2, a_3) {
  arg <- (a_2*a_3)/((1 + sqrt(1 + a_2^2))*(a_2 + sqrt(a_2^2 + a_3^2)))
  return(4/pi * atan(arg))
}

# compute Xi for a 3sp consumer-resource system
# alpha is the attack rate (efficiency is 1)
get_Xi_3CR <- function(alpha) {
  4/pi*atan((sqrt(2)*alpha*sqrt(1 + alpha^2))/
              (sqrt(2) + sqrt(1 + alpha^2) + alpha^2*(sqrt(2) + 2*sqrt(1 + alpha^2))))
}

## Calculation of feasibility (Xi) w simulations -------------------------------------

#Just wrapper around feasibility() for easier use of pmap.
get_Xi_nCR <- function(A_CR, ...){
  return(feasibility(A_CR))
}

#Just combines making an A for a food chain and computing feasibility for easier use of pmap.
#s_r, s_c, s_p are number of resources, consumers, and predators; 
#A_CR and A_PC are matrices with consumer-resource and predator-consumer int's.
get_Xi_nchain <- function(s_r, s_c, s_p, A_chain, ...){
  #Compute feasibility, only accounting for the right orthant, using Chuliang's function 
  B        <- diag(c(rep(-1, s_r), rep(1, s_c+s_p))) #constraints
  return(2^(s_r+s_c+s_p)*calculate_omega_constraint(A=A_chain, B=B)^(s_r+s_c+s_p))
}

# compute feasibility empirically in a 3sp system, 
# by evaluating species richness at various combinations of r.
# Interpret these various combos as different locations, 
# and allow them to be connected by dispersal.
# A_spatial, R_spatial, D are matrices of sp. ints, intrinsic growth rates, and dispersal
get_Xi_3_spatial <- function(n, A_spatial, R_spatial, D, ...) {
  p                         <- length(R_spatial)/n #infer nr of patches
  initial                   <- -solve(A_spatial)%*%R_spatial#initial guess of equilibrium = eq. w/o dispersal, ...
  initial[which(initial<0)] <-1e-2 #...setting negative values to very small
  local_richness            <- get_richness_spatial(n=3, A=A_spatial, D=D, R=R_spatial) #compute local richness
  return(sum(local_richness$local_richness==3)/p)
}

## Calculation of dXi/de -------------------------------------

#gives the derivative (with respect to epsilon, the env. change factor) of 
#Xi for a 3sp food chain (sp1 is the resource; sp.3 is top predator), 
#where the intrinsic growth rates and consumption rates are affected in a multiplicative way
#by the env. change factor epsilon.
#Pars are:
#a = attack rate (same for both)
#r_i = response of sp_i's intrinsic growth rate (growth for resources, mortality for consumer) to epsilon
#r_2g and r_3g = response of sp 2 and 3's attack rate to epsilon 
get_dXi_3Chain_de <- function(a_2=0.5, a_3=0.5, r_1=-0.1, r_2=-0.1, r_3=0.1, 
                                     r_2g=-0.1, r_3g=-0.1) {
  delta <- r_1 - r_2 + r_2g
  delta2<- r_1 - r_2g - r_3 + r_3g  
  num   <- (2*a_2*a_3)*(sqrt(a_2^2+a_3^2)*delta + a_2*sqrt(1+a_2^2)*delta2)
  denom <- (a_3^2+a_2*(a_2+a_2^3+a_2*a_3^2+sqrt((1+a_2^2)*(a_2^2+a_3^2))))*pi
  return(num/denom)
}

#gives the derivative (with respect to epsilon, the env. change factor) of 
#Xi for a 3sp consumer-resource system (sp1 and 2 are resources; sp.3 is a consumer), 
#where the intrinsic growth rates and consumption rates are affected in a multiplicative way
#by the env. change factor epsilon.
#Pars are:
#a = attack rate
#r_i = response of sp_i's intrinsic growth rate (growth for resources, mortality for consumer) to epsilon
#r_3g = response of sp_3's attack rate to epsilon 
get_dXi_3CR_de <- function(a, r_1=-0.1, r_2=-0.1, r_3=0.1, r_3g=-0.1) {
  return(2*a*(sqrt(2)+2*sqrt(1+a^2)-sqrt(2)*a^2*(1+2*a^2))*(r_1+r_2-2*r_3+2*r_3g)/
    ((1+3*a^2+2*a^4)*(3+2*a^2+2*sqrt(2)*sqrt(1+a^2))*pi))
}

## Related to LV simulations -------------------------------------

#Cut-off function to limit the values of the state to enhance 
#numerical stability when solving the odes
cutoff <- function(n, a) {
  return(1*(n<a)*(3*(n/a)^2-2*(n/a)^3)*ifelse(n<0, 0, 1)+ifelse(n<a, 0, 1))
}

#the LV model; R= a vector of max. intrinsic growth rates; 
#A= matrix with linear species interactions, global var 
run_LV <- function(time, state, pars) {
  return(list(as.numeric(
    ((state*(pars$R+pars$A%*%state)))*cutoff(state, 1e-10)
  )))
}

#Spatial LV model; R= a vector of max. intrinsic growth rates; 
#A= matrix with linear species interactions, D=dispersal matrix
run_LV_spatial <- function(time, state, pars) {
  return(list(as.numeric(
    ((state*((pars$R-colSums(pars$D))+pars$A%*%state))+pars$D%*%state)*
      cutoff(state, 1e-10)
  )))
}

#compute richness from spatial LV simulation
#n = nr of sp
#R, A, D: as in run_LV_spatial
#max_time = max simulation time
get_richness_spatial <- function(n, R, A, D, max_time=100, 
                         extinction_threshold=1e-3){
  p         <- length(R)/n #infer nr of patches
  densities <- ode(func=run_LV_spatial, y=rep(0.1,n*p), 
                   parms=list(R=R, A=A, D=D), 
                   times=seq(1,max_time,0.1))
  locations <- NULL; for (i in 1:p) {locations <- c(locations, rep(i, n))}
  return(as_tibble(cbind(densities=tail(densities, 1)[-1], locations)) %>%
    group_by(locations) %>%
    summarize(local_richness=sum(densities>extinction_threshold)))
}

## Misc. ---------------

# normalise a vector to unit length
normalize_vector <- function(g) {
  g/sqrt(sum(g^2))
}


  



