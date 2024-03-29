library(momentuHMM)
library(plot3D)

SVtracks <- read.csv("../data/trips_SV_200_steps.csv", stringsAsFactors = FALSE)

## select few tracks
colnames(SVtracks) <- c("ID", "lon", "lat")
select = unique(SVtracks$ID)
tracks <- SVtracks[!is.na(match(SVtracks$ID, select)),]

## PROJECTION geostrophic
colony <- c(-77.264,-11.773)
R <- 6371

tracks$x <- pi * R * (tracks$lon - colony[1])/180
tracks$y <- pi * R * (tracks$lat - colony[2])/180

center <- matrix(c(0,0),nrow=1,dimnames=list("colony"))
data <- prepData(data = tracks, type = 'UTM', centers = center)


### Model 1 : no covariates
stateNames <- c("search","forage","inbound")

# initial parameters
stepPar0 <- c(0.7, 0.3, 1,
              0.1,0.2,0.05)
anglePar0 <- c(11,0.4,15)

# constrain transition probabilities
fixbeta <- matrix(c(NA,-100,NA,NA,-100,-100),
                  nrow=1)

m1 <- fitHMM(data=data, nbStates=3, dist=list(step="gamma",angle="vm"),
             Par0=list(step=stepPar0, angle=anglePar0),
             fixPar=list(beta=fixbeta), stateNames = stateNames)

# plot(m1)


### Model 2 : Biased random walks are used to model the movement in states 1 and 3, with
# repulsion/attraction towards the colony

angleFormula <- ~ state3(colony.angle)
fixPar <- list(angle=c(1,NA,NA,NA),beta=fixbeta)

Par0 <- getPar0(model=m1, nbStates=3,
                DM=list(angle=list(mean=angleFormula, concentration=~1)),
                estAngleMean=list(angle=TRUE),
                circularAngleMean=list(angle=0))


m2 <- fitHMM(data=data, nbStates=3, dist=list(step="gamma",angle="vm"),
             Par0=list(step=Par0$Par$step, angle=Par0$Par$angle),
             beta0=Par0$beta, fixPar=fixPar,
             DM=list(angle=list(mean=angleFormula, concentration=~1)),
             estAngleMean=list(angle=TRUE), circularAngleMean=list(angle=0),
             stateNames = stateNames)

# plot(m2)

### Model 3 : Biased random walks are used to model the movement in states 1 and 3, with
# repulsion/attraction towards the colony + transition probability effect

# time spent since left colony
time <- NULL
for(id in unique(data$ID)) {
  nbSubObs <- length(which(data$ID==id))
  # approximately in months for interval = 9.6h
  time <- c(time, (1:nbSubObs)/60)
}
data$time <- time

# compute time since departure and include in formula below
formula <- ~ time

fixbeta <- matrix(c(NA,-100,NA,NA,-100,-100,
                    0, 0, 0, NA, 0, 0),
                  nrow=2,byrow=TRUE)

# angleFormula <- ~ state3(colony.angle)
# fixPar <- list(angle=c(1,NA,NA,NA),beta=fixbeta)

Par0 <- getPar0(model=m1, nbStates=3,
                # DM=list(angle=list(mean=angleFormula, concentration=~1)),
                # estAngleMean=list(angle=TRUE),
                circularAngleMean=list(angle=0),
                formula=formula)

m3 <- fitHMM(data=data, nbStates=3, dist=list(step="gamma",angle="vm"),
             Par0=list(step=Par0$Par$step, angle=Par0$Par$angle),
             beta0=Par0$beta, 
             # fixPar=fixPar,
             # DM=list(angle=list(mean=angleFormula, concentration=~1)),
             # estAngleMean=list(angle=TRUE),
             circularAngleMean=list(angle=0),
             stateNames = stateNames, 
             formula=formula)


# plot(m3)

### Model 4 : both transition prob + biaised RW

Par0 <- getPar0(model=m3, nbStates=3,
                DM=list(angle=list(mean=angleFormula, concentration=~1)),
                estAngleMean=list(angle=TRUE),
                circularAngleMean=list(angle=0),
                formula=formula)

fixPar <- list(angle=c(1,NA,NA,NA), beta=Par0$beta)

m4 <- fitHMM(data=data, nbStates=3, dist=list(step="gamma",angle="vm"),
             Par0=list(step=Par0$Par$step, angle=Par0$Par$angle),
             beta0=Par0$beta, 
             fixPar=fixPar,
             DM=list(angle=list(mean=angleFormula, concentration=~1)),
             estAngleMean=list(angle=TRUE),
             circularAngleMean=list(angle=0),
             stateNames = stateNames, 
             formula=formula)


# plot(m4)

viterbi(m4)

AIC(m1,m2,m3,m4)

### LOAD GEBCO  

library(ncdf4)

getBathy <- function(newdata){
  
  gebco <- nc_open("/home/amdroy/MEGA/DATA/SATELLITE/GEBCO/gebco_2019_pescadores.nc")
  lon <- ncvar_get(gebco, 'lon')
  lat <- ncvar_get(gebco, 'lat')
  
  x <- pi * R * (lon - colony[1])/180
  y <- pi * R * (lat - colony[2])/180
  
  X <- newdata$x
  Y <- newdata$y
  
  elevation <- sapply(1:nrow(newdata), function(i){
    idx_x = which.min(abs(x-X[i]))
    idx_y = which.min(abs(y-Y[i]))
    
    elev <- ncvar_get(gebco, 'elevation', start = c(idx_x, idx_y), count = c(1,1))
  })
  
  return(elevation)
}


### SIMULATE TRAJECTORY

## initial direction
# id = data$ID[1]
theta_all <- NULL
for (id in unique(data$ID)){
  tt <- data[data$ID == id,]
  
  theta <- atan2((tt$y[2] - tt$y[1]),(tt$x[2] - tt$x[1]))%%2*pi
  theta_all <- c(theta_all, theta)
}


gebco <- nc_open("/home/amdroy/MEGA/DATA/SATELLITE/GEBCO/gebco_2019_pescadores.nc")
lon <- ncvar_get(gebco, 'lon')
lat <- ncvar_get(gebco, 'lat')

TRAJ <- NULL

kk <- 1
while (kk <= 100){

theta <-  sample(theta_all, 1)
# theta <- runif(1, pi/2, 3*pi/2)

Par <- getPar0(m4)
newdata <- simData(nbAnimals = 1, 
                   nbStates = 3,
                   dist = list(step = 'gamma', angle = 'vm'),
                   Par = Par$Par,
                   beta = Par$beta,
                   DM = list(angle = list(mean = ~state3(colony.angle), concentration = ~1)),
                   delta = c(1-1e-10, 1e-10/2, 1e-10/2),
                   formula = m4$conditions$formula,
                   zeroInflation = m4$conditions$zeroInflation,
                   oneInflation = m4$conditions$oneInflation,
                   circularAngleMean = m4$conditions$circularAngleMean,
                   covs = data.frame(time = seq(1/60, 2, by = 1/60)),
                   initialPosition = c(0, 0),
                   centers = matrix(c(0,0),nrow=1,dimnames=list("colony")),
                   stateNames = stateNames,
                   states = TRUE,
                   obsPerAnimal = 200)

X <- newdata$x * cos(theta) - newdata$y * sin(theta)
Y <- newdata$x * sin(theta) + newdata$y * cos(theta)

traj <- newdata
traj$x <- X
traj$y <- Y

nb <- 0
while(sum(getBathy(traj) > 0) > 0){
  
  ii <- which(getBathy(traj) > 0)[1]
  traj <- traj[1:(ii-1),]
  
  a <- as.double(traj[ii-2, c("x", "y")])
  b <- as.double(traj[ii-1, c("x", "y")])
  theta <- atan2((b[2] - a[2]),(b[1] - a[1]))%%2*pi
  
  # initial transformation
  initial_position <- c(0,0)
  initial_position[1] <- traj$x[ii-1] * cos(theta) + traj$y[ii-1] * sin(theta)
  initial_position[2] <- traj$x[ii-1] * -sin(theta) + traj$y[ii-1] * cos(theta)
            
  newdata <- simData(nbAnimals = 1, 
                      nbStates = 3,
                      dist = list(step = 'gamma', angle = 'vm'),
                      Par = Par$Par,
                      beta = Par$beta,
                      DM = list(angle = list(mean = ~state3(colony.angle), concentration = ~1)),
                      delta = c(1-1e-10, 1e-10/2, 1e-10/2),
                      formula = m4$conditions$formula,
                      zeroInflation = m4$conditions$zeroInflation,
                      oneInflation = m4$conditions$oneInflation,
                      circularAngleMean = m4$conditions$circularAngleMean,
                      covs = data.frame(time = seq(1/60*ii, 5, by = 1/60)),
                      initialPosition = initial_position,
                      centers = matrix(c(0,0),nrow=1,dimnames=list("colony")),
                      stateNames = stateNames,
                      states = TRUE,
                      obsPerAnimal = 200-ii+2)
  
  X <- newdata$x * cos(theta) - newdata$y * sin(theta)
  Y <- newdata$x * sin(theta) + newdata$y * cos(theta)
  
  newtraj <- newdata
  newtraj$x <- X
  newtraj$y <- Y
  
  traj <- rbind(traj, newtraj[2:nrow(newtraj),])
  
  nb <- nb + 1
  if ( nb > 10){
    break
  }
}

if ( nb < 10){
  ### PLOT TRAJECTORY
  x <- pi * R * (lon - colony[1])/180
  y <- pi * R * (lat - colony[2])/180
  bathy <-  ncvar_get(gebco, "elevation")
  bathy[bathy < 0] = 0
  bathy[bathy > 0] = 1
  
  
  image2D(bathy, x, y, xlim = c(-50, 50), ylim = c(-50, 50), main = theta)
  points(traj$x, traj$y, col = traj$states + 1)
  
  traj$ID <- kk
  TRAJ <- rbind(TRAJ, traj)
  kk <- kk + 1
}

}

write.csv(TRAJ, file = '../results/sim_HMM_SV.csv')
