#setwd("~/Desktop/GitHub/lstrachan_honeybee_sim/YearCycleSimulation")
# Clean workspace

rm(list = ls())
getwd()


# Define functions
maintainIrelandSize <- function(age0 = NULL, age1 = NULL) {
  if ((nColonies(age0) + nColonies(age1)) > IrelandSize) { # check if the sum of all colonies is greater than population size
    IDsplits <- getId(age0)[hasSplit(age0)] # get the IDs of age 0 that are splits
    splits0 <- pullColonies(age0, ID = IDsplits) # pull the splits out of age 0
    age0split <- splits0$pulled # create an object for age 0 splits
    age0swarm <- splits0$remnant # create an object for swarms and superseded colonies
    splitsI<-pullColonies(age0split, ID=IdImportColonies) #pull the imports out of split
    age0splitImport <- splitsI$pulled # create an object for imported splits
    age0splitMel <- splitsI$remnant # create an object for non imported splits
    age0needed <- IrelandSize - nColonies(age1) # calculate the number of age 0 colonies that are needed to fill up the apiary
    if (age0needed <= nColonies(age0swarm)) { # check if the number of age 0 colonies needed is lower or equal to age 0 swarms
      swarmID <- sample(getId(age0swarm), age0needed) # if yes, select the ids of swarms that will stay in apiary
      swarmTMP <- pullColonies(age0swarm, ID = swarmID) # pull out those selected age0 swarms
      age0 <- swarmTMP$pulled # put selected swarms to age 0 object
    } else if (age0needed > nColonies(age0swarm)) { # in case when age 0 needed is greater than number of importedsplits select splits
      nSplitNeeded <- age0needed - nColonies(age0swarm) # calculate the number of splits needed
      if (nSplitNeeded>nColonies(age0splitImport)){ #check if the number of importedsplits needed is greater than the number of importedsplits we have
        age0<- c(age0swarm,age0splitImport) #add all the imported splits to age 0 object
        nSplitMelNeeded <- age0needed - nColonies(age0) # calculate the number of not imported splits needed
        splitMelId <- sample(getId(age0splitMel), nSplitMelNeeded) # select ids of not imported split 
        splitMelTmp <- pullColonies(age0splitMel, ID = splitMelId) # pull the splits
        splitsMel<- splitMelTmp$pulled # select pulled splits
        age0<-c(age0,splitsMel) #add them to age0 object
      } else {  #if the imported splits needed are lower than the imported splits we have
        splitImportId <- sample(getId(age0splitImport), nSplitNeeded) # select ids of split import
        ImportTmp <- pullColonies(age0splitImport, ID = splitImportId) # pull the split import
        splitImports<- ImportTmp$pulled # select pulled split import
        age0 <- c(age0swarms, splitImports) # combine imported splits and swarms in age 0 object
      }
    }
    return(age0)
  }
}
maintainCarSize <- function(age0 = NULL, age1 = NULL) {
  if ((nColonies(age0) + nColonies(age1)) > CarSize) { # check if the sum of all colonies is greater than apiary size
    IDsplits <- getId(age0)[hasSplit(age0)] # get the IDs of age 0 that are splits
    splits0 <- pullColonies(age0, ID = IDsplits) # pull the splits out of age 0
    age0split <- splits0$pulled # create an object for age 0 splits
    age0swarm <- splits0$remnant # create an object for swarms and superseded colonies
    age0needed <- CarSize - nColonies(age1) # calculate the number of age 0 colonies that are needed to fill up the apiary
    splitsNeeded <- age0needed - nColonies(age0swarm) # calculate the number of splits needed
    if (age0needed <= nColonies(age0swarm)) { # check if the number of age 0 colonies needed is lower or equal to age 0 swarms
      swarmID <- sample(getId(age0swarm), age0needed) # if yes, select the ids of swarms that will stay in apiary
      swarmTMP <- pullColonies(age0swarm, ID = swarmID) # pull out those selected age0 swarms
      age0 <- swarmTMP$pulled # put selected swarms to age 0 object
    } else if (age0needed > nColonies(age0swarm)) { # in case when age 0 needed is grater than number of swarm select splits
      nSplitsNeeded <- age0needed - nColonies(age0swarm) # calculate the number of splits needed
      splitId <- sample(getId(age0split), nSplitsNeeded) # select ids of splits
      splitTmp <- pullColonies(age0split, ID = splitId) # pull the splits
      splits <- splitTmp$pulled # select pulled splits
      age0 <- c(age0swarm, splits) # combine splits and swarms in age 0 object
    }
    return(age0)
  }
}


# Load packages
library(AlphaSimR)
library(ggplot2)
library(tictoc)
library(R6)
library(nadiv)
library(Matrix)
library(SIMplyBee)
library(dplyr)
library(tidyr)
# TODO: replace with devtools installation from Github once the package is operational
# Source the development version of AlphaSimR
getwd()

# Founder population parameters -------------------------------------------------------------------
nMelN = 450                   # Number of Mellifera 450
nCar = 150                    # Number of Carnica 150
nChr = 1                     # Number of chromomsome
nDronesPerQueen = 50
nSegSites = 100              # Number of segregating sites

# Population parameters -------------------------------------------------------------------
nRep <- 1                     # Number of repeats
nYear <- 10                   # Number of years
IrelandSize<-300              #Ireland population size 300
CarSize<-100                  #Carnica pop size 100
nWorkers <- 10                # Number of workers in a full colony
nDrones <- 50                 # Number of drones in a full colony (typically nWorkers * 0.2 (not in the example))
pFathers <- nFathersPoisson   # Number of drones the queen mates with (could also be a function)
nVirginQueens <- 1            # Number of created virgin queens

# Period parameters -------------------------------------------------------------------
# Period 1 (spring)
p1swarm <- 0.05              # Percentage of colonies that swarm in period 1
p1supersede <- 0.05          # Percentage of colonies that supersede in period 1
p1collapse <- 0.10           # Percentage of colonies that  collapse in period 1

# Period2 (summer)
p2swarm <- 0.01              # Percentage of colonies that swarm in period 2
p2supersede <- p1supersede   # Percentage of colonies that supersede in period 2
p2collapse <- p1collapse     # Percentage of colonies that collapse in period 2

# Period3 (winter)
p3collapseAge0 <- 0.25      # Percentage of age 0 colonies that collapse in period 3
p3collapseAge1 <- 0.3       # Percentage of age 2 colonies that collapse in period 3

#Import parameters -------------------------------------------------------------------
pImport <- 0.3              # Percentage import from carnica to mellifera

# Create data frames for recording the number of age0 and age1 colonies, csd variability and for recording cpu time
loopTime <- data.frame(Rep = NA, tic = NA, toc = NA, msg = NA, time = NA)

# Prepare recording function
data_rec <- function(datafile, colonies, year, population) {
  queens = mergePops(getQueen(colonies))
  IBDh = apply(getIbdHaplo(queens),MARGIN = 1, FUN =  function(X) sum(X %in% 1:(nMelN*2)/length(X)))
  datafile = rbind(datafile,
                   data.frame(colonies             = deparse(substitute(colonies)),
                              population           = population,
                              year                 = year,
                              Id                   = queens@id,
                              MId                  = queens@mother,
                              FId                  = queens@father,
                              nFathers             = nFathers(queens),
                              nDPQ                 = sapply(getFathers(queens), function(x) length(unique(x@mother))),
                              nCsdAlColony         = sapply(colonies@colonies, function(x) nCsdAlleles(x, collapse = TRUE)),
                              nCsdApiary           = rep(nCsdAlleles(colonies, collapse = TRUE), queens@nInd),
                              pHomBrood            = calcQueensPHomBrood(queens),
                              gvQueens_QueenHoneyYield  = sapply(getGv(colonies, caste = "queen"), function(x) x[1,1]),
                              gvQueens_QueenFitness = sapply(getGv(colonies, caste = "queen"), function(x) x[1,2]),
                              IBD = sapply(seq(1,length(IBDh),2), FUN = function(z) sum(IBDh[z:(z+1)])/2)
                              
                   ))}
colonyRecords = NULL

# Start of the rep-loop ---------------------------------------------------------------------
#for (Rep in 1:nRep) {
Rep <- 1 #(you can use this to check if your code is running alright for one whole loop)
cat(paste0("Rep: ", Rep, "/", nRep, "\n"))
tic(paste0(nYear, 'y loop'))         # Measure cpu time
Rprof()                              # Start profiling


# Founder population ---------------------------------------------------------
# STEP 1:  Create a founder population of A. m. mellifera and A. m. carnica bees (un-# the one you want to use)

#load("C:/Users/Usuario/Desktop/Máster/Bees/FounderGenomes_ThreePop_16chr.RData")

# quick haplo to get the founder genomes for now.
founderGenomes<- quickHaplo(sum(nMelN,nCar),4,segSites = 1000)

# STEP 2: Create SP object and write in the global simulation/population parameters
SP <- SimParamBee$new(founderGenomes, csdChr = ifelse(nChr >= 3, 3, 1), nCsdAlleles = 128)
SP$nWorkers <- nWorkers
SP$nDrones <- nDrones
SP$nFathers <- pFathers
SP$nVirginQueens <- nVirginQueens
SP$swarmP <- 0.5                # Probability of swarming? (TODO: double check this)
SP$splitP <- 0.3
SP$setTrackPed(TRUE)            # Track the pedigree
SP$setTrackRec(TRUE)            # Track the recombination

SP$addSnpChip(nSnpPerChr = 10)   # Add a SNP chip with 3 SNPs per chromosome
csdChr <- SP$csdChr             # define csd chromomsome

# Add traits - taken from the QuantGen vignette 
mean <- c(0, 0)
varA <- c(0.25, 0.1)
corA <- matrix(data = c( 1.0, 0,
                         0,  1.0), nrow = 2, byrow = TRUE)
SP$addTraitA(nQtlPerChr = 100, mean = mean, var = varA, corA = corA,
             name = c("QueenHoneyYield", "QueenFitness"))

varE <- c(0.75, 0.9 / SP$nWorkers)

# TODO: what is a reasonable environmental correlation between queen and worker effects?
corE <- matrix(data = c(1.0, 0,
                        0, 1.0), nrow = 2, byrow = TRUE)
SP$setVarE(varE = varE, corE = corE)

# STEP 3: Set up your base population
# Create a base population for A. m. mellifera, A. m. mellifera cross, and A. m. carnica (400 of each)
virginQueens <- list(Mel = createVirginQueens(x = founderGenomes[1:(nMelN)]),
                     Car = createVirginQueens(x = founderGenomes[(nMelN+1):(nMelN + nCar)]))

# Create drones for A. m. mellifera, A. m. mellifera cross, and A. m. carnica
drones <- list(Mel = createDrones(x = virginQueens$Mel[(IrelandSize+1):(nMelN)], nInd = nDronesPerQueen),
               Car = createDrones(x = virginQueens$Car[(CarSize+1):nCar], nInd = nDronesPerQueen))


# Get fathers for Mel, MelCross and Car
fathersMel <- pullDroneGroupsFromDCA(drones$Mel, n = nInd(virginQueens$Mel[1:IrelandSize]), nDrones = nFathersPoisson) #crea 300 grupos de drones con 10 drones por grupo
fathersCar <- pullDroneGroupsFromDCA(drones$Car, n = nInd(virginQueens$Car[1:CarSize]), nDrones = nFathersPoisson)

# Mate virgin queens with fathers to make them queens
queens <- list(Mel = SIMplyBee::cross(x = virginQueens$Mel[1:IrelandSize], drones = fathersMel),
               Car = SIMplyBee::cross(x = virginQueens$Car[1:CarSize], drones = fathersCar))

year=1

# Start the year-loop ------------------------------------------------------------------
#for (year in 1:nYear) {
print("Starting the cycle")
#year <- 1 (Use this to check that things are working without setting the whole for loop off )
#year <- year + 1
cat(paste0("Year: ", year, "/", nYear, "\n"))

# If this is the first year, create some colonies to start with
if (year == 1) {
  print("Creating initial colonies")
  age1 <- list(Mel = createMultiColony(x = queens$Mel, n = IrelandSize),
               Car = createMultiColony(x = queens$Car, n = CarSize))
  #aqui empezamos a hacer lo de spatial
  
  #age1$Mel <-  setLocation(age1$Mel, 
  #                      location = Map(c, runif(nColonies(age1$Mel), 0, 6), runif(nColonies(age1$Mel), 0, 6)))
  nCol<-nColonies(age1$Mel)
  age1i<-pullColonies(age1$Mel, n= nCol/2)
  age1dub<-age1i$pulled
  age1rest<-age1i$remnant
  
  age1dub<-setLocation(age1dub, location = Map(c, runif(nColonies(age1dub), 3, 6), runif(nColonies(age1dub), 0, 3)))
  
  nColo<-nColonies(age1rest)
  age1o<-pullColonies(age1rest, n= nColo/3)
  haf<-age1o$pulled
  rest<-age1o$remnant
  
  haf<-setLocation(haf, location = Map(c, runif(nColonies(haf), 0, 3), runif(nColonies(haf), 0, 3)))
  rest<-setLocation(rest, location = Map(c, runif(nColonies(rest), 0, 6), runif(nColonies(haf), 3, 6)))
  locationsDF <- data.frame(Location = getLocation(c(age1dub, haf, rest), collapse = TRUE),
                            Beekeeper = c(rep("Beekeeper1", nColonies(age1dub)),
                                          rep("Beekeeper2", nColonies(haf)),
                                          rep("Beekeeper3", nColonies(rest))))
  
  age1$Mel<-c(age1dub,haf,rest)
  age1
  locationsDF <- data.frame(Location = getLocation(c(age1$Mel), collapse = TRUE),
                            Beekeeper = c(rep("Beekeeper1", nColonies(age1$Mel))))
  ggplot(data = locationsDF, aes(x = Location.1, y = Location.2, colour = Beekeeper)) + 
    geom_point()
  
  
  # If not, promote the age0 to age1, age1 to age2 and remove age2 colonies
} else {
  age2 <- list(Mel = age1$Mel, Car = age1$Car)
  age1 <- list(Mel = age0$Mel, Car = age0$Car)
  age0 <- list(Mel = NULL, Car = NULL) 
  age0p1 <- list(Mel = NULL, Car = NULL)
  age0p2 <- list(Mel = NULL, Car = NULL) 
}

# Period1 ------------------------------------------------------------------
# Build-up the colonies
print(paste0("Building up the colonies to ", nWorkers, " and ", nDrones))
print(Sys.time())
age1 <- list(Mel = buildUp(age1$Mel),
             Car = buildUp(age1$Car))
if (year > 1) {
  age2 <- list(Mel = buildUp(age2$Mel),
               Car = buildUp(age2$Car))
}

# Split all age1 colonies
print("Splitting the colonies")
print(Sys.time())
tmp <- list(Mel = split(age1$Mel),
            Car = split(age1$Car))
age1 <- list(Mel = tmp$Mel$remnant,
             Car = tmp$Car$remnant)

# The queens of the splits are 0 years old
age0p1 <- list(Mel = tmp$Mel$split, Car = tmp$Car$split) 

#Set the location of splits to a location near where the original colony is
x <- sapply(getLocation(age1$Mel), function(X) X[1])
y <- sapply(getLocation(age1$Mel), function(X) X[2])
age0p1$Mel<-setLocation(age0p1$Mel, location = Map(c, (x+runif(length(x), min = -0.1, max = 0.1)), (y+runif(length(y), min = -0.1, max = 0.1))))

if (year > 1) {
  # Split all age2 colonies
  tmp <- list(Mel = split(age2$Mel),
              Car = split(age2$Car))
  
  age2 <- list(Mel = tmp$Mel$remnant,
               Car = tmp$Car$remnant)
  
  #Set the location of splits to a location near where the original colony is
  x <- sapply(getLocation(age2$Mel), function(X) X[1])
  y <- sapply(getLocation(age2$Mel), function(X) X[2])
  tmp$Mel$split<-setLocation(tmp$Mel$split, location = Map(c, (x+runif(length(x), min = -0.1, max = 0.1)), (y+runif(length(y), min = -0.1, max = 0.1))))
  # The queens of the splits are 0 years old
  age0p1 <- list(Mel = c(age0p1$Mel, tmp$Mel$split),
                 Car = c(age0p1$Car, tmp$Car$split))
}


print("Create virgin queens, period 1")
print(Sys.time())
virginDonor <- list(Mel = sample.int(n = nColonies(age1$Mel), size = 1),
                    Car = sample.int(n = nColonies(age1$Car), size = 1))
# Virgin queens for splits!
pImport<-0.3
all <- sample(getId(age0p1$Mel)[sapply(getLocation(age0p1$Mel), function(coords) coords[1] >= 3 & coords[2] <= 3)],nColonies(age0p1$Mel)*pImport)
idstopull<-all
tmp <- (Mel = pullColonies(age0p1$Mel,ID=idstopull))
IdImportColonies<-getId(tmp$pulled)
age0p1 <- list(Mel = tmp$remnant,
               MelImport = tmp$pulled,
               Car = c(age0p1$Car, tmp$Car$split))

virginQueens <- list(Mel = createVirginQueens(age1$Mel[[virginDonor$Mel]], nInd = nColonies(age0p1$Mel)),
                     Car = createVirginQueens(age1$Car[[virginDonor$Car]], nInd = nColonies(age0p1$Car)+nColonies(age0p1$MelImport)))


# Requeen the splits --> queens are now 0 years old

nColoniesMelImport<-nColonies(age0p1$MelImport)
nColoniesCar<-nColonies(age0p1$Car)+nColonies(age0p1$MelImport)

age0p1 <- list(Mel = c(reQueen(age0p1$Mel, queen = (virginQueens$Mel)) ,
                       reQueen(age0p1$MelImport, queen = c((virginQueens$Car)[1:nColoniesMelImport]))),       #,(virginQueens$Lig)[1:(nColoniesMelImport/2)]))),
               Car = reQueen(age0p1$Car, queen = virginQueens$Car[(nColoniesMelImport+1):nColoniesCar]))

# Swarm a percentage of age1 colonies
print("Swarm colonies, P1")
print(Sys.time())
tmp <- list(Mel = pullColonies(age1$Mel, p = p1swarm),
            Car = pullColonies(age1$Car, p = p1swarm))

age1 <- list(Mel = tmp$Mel$remnant,
             Car = tmp$Car$remnant)

tmp <- list(Mel = swarm(tmp$Mel$pulled, sampleLocation = T, radius = 0.3),
            Car = swarm(tmp$Car$pulled, sampleLocation = T, radius = 0.3))

age0p1 <- list(Mel = c(age0p1$Mel, tmp$Mel$remnant),
               Car = c(age0p1$Car, tmp$Car$remnant))

age1 <- list(Mel = c(age1$Mel, tmp$Mel$swarm),
             Car = c(age1$Car, tmp$Car$swarm))


if (year > 1) {
  # Swarm a percentage of age2 colonies
  tmp <- list(Mel = pullColonies(age2$Mel, p = p1swarm),
              Car = pullColonies(age2$Car, p = p1swarm))
  
  age2 <- list(Mel = tmp$Mel$remnant,     
               Car = tmp$Car$remnant)
  
  tmp <- list(Mel = swarm(tmp$Mel$pulled, sampleLocation = T, radius = 0.3),
              Car = swarm(tmp$Car$pulled, sampleLocation = T, radius = 0.3))
  
  age0p1 <- list(Mel = c(age0p1$Mel, tmp$Mel$remnant),
                 Car = c(age0p1$Car, tmp$Car$remnant))
  
  age2 <- list(Mel = c(age2$Mel, tmp$Mel$swarm),
               Car = c(age2$Car, tmp$Car$swarm))
  
}

# Supersede age1 colonies
print("Supersede colonies, P1")
print(Sys.time())
tmp <- list(Mel = pullColonies(age1$Mel, p = p1supersede),
            Car = pullColonies(age1$Car, p = p1supersede))

age1 <- list(Mel = tmp$Mel$remnant,
             Car = tmp$Car$remnant)

tmp <- list(Mel = supersede(tmp$Mel$pulled),
            Car = supersede(tmp$Car$pulled))

age0p1 <- list(Mel = c(age0p1$Mel, tmp$Mel),
               Car = c(age0p1$Car, tmp$Car))

if (year > 1) {
  # Supersede age2 colonies
  tmp <- list(Mel = pullColonies(age2$Mel, p = p1supersede),
              Car = pullColonies(age2$Car, p = p1supersede))
  
  age2 <- list(Mel = tmp$Mel$remnant,
               Car = tmp$Car$remnant)
  
  tmp <- list(Mel = supersede(tmp$Mel$pulled),
              Car = supersede(tmp$Car$pulled))
  
  age0p1 <- list(Mel = c(age0p1$Mel, tmp$Mel),
                 Car = c(age0p1$Car, tmp$Car))
}

# Mate the split colonies

print("Mate split colonies, P1")
print(Sys.time())
if (year == 1) {
  age0p1$Mel <- cross(age0p1$Mel, droneColonies = age1$Mel, crossPlan= "create", spatial= T, radius= 2, nDrones= nDronesPoisson, checkCross = "warning")
  DCACar <- createDCA(age1$Car)
  age0p1$Car <- cross(age0p1$Car, drones = pullDroneGroupsFromDCA(DCA = DCACar, n = nColonies(age0p1$Car), nDrones = nFathersPoisson))
  
} else {
  age0p1$Mel <- cross(age0p1$Mel, droneColonies = c(age1$Mel,age2$Mel), crossPlan= "create", spatial= T, radius= 2, nDrones= nDronesPoisson, checkCross = "warning")
  DCACar <- createDCA(c(age1$Car, age2$Car))
  age0p1$Car <- cross(age0p1$Car, drones = pullDroneGroupsFromDCA(DCA = DCACar, n = nColonies(age0p1$Car), nDrones = nFathersPoisson))
  
}


# Collapse
print("Collapse colonies, P1")
print(Sys.time())
age1 <- list(Mel = selectColonies(age1$Mel, p = 1 - p1collapse),
             Car = selectColonies(age1$Car, p = 1 - p1collapse))

if (year > 1) {
  age2 <- list(Mel = selectColonies(age2$Mel, p = 1 - p1collapse),
               Car = selectColonies(age2$Car, p = 1 - p1collapse))
  
}


# Period2 ------------------------------------------------------------------
print("PERIOD 2")
# Swarm a percentage of age1 colonies
# Mellifera
print("Swarm colonies, P2")
print(Sys.time())
tmp <- list(Mel = pullColonies(age1$Mel, p = p2swarm),
            Car = pullColonies(age1$Car, p = p2swarm))

age1 <- list(Mel = tmp$Mel$remnant,
             Car = tmp$Car$remnant)

tmp <- list(Mel = swarm(tmp$Mel$pulled, sampleLocation = T, radius = 0.3),
            Car = swarm(tmp$Car$pulled, sampleLocation = T, radius = 0.3))

# The queens of the remnant colonies are of age 0
age0p2 <- list(Mel = tmp$Mel$remnant,
               Car = tmp$Car$remnant)

age1 <- list(Mel = c(age1$Mel, tmp$Mel$swarm),
             Car = c(age1$Car, tmp$Car$swarm))



if (year > 1) {
  # Swarm a percentage of age2 colonies
  tmp <- list(Mel = pullColonies(age2$Mel, p = p2swarm),
              Car = pullColonies(age2$Car, p = p2swarm))
  
  age2 <- list(Mel = tmp$Mel$remnant,
               Car = tmp$Car$remnant)
  
  tmp <- list(Mel = swarm(tmp$Mel$pulled, sampleLocation = T, radius = 0.3),
              Car = swarm(tmp$Car$pulled, sampleLocation = T, radius = 0.3))
  
  # The queens of the remnant colonies are of age 0
  age0p2 <- list(Mel = tmp$Mel$remnant,
                 Car = tmp$Car$remnant)
  
  age2 <- list(Mel = c(age2$Mel, tmp$Mel$swarm),
               Car = c(age2$Car, tmp$Car$swarm))
  
}


# Supersede a part of age1 colonies
print("Supersede colonies, P2")
print(Sys.time())

tmp <- list(Mel = pullColonies(age1$Mel, p = p2supersede),
            Car = pullColonies(age1$Car, p = p2supersede))

age1 <- list(Mel = tmp$Mel$remnant,
             Car = tmp$Car$remnant)

tmp <- list(Mel = supersede(tmp$Mel$pulled),
            Car = supersede(tmp$Car$pulled))

# The queens of superseded colonies are of age 0
age0p2 <- list(Mel = c(age0p2$Mel, tmp$Mel),
               Car = c(age0p2$Car, tmp$Car))

if (year > 1) {
  # Supersede a part of age2 colonies
  tmp <- list(Mel = pullColonies(age2$Mel, p = p2supersede),
              Car = pullColonies(age2$Car, p = p2supersede))
  
  age2 <- list(Mel = tmp$Mel$remnant,
               Car = tmp$Car$remnant)
  
  tmp <- list(Mel = supersede(tmp$Mel$pulled),
              Car = supersede(tmp$Car$pulled))
  
  # The queens of superseded colonies are of age 0
  age0p2 <- list(Mel = c(age0p2$Mel, tmp$Mel),
                 Car = c(age0p2$Car, tmp$Car))
  
}

# Replace all the drones
print("Replace Drones, P2")
print(Sys.time())

age1$Mel <- replaceDrones(age1$Mel)
age1$Car <- replaceDrones(age1$Car)

if (year > 1) {
  age2$Mel <- replaceDrones(age2$Mel)
  age2$Car <- replaceDrones(age2$Car)
  
}

# Mate the colonies
# Import p percentage of carnica colonies into mellifera DCA
print("Mate colonies, P2")
print(Sys.time())


if (year == 1) {
  age0p2$Mel <- cross(age0p2$Mel, droneColonies = age1$Mel, crossPlan= "create", spatial= T, radius= 5, nDrones= nDronesPoisson, checkCross = "warning")
  DCACar <- createDCA(age1$Car)
  age0p2$Car <- cross(age0p2$Car, drones = pullDroneGroupsFromDCA(DCA = DCACar, n = nColonies(age0p2$Car), nDrones = nFathersPoisson))
  
} else {
  age0p2$Mel <- cross(age0p2$Mel, droneColonies = c(age1$Mel,age2$Mel), crossPlan= "create", spatial= T, radius= 5, nDrones= nDronesPoisson, checkCross = "warning")
  DCACar <- createDCA(c(age1$Car, age2$Car))
  fathersCar <-  pullDroneGroupsFromDCA(DCA = DCACar, n = nColonies(age0p2$Car), nDrones = nFathersPoisson)
  age0p2$Car <- cross(age0p2$Car, drones = fathersCar)
}

# Collapse
age1 <- list(Mel = selectColonies(age1$Mel, p = 1 - p2collapse),
             Car = selectColonies(age1$Car, p = 1 - p2collapse))

if (year > 1) {
  age2 <- list(Mel = selectColonies(age2$Mel, p = 1 - p2collapse),
               Car = selectColonies(age2$Car, p = 1 - p2collapse))
  
}

# Merge all age 0 colonies (from both periods)
age0 <- list(Mel = c(age0p1$Mel, age0p2$Mel),
             Car = c(age0p1$Car, age0p2$Car))

# Period3 ------------------------------------------------------------------
# Collapse age0 queens
print("PERIOD 3")
print("Collapse colonies, P3")
print(Sys.time())

age0 <- list(Mel = selectColonies(age0$Mel, p = (1 - p3collapseAge0)),
             Car = selectColonies(age0$Car, p = (1 - p3collapseAge0)))

age1 <- list(Mel = selectColonies(age1$Mel, p = (1 - p3collapseAge1)),
             Car = selectColonies(age1$Car, p = (1 - p3collapseAge1)))

age2 <- list(Mel = NULL, Car = NULL)#We don't need this but just to show the workflow!!!



# Maintain the number of colonies ------------------------------------------
# Keep all of age1, age0 swarmed so we build it up with some splits, while we remove (sell) the other splits
print("Maintain the number, P2")
print(Sys.time())

age0$Mel <- maintainIrelandSize(age0 = age0$Mel, age1 = age1$Mel)
age0$Car <- maintainCarSize(age0 = age0$Car, age1 = age1$Car)


for (subspecies in c("Mel", "Car")) {     #,"Lig"
  if ((nColonies(age0[[subspecies]]) + nColonies(age1[[subspecies]])) == IrelandSize
      | (nColonies(age0[[subspecies]]) + nColonies(age1[[subspecies]])) == CarSize)
  {
  } 
  else 
  {stop(paste0("The number of colonies for ", subspecies, " does not match the population size!"))}
}


} # Year-loop

locationsDF <- data.frame(Location = getLocation(c(age1$Mel, age0$Mel), collapse = TRUE),
                          Colony = c(rep("Age1", nColonies(age1$Mel)),
                                     c(rep("Age0", nColonies(age0$Mel)))))
ggplot(data = locationsDF, aes(x = Location.1, y = Location.2, colour = Colony)) + 
  geom_point()



a <- toc()
loopTime <- rbind(loopTime, c(Rep, a$tic, a$toc, a$msg, (a$toc - a$tic)))

} # Rep-loop


#I would try with raster plots for what I want to show although the heat map is acting like a raster, try both options.


# Load necessary library
library(ggplot2)

# Sample data
data <- expand.grid(x = 1:10, y = 1:10)
data$z <- runif(100, 1, 100)

# Create the heat map
ggplot(data, aes(x = x, y = y, fill = z)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "red") +
  theme_minimal()
#z seria mi introgresion