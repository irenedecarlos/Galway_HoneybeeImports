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
#nLig = 150
nChr = 1                     # Number of chromomsome
nDronesPerQueen = 50
nSegSites = 100              # Number of segregating sites

# Population parameters -------------------------------------------------------------------
nRep <- 1                     # Number of repeats
nYear <- 10                   # Number of years
IrelandSize<-300              #Ireland population size 300
CarSize<-100                  #Carnica pop size 100
#LigSize<-100                  #Ligustica pop size
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
getIbdHaplo(drones$Mel)[7500,]#los drones van hasta 900
#,Lig = createDrones(x = virginQueens$Car[(LigSize+1):nLig], nInd = nDronesPerQueen)
# Get fathers for Mel, MelCross and Car
fathersMel <- pullDroneGroupsFromDCA(drones$Mel, n = nInd(virginQueens$Mel[1:IrelandSize]), nDrones = nFathersPoisson) #crea 300 grupos de drones con 10 drones por grupo
fathersCar <- pullDroneGroupsFromDCA(drones$Car, n = nInd(virginQueens$Car[1:CarSize]), nDrones = nFathersPoisson)
#fathersLig <- pullDroneGroupsFromDCA(drones$Lig, n = nInd(virginQueens$Lig[1:LigSize]), nDrones = nFathersPoisson)
# Mate virgin queens with fathers to make them queens
queens <- list(Mel = SIMplyBee::cross(x = virginQueens$Mel[1:IrelandSize], drones = fathersMel),
               Car = SIMplyBee::cross(x = virginQueens$Car[1:CarSize], drones = fathersCar))
#,Lig = SIMplyBee::cross(x = virginQueens$Lig[1:LigSize], drones = fathersLig)


#Set allele frequency for queens
tmp <- c(virginQueens$Mel, virginQueens$Car) #, virginQueens$Lig

alleleFreqBaseQueens <- calcBeeAlleleFreq(x = getSegSiteGeno(tmp),
                                          sex = tmp@sex)

alleleFreqBaseQueensCar <- calcBeeAlleleFreq(x = getSegSiteGeno(virginQueens$Car),
                                             sex = virginQueens$Car@sex)

#alleleFreqBaseQueensLig <- calcBeeAlleleFreq(x = getSegSiteGeno(virginQueens$Lig),
#sex = virginQueens$Lig@sex)

alleleFreqBaseQueensMel <- calcBeeAlleleFreq(x = getSegSiteGeno(virginQueens$Mel),
                                             sex = virginQueens$Mel@sex)

#Get allele freq for csd locus
csdLocus <- paste0(SP$csdChr, "_", SP$csdPosStart:SP$csdPosStop)
alleleFreqCsdLocusBaseQueens <- alleleFreqBaseQueens[csdLocus]
alleleFreqCsdLocusBaseCar <- alleleFreqBaseQueensCar[csdLocus]
#alleleFreqCsdLocusBaseLig <- alleleFreqBaseQueensLig[csdLocus]
alleleFreqCsdLocusBaseMel <- alleleFreqBaseQueensMel[csdLocus]

#Get allele freq for csd Chromosome - this pulls out only the 3rd chromosome
alleleFreqCsdChrBaseQueens <- t(as.data.frame(alleleFreqBaseQueens))[, grepl(pattern = paste0("^", csdChr, "_"), x = colnames(t(as.data.frame(alleleFreqBaseQueens))))] %>% t()
alleleFreqCsdChrBaseCar <- t(as.data.frame(alleleFreqBaseQueensCar))[, grepl(pattern = paste0("^", csdChr, "_"), x = colnames(t(as.data.frame(alleleFreqBaseQueensCar))))] %>% t()
#alleleFreqCsdChrBaseLig <- t(as.data.frame(alleleFreqBaseQueensLig))[, grepl(pattern = paste0("^", csdChr, "_"), x = colnames(t(as.data.frame(alleleFreqBaseQueensLig))))] %>% t()
alleleFreqCsdChrBaseMel <- t(as.data.frame(alleleFreqBaseQueensMel))[, grepl(pattern = paste0("^", csdChr, "_"), x = colnames(t(as.data.frame(alleleFreqBaseQueensMel))))] %>% t()

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
  first<-getId(age1$Mel[1:150])
  second<-getId(age1$Mel[151:300])
  age1$Mel[1:150] <-  setLocation(age1$Mel[1:150], 
                             location = Map(c, runif(nColonies(age1$Mel[1:150]), 0, pi), runif(nColonies(age1$Mel[1:150]), 0, 2*pi)))
  age1$Mel[151:300] <-  setLocation(age1$Mel[151:300], 
                                  location = Map(c, runif(nColonies(age1$Mel[151:300]), pi, 2*pi), runif(nColonies(age1$Mel[151:300]), 0, 2*pi)))

  
  locationsDF <- data.frame(Location = getLocation(c(age1$Mel[1:150], age1$Mel[151:300]), collapse = TRUE),
                            Beekeeper = c(rep("Beekeeper1", nColonies(age1$Mel[1:150])),
                                          rep("Beekeeper3", nColonies(age1$Mel[151:300]))))
  
  ggplot(data = locationsDF, aes(x = Location.1, y = Location.2, colour = Beekeeper)) + 
    geom_point()
  
  
  print("Record initial colonies")
  colonyRecords <- data_rec(datafile = colonyRecords, colonies = age1$Mel, year = year, population = "Mel")
  colonyRecords <- data_rec(datafile = colonyRecords, colonies = age1$Car, year = year, population = "Car")
  #colonyRecords <- data_rec(datafile = colonyRecords, colonies = age1$Lig, year = year, population = "Lig")
  
  # If not, promote the age0 to age1, age1 to age2 and remove age2 colonies
} else {
  age2 <- list(Mel = age1$Mel, Car = age1$Car) #, Lig = age1$Lig
  age1 <- list(Mel = age0$Mel, Car = age0$Car) #, Lig = age0$Lig
  age0 <- list(Mel = NULL, Car = NULL) #, Lig = NULL
  age0p1 <- list(Mel = NULL, Car = NULL) #, Lig = NULL
  age0p2 <- list(Mel = NULL, Car = NULL) #, Lig = NULL
}

# Period1 ------------------------------------------------------------------
# Build-up the colonies
print(paste0("Building up the colonies to ", nWorkers, " and ", nDrones))
print(Sys.time())
age1 <- list(Mel = buildUp(age1$Mel),
             Car = buildUp(age1$Car))
#,Lig = buildUp(age1$Lig))
if (year > 1) {
  age2 <- list(Mel = buildUp(age2$Mel),
               Car = buildUp(age2$Car))
  #,Lig = buildUp(age2$Lig))
}

# Split all age1 colonies
print("Splitting the colonies")
print(Sys.time())
tmp <- list(Mel = split(age1$Mel),
            Car = split(age1$Car))
#,Lig = split(age1$Lig))
age1 <- list(Mel = tmp$Mel$remnant,
             Car = tmp$Car$remnant)
#,Lig = tmp$Lig$remnant)

# The queens of the splits are 0 years old
age0p1 <- list(Mel = tmp$Mel$split, Car = tmp$Car$split) #,Lig = tmp$Lig$split

idfirst<-getId(age0p1$Mel[1:150])
idsecond<-getId(age0p1$Mel[151:300])
age0p1$Mel[1:150] <-  setLocation(age0p1$Mel[1:150], 
                                location = Map(c, runif(nColonies(age0p1$Mel[1:150]), 0, pi), runif(nColonies(age0p1$Mel[1:150]), 0, 2*pi)))
age0p1$Mel[151:300] <-  setLocation(age0p1$Mel[151:300], 
                                  location = Map(c, runif(nColonies(age0p1$Mel[151:300]), pi, 2*pi), runif(nColonies(age0p1$Mel[151:300]), 0, 2*pi)))


locationsDFp1 <- data.frame(Location = getLocation(c(age0p1$Mel[1:150], age0p1$Mel[151:300]), collapse = TRUE),
                          Beekeeper = c(rep("Beekeeper1", nColonies(age0p1$Mel[1:150])),
                                        rep("Beekeeper3", nColonies(age0p1$Mel[151:300]))))

ggplot(data = locationsDFp1, aes(x = Location.1, y = Location.2, colour = Beekeeper)) + 
  geom_point()




 if (year > 1) {
  # Split all age2 colonies
  tmp <- list(Mel = split(age2$Mel),
              Car = split(age2$Car))
  #,Lig = split(age2$Lig))
  age2 <- list(Mel = tmp$Mel$remnant,
               Car = tmp$Car$remnant
  ) #,Lig = tmp$Lig$remnant
  # The queens of the splits are 0 years old
  age0p1 <- list(Mel = c(age0p1$Mel, tmp$Mel$split),
                 Car = c(age0p1$Car, tmp$Car$split))
  #,Lig = c(age0p1$Lig, tmp$Lig$split))
}


print("Create virgin queens, period 1")
print(Sys.time())
virginDonor <- list(Mel = sample.int(n = nColonies(age1$Mel), size = 1),
                    Car = sample.int(n = nColonies(age1$Car), size = 1))
#Lig = sample.int(n = nColonies(age1$Lig), size = 1))
# Virgin queens for splits!
pImport<-0.3

tmp <- (Mel = pullColonies(age0p1$Mel, p=pImport))
IdImportColonies<-getId(tmp$pulled)
age0p1 <- list(Mel = tmp$remnant,
               MelImport = tmp$pulled,
               Car = c(age0p1$Car, tmp$Car$split))
#,Lig = c(age0p1$Lig, tmp$Lig$split))
virginQueens <- list(Mel = createVirginQueens(age1$Mel[[virginDonor$Mel]], nInd = nColonies(age0p1$Mel)),
                     Car = createVirginQueens(age1$Car[[virginDonor$Car]], nInd = nColonies(age0p1$Car)+nColonies(age0p1$MelImport)))
#,Lig = createVirginQueens(age1$Lig[[virginDonor$Lig]], nInd = nColonies(age0p1$Lig)+(nColonies(age0p1$MelImport)/2)))

# Requeen the splits --> queens are now 0 years old


nColoniesMelImport<-nColonies(age0p1$MelImport)
nColoniesCar<-nColonies(age0p1$Car)+nColonies(age0p1$MelImport)

age0p1 <- list(Mel = c(reQueen(age0p1$Mel, queen = (virginQueens$Mel)) ,
                       reQueen(age0p1$MelImport, queen = c((virginQueens$Car)[1:nColoniesMelImport]))),       #,(virginQueens$Lig)[1:(nColoniesMelImport/2)]))),
               Car = reQueen(age0p1$Car, queen = virginQueens$Car[(nColoniesMelImport+1):nColoniesCar]))

caca<-pullColonies(age0p1$Mel, ID=idfirst)
firsthalf<-caca$pulled
secondhalf<-caca$remnant
locationsDFp1 <- data.frame(Location = getLocation(c(firsthalf, secondhalf), collapse = TRUE),
                            Beekeeper = c(rep("Beekeeper1", nColonies(firsthalf)),
                                          rep("Beekeeper3", nColonies(secondhalf))))

ggplot(data = locationsDFp1, aes(x = Location.1, y = Location.2, colour = Beekeeper)) + 
  geom_point()


#ggplot(data = locationsDFp1, aes(x = Location.1, y = Location.2, colour = Beekeeper)) + 
#  geom_point()
# Swarm a percentage of age1 colonies
print("Swarm colonies, P1")
print(Sys.time())
tmp <- list(Mel = pullColonies(age1$Mel, p = p1swarm),
            Car = pullColonies(age1$Car, p = p1swarm))
lower150<-which((getId(tmp$Mel$pulled)<=150))  #with this I am selecting the ids of the colonies that swarm that are lower than 150 and greater
greater150<-which((getId(tmp$Mel$pulled)>150)) 
age1 <- list(Mel = tmp$Mel$remnant,
             Car = tmp$Car$remnant)
#,Lig = tmp$Lig$remnant)
tmp <- list(Mel = swarm(tmp$Mel$pulled, sampleLocation = T, radius = 0.3),
            Car = swarm(tmp$Car$pulled, sampleLocation = T, radius = 0.3))
swarms2<-getId(tmp$Mel$swarm)
swarms3<-getId(tmp$Mel$remnant)
#,Lig = swarm(tmp$Lig$pulled))
age0p1 <- list(Mel = c(age0p1$Mel, tmp$Mel$remnant),
               Car = c(age0p1$Car, tmp$Car$remnant))
#,Lig = c(age0p1$Lig, tmp$Lig$remnant))
age1 <- list(Mel = c(age1$Mel, tmp$Mel$swarm),
             Car = c(age1$Car, tmp$Car$swarm))

age1idfirst<-c(first,swarms2[lower150])      # and here im adding those colonies to the previous list of IDs for age1
age1idsecond<-c(second,swarms2[greater150])
#to plot their location
#first2<-pullColonies(age1$Mel, ID = c(first,swarms2[lower150]))
#first3<-first2$pulled
#second2<-first2$remnant
age0p1idfirst<-c(idfirst,swarms3[lower150]) #and here for age0p1
age0p1idsecond<-c(idsecond,swarms3[greater150])
#idfirst2<-pullColonies(age0p1$Mel, ID = c(idfirst, swarms3[lower150]))
#idfirst3<-idfirst2$pulled
#idsecond2<-idfirst2$remnant

#firsthalf y secondhalf son age0p1 y les has añadido algo de age1. hay que añadirles las remnant
#locationsDF1 <- data.frame(Location = getLocation(c(idfirst3, idsecond2), collapse = TRUE),
#                            Beekeeper = c(rep("Beekeeper1", nColonies(idfirst3)),
#                                          rep("Beekeeper3", nColonies(idsecond2))))

#ggplot(data = locationsDF1, aes(x = Location.1, y = Location.2, colour = Beekeeper)) + 
#  geom_point()

if (year > 1) {
  # Swarm a percentage of age2 colonies
  tmp <- list(Mel = pullColonies(age2$Mel, p = p1swarm),
              Car = pullColonies(age2$Car, p = p1swarm))
  #, Lig = pullColonies(age2$Lig, p = p1swarm))
  age2 <- list(Mel = tmp$Mel$remnant,     
               Car = tmp$Car$remnant)
  #,Lig = tmp$Lig$remnant)
  tmp <- list(Mel = swarm(tmp$Mel$pulled),
              Car = swarm(tmp$Car$pulled))
  #,Lig = swarm(tmp$Lig$pulled))
  age0p1 <- list(Mel = c(age0p1$Mel, tmp$Mel$remnant),
                 Car = c(age0p1$Car, tmp$Car$remnant))
  #,Lig = c(age0p1$Lig, tmp$Lig$remnant))
  age2 <- list(Mel = c(age2$Mel, tmp$Mel$swarm),
               Car = c(age2$Car, tmp$Car$swarm))
  #,Lig = c(age2$Lig, tmp$Lig$swarm))
}

# Supersede age1 colonies
print("Supersede colonies, P1")
print(Sys.time())
tmp <- list(Mel = pullColonies(age1$Mel, p = p1supersede),
            Car = pullColonies(age1$Car, p = p1supersede))
lower150s<-(getId(tmp$Mel$pulled)[getId(tmp$Mel$pulled)<=150])
greater150s<-(getId(tmp$Mel$pulled)[151:300])

age1 <- list(Mel = tmp$Mel$remnant,
             Car = tmp$Car$remnant)

tmp <- list(Mel = supersede(tmp$Mel$pulled),
            Car = supersede(tmp$Car$pulled))

age0p1 <- list(Mel = c(age0p1$Mel, tmp$Mel),
               Car = c(age0p1$Car, tmp$Car))

age0p1idfirst<-c(age0p1idfirst,lower150s)
age0p1idsecond<-c(age0p1idsecond,greater150s)

patufo<-pullColonies(age0p1$Mel, ID = age0p1idfirst)
patufo2<-patufo$pulled
patufo3<-patufo$remnant

locationsDFp1 <- data.frame(Location = getLocation(c(patufo2, patufo3), collapse = TRUE),
                            Beekeeper = c(rep("Beekeeper1", nColonies(patufo2)),
                                          rep("Beekeeper3", nColonies(patufo3))))

ggplot(data = locationsDFp1, aes(x = Location.1, y = Location.2, colour = Beekeeper)) + 
  geom_point()

if (year > 1) {
  # Supersede age2 colonies
  tmp <- list(Mel = pullColonies(age2$Mel, p = p1supersede),
              Car = pullColonies(age2$Car, p = p1supersede))
  #,Lig = pullColonies(age2$Lig, p = p1supersede))
  age2 <- list(Mel = tmp$Mel$remnant,
               Car = tmp$Car$remnant)
  #,Lig = tmp$Lig$remnant)
  tmp <- list(Mel = supersede(tmp$Mel$pulled),
              Car = supersede(tmp$Car$pulled))
  #,Lig = supersede(tmp$Lig$pulled))
  age0p1 <- list(Mel = c(age0p1$Mel, tmp$Mel),
                 Car = c(age0p1$Car, tmp$Car))
  #,Lig = c(age0p1$Lig, tmp$Lig))
}


# Mate the split colonies
print("Mate split colonies, P1")
print(Sys.time())
if (year == 1) {
  DCAMel <- createDCA(age1$Mel)
  age0p1$Mel <- cross(age0p1$Mel, drones = pullDroneGroupsFromDCA(DCA = DCAMel, n = nColonies(age0p1$Mel), nDrones = nFathersPoisson))
  DCACar <- createDCA(age1$Car)
  age0p1$Car <- cross(age0p1$Car, drones = pullDroneGroupsFromDCA(DCA = DCACar, n = nColonies(age0p1$Car), nDrones = nFathersPoisson))
  #DCALig <- createDCA(age1$Lig)
  #age0p1$Lig <- cross(age0p1$Lig, drones = pullDroneGroupsFromDCA(DCA = DCALig, n = nColonies(age0p1$Lig), nDrones = nFathersPoisson))
} else {
  DCAMel <- createDCA(c(age1$Mel, age2$Mel))
  age0p1$Mel <- cross(age0p1$Mel, drones = pullDroneGroupsFromDCA(DCA = DCAMel, n = nColonies(age0p1$Mel), nDrones = nFathersPoisson))
  DCACar <- createDCA(c(age1$Car, age2$Car))
  age0p1$Car <- cross(age0p1$Car, drones = pullDroneGroupsFromDCA(DCA = DCACar, n = nColonies(age0p1$Car), nDrones = nFathersPoisson))
  #DCALig <- createDCA(c(age1$Lig, age2$Lig)) 
  #age0p1$Lig <- cross(age0p1$Lig, drones = pullDroneGroupsFromDCA(DCA = DCALig, n = nColonies(age0p1$Lig), nDrones = nFathersPoisson))
}


# Collapse
print("Collapse colonies, P1")
print(Sys.time())
age1 <- list(Mel = selectColonies(age1$Mel, p = 1 - p1collapse),
             Car = selectColonies(age1$Car, p = 1 - p1collapse))
#,Lig = selectColonies(age1$Lig, p = 1 - p1collapse))
if (year > 1) {
  age2 <- list(Mel = selectColonies(age2$Mel, p = 1 - p1collapse),
               Car = selectColonies(age2$Car, p = 1 - p1collapse))
  #,Lig = selectColonies(age2$Lig, p = 1 - p1collapse))
}


# Period2 ------------------------------------------------------------------
print("PERIOD 2")
# Swarm a percentage of age1 colonies
# Mellifera
print("Swarm colonies, P2")
print(Sys.time())
tmp <- list(Mel = pullColonies(age1$Mel, p = p2swarm),
            Car = pullColonies(age1$Car, p = p2swarm))
#,Lig = pullColonies(age1$Lig, p = p2swarm))
age1 <- list(Mel = tmp$Mel$remnant,
             Car = tmp$Car$remnant)
#,Lig = tmp$Lig$remnant)
tmp <- list(Mel = swarm(tmp$Mel$pulled),
            Car = swarm(tmp$Car$pulled))
#,Lig = swarm(tmp$Lig$pulled))
# The queens of the remnant colonies are of age 0
age0p2 <- list(Mel = tmp$Mel$remnant,
               Car = tmp$Car$remnant)
#,Lig = tmp$Lig$remnant)
age1 <- list(Mel = c(age1$Mel, tmp$Mel$swarm),
             Car = c(age1$Car, tmp$Car$swarm))
#,Lig = c(age1$Lig, tmp$Lig$swarm))

if (year > 1) {
  # Swarm a percentage of age2 colonies
  tmp <- list(Mel = pullColonies(age2$Mel, p = p2swarm),
              Car = pullColonies(age2$Car, p = p2swarm))
  #,Lig = pullColonies(age2$Lig, p = p2swarm))
  age2 <- list(Mel = tmp$Mel$remnant,
               Car = tmp$Car$remnant)
  #,Lig = tmp$Lig$remnant)
  tmp <- list(Mel = swarm(tmp$Mel$pulled),
              Car = swarm(tmp$Car$pulled))
  #,Lig = swarm(tmp$Lig$pulled))
  # The queens of the remnant colonies are of age 0
  age0p2 <- list(Mel = tmp$Mel$remnant,
                 Car = tmp$Car$remnant)
  #,Lig = tmp$Lig$remnant)
  age2 <- list(Mel = c(age2$Mel, tmp$Mel$swarm),
               Car = c(age2$Car, tmp$Car$swarm))
  #,Lig = c(age2$Lig, tmp$Lig$swarm))
}

# Supersede a part of age1 colonies
print("Supersede colonies, P2")
print(Sys.time())

tmp <- list(Mel = pullColonies(age1$Mel, p = p2supersede),
            Car = pullColonies(age1$Car, p = p2supersede))
#,Lig = pullColonies(age1$Lig, p = p2supersede))
age1 <- list(Mel = tmp$Mel$remnant,
             Car = tmp$Car$remnant)
#,Lig = tmp$Lig$remnant)
tmp <- list(Mel = supersede(tmp$Mel$pulled),
            Car = supersede(tmp$Car$pulled))
#,Lig = supersede(tmp$Lig$pulled))
# The queens of superseded colonies are of age 0
age0p2 <- list(Mel = c(age0p2$Mel, tmp$Mel),
               Car = c(age0p2$Car, tmp$Car))
#,Lig = c(age0p2$Lig, tmp$Lig))

if (year > 1) {
  # Supersede a part of age2 colonies
  tmp <- list(Mel = pullColonies(age2$Mel, p = p2supersede),
              Car = pullColonies(age2$Car, p = p2supersede))
  #,Lig = pullColonies(age2$Lig, p = p2supersede))
  age2 <- list(Mel = tmp$Mel$remnant,
               Car = tmp$Car$remnant)
  #Lig = tmp$Lig$remnant)
  tmp <- list(Mel = supersede(tmp$Mel$pulled),
              Car = supersede(tmp$Car$pulled))
  #Lig = supersede(tmp$Lig$pulled))
  # The queens of superseded colonies are of age 0
  age0p2 <- list(Mel = c(age0p2$Mel, tmp$Mel),
                 Car = c(age0p2$Car, tmp$Car))
  #, Lig = c(age0p2$Lig, tmp$Lig))
}

# Replace all the drones
print("Replace Drones, P2")
print(Sys.time())

age1$Mel <- replaceDrones(age1$Mel)
age1$Car <- replaceDrones(age1$Car)
#age1$Lig <- replaceDrones(age1$Lig)
if (year > 1) {
  age2$Mel <- replaceDrones(age2$Mel)
  age2$Car <- replaceDrones(age2$Car)
  #age2$Lig <- replaceDrones(age2$Lig)
}

# Mate the colonies
# Import p percentage of carnica colonies into mellifera DCA
print("Mate colonies, P2")
print(Sys.time())


if (year == 1) {
  DCAMel <- createDCA(age1$Mel)
  age0p2$Mel <- cross(age0p2$Mel, drones = pullDroneGroupsFromDCA(DCA = DCAMel, n = nColonies(age0p2$Mel), nDrones = nFathersPoisson))
  DCACar <- createDCA(age1$Car)
  age0p2$Car <- cross(age0p2$Car, drones = pullDroneGroupsFromDCA(DCA = DCACar, n = nColonies(age0p2$Car), nDrones = nFathersPoisson))
  #DCALig <- createDCA(age1$Lig)
  #age0p2$Lig <- cross(age0p2$Lig, drones = pullDroneGroupsFromDCA(DCA = DCALig, n = nColonies(age0p2$Lig), nDrones = nFathersPoisson))
} else {
  DCAMel <- createDCA(c(age1$Mel, age2$Mel))
  fathersMel <- pullDroneGroupsFromDCA(DCA = DCAMel, n = nColonies(age0p2$Mel), nDrones = nFathersPoisson)
  fathersMel[[1]] <- c(fathersMel[[1]], createDrones(age1$Mel[[1]], nInd = 2))
  age0p2$Mel <- cross(age0p2$Mel, drones = fathersMel)
  DCACar <- createDCA(c(age1$Car, age2$Car))
  fathersCar <-  pullDroneGroupsFromDCA(DCA = DCACar, n = nColonies(age0p2$Car), nDrones = nFathersPoisson)
  fathersCar[[1]] <- c(fathersCar[[1]], createDrones(age1$Car[[1]], nInd = 2))
  age0p2$Car <- cross(age0p2$Car, drones = fathersCar)
  #DCALig <- createDCA(c(age1$Lig, age2$Lig))
  #fathersLig <-  pullDroneGroupsFromDCA(DCA = DCALig, n = nColonies(age0p2$Lig), nDrones = nFathersPoisson)
  #fathersLig[[1]] <- c(fathersLig[[1]], createDrones(age1$Lig[[1]], nInd = 2))
  #age0p2$Lig <- cross(age0p2$Lig, drones = fathersLig)
}

# Collapse
age1 <- list(Mel = selectColonies(age1$Mel, p = 1 - p2collapse),
             Car = selectColonies(age1$Car, p = 1 - p2collapse))
#,Lig = selectColonies(age1$Lig, p = 1 - p2collapse))
if (year > 1) {
  age2 <- list(Mel = selectColonies(age2$Mel, p = 1 - p2collapse),
               Car = selectColonies(age2$Car, p = 1 - p2collapse))
  #,Lig = selectColonies(age2$Lig, p = 1 - p2collapse))
}

# Merge all age 0 colonies (from both periods)
age0 <- list(Mel = c(age0p1$Mel, age0p2$Mel),
             Car = c(age0p1$Car, age0p2$Car))

#,Lig = c(age0p1$Lig, age0p2$Lig))
colonyRecords <- data_rec(datafile = colonyRecords, colonies = age0$Mel, year = year, population = "Mel")
colonyRecords <- data_rec(datafile = colonyRecords, colonies = age0$Car, year = year, population = "Car")
#colonyRecords <- data_rec(datafile = colonyRecords, colonies = age0$Lig, year = year, population = "Lig")
# Period3 ------------------------------------------------------------------
# Collapse age0 queens
print("PERIOD 3")
print("Collapse colonies, P3")
print(Sys.time())

age0 <- list(Mel = selectColonies(age0$Mel, p = (1 - p3collapseAge0)),
             Car = selectColonies(age0$Car, p = (1 - p3collapseAge0)))
#,Lig = selectColonies(age0$Lig, p = (1 - p3collapseAge0)))
age1 <- list(Mel = selectColonies(age1$Mel, p = (1 - p3collapseAge1)),
             Car = selectColonies(age1$Car, p = (1 - p3collapseAge1)))
#,Lig = selectColonies(age1$Lig, p = (1 - p3collapseAge1)))
age2 <- list(Mel = NULL, Car = NULL) #,Lig=NULL)#We don't need this but just to show the workflow!!!




# Maintain the number of colonies ------------------------------------------
# Keep all of age1, age0 swarmed so we build it up with some splits, while we remove (sell) the other splits
print("Maintain the number, P2")
print(Sys.time())

age0$Mel <- maintainIrelandSize(age0 = age0$Mel, age1 = age1$Mel)
age0$Car <- maintainCarSize(age0 = age0$Car, age1 = age1$Car)
#age0$Lig <- maintainCarSize(age0 = age0$Lig, age1 = age1$Lig)

for (subspecies in c("Mel", "Car")) {     #,"Lig"
  if ((nColonies(age0[[subspecies]]) + nColonies(age1[[subspecies]])) == IrelandSize
      | (nColonies(age0[[subspecies]]) + nColonies(age1[[subspecies]])) == CarSize)
  {
  } 
  else 
  {stop(paste0("The number of colonies for ", subspecies, " does not match the population size!"))}
}
Carsurvqueens = mergePops(getQueen(age0$Mel))
IBD0<-apply(getIbdHaplo(Carsurvqueens),MARGIN = 1, FUN =  function(X) sum(X %in% 1:(nMelN*2)/length(X)))
IBDalive<-sapply(seq(1,length(IBD0),2), FUN = function(z) sum(IBDh[z:(z+1)])/2)
#function that looks at IBDalive and returns how many values are 0
Carsurviving<-sum(IBDalive==0)


} # Year-loop

a <- toc()
loopTime <- rbind(loopTime, c(Rep, a$tic, a$toc, a$msg, (a$toc - a$tic)))

} # Rep-loop

print("Saving image data")
save.image("SpringerSimulation_import.RData")

queens = mergePops(getQueen(age0$Mel))
getIbdHaplo(queens)
IBDh<-apply(getIbdHaplo(queens),MARGIN = 1, FUN =  function(X) sum(X %in% 1:(nMelN*2)/length(X)))
IBD<-sapply(seq(1,length(IBDh),2), FUN = function(z) sum(IBDh[z:(z+1)])/2)
IBD
getIbdHaplo(queens)
mean(IBD)
var(IBD)

queen1<- mergePops(getQueen(age1$Mel))
getIbdHaplo(queen1)[294,]
IBDh<-apply(getIbdHaplo(queen1),MARGIN = 1, FUN =  function(X) sum(X %in% 1:(nMelN*2)/length(X)))
IBD<-sapply(seq(1,length(IBDh),2), FUN = function(z) sum(IBDh[z:(z+1)])/2)
IBD
sum(IBD==1)
workers1<- mergePops(getWorkers(age0$Mel))
workers2<- mergePops(getWorkers(age0$Mel[1,]))
getIbdHaplo(workers1)[980,]


# hacer esto con los wrokers para ver cuantos haplocodes salen, quiero ver si sale el uno y dos o que ostias esta pasando
worker<- mergePops(getWorkers(age1$Mel))
workersHaploCodes<-unique(getIbdHaplo(worker))
hapl<-as.data.frame.array(workersHaploCodes,row.names=NULL)
HaploCodes<-hapl[,1]
Haploorden<-sort(HaploCodes)
IBDh<-apply(getIbdHaplo(worker),MARGIN = 1, FUN =  function(X) sum(X %in% 1:(nMelN*2)/length(X)))
IBD<-sapply(seq(1,length(IBDh),2), FUN = function(z) sum(IBDh[z:(z+1)])/2)
IBD
#contar cuantos elementos son igual a 0.5
sum(IBD==0.5)
length(IBD)
sum(IBD==1)
#A ver lo mismo con drones:
dron<- mergePops(getDrones(age0$Car))
getIbdHaplo(dron)[1155,]
dronHaploCodes<-unique(getIbdHaplo(dron))
hapl<-as.data.frame.array(dronHaploCodes,row.names=NULL)
HaploCodes<-hapl[,1]
Haploorden<-sort(HaploCodes)
IBDh<-apply(getIbdHaplo(dron),MARGIN = 1, FUN =  function(X) sum(X %in% 1:(nMelN*2)/length(X)))
IBD<-sapply(seq(1,length(IBDh),2), FUN = function(z) sum(IBDh[z:(z+1)])/2)
IBD
# I think it is okey if we put 1:(nMelN)*2 because all queens of mellifera will 
#have an haplotype from 1 to twice the founder genomes o Mel (nMelN*2)
colonyRecords[800:900,]
