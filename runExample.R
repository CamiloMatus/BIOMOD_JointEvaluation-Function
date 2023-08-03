#Example based on two virtual species.
#Models are fitted for each species using five different algorithms using biomod2 package.
#In addition, two types of ensemble models are included.

    
    library(biomod2)
    require(PresenceAbsence)
    require(stringr)
    library(dplyr)
    library(purrr)
    source("BIOMOD_JointEvaluation.R")
    load("ws.RData")
    

    #Virtual species (two species)
    sp_occ 
    
    #Virtual environment variables
    myExpl


    
    #Loop forward by species to modeling using biomod2 package

    indModels <- list() #in this list are saved the biomod2 modeling results by species
    for (i in c(1:2)) {

    # the name of species to modeling
    myRespName.h <- myRespName[i]

    #The presence/absences data for our species
    myResp.h <- as.numeric(sp_occ[,i])


    #Formatting Data
    myBiomodData <- BIOMOD_FormatingData(resp.var = myResp.h,
                                         expl.var = myExpl,
                                         resp.xy = myRespXY,
                                         resp.name = myRespName.h)

    #Modeling through five different algorithm
    myBiomodModelOut <- BIOMOD_Modeling(myBiomodData,
                                        models=  c("GLM", "FDA", "RF", "MARS","SRE"), ##These are the algorithms.
                                        NbRunEval=2, ## It should be at least 10. Here only 2 to make the process faster
                                        DataSplit=80,
                                        modeling.id='modelNo9999')

    #Saving models in a R list object
    indModels[[i]] <- assign(paste("modelsAjus",i,myRespName.h,sep = "."),myBiomodModelOut)}



    #Doing Ensemble Modelling

    #Loop forward by species
    ensModels <- list() #in this list are saved the biomod2 ensembles results by species
    for (i in c(1:2)) {

    # the name of species to modeling
    myRespName.h <- myRespName[i]


    myBiomodModelOut <-   get(paste("modelsAjus",i,myRespName.h,sep = "."))

    #Making ensemble through two types of ensemble technique
    myBiomodEM <- BIOMOD_EnsembleModeling(modeling.output = myBiomodModelOut,
                                          prob.mean = T,
                                          committee.averaging = T)

    #Saving ensemble models in R list object
    ensModels[[i]] <- assign(paste("modelsEnemble",i,myRespName.h,sep = "."),myBiomodEM)}

    
    
    
    #Here we make the joint evaluation
    
    BIOMOD_JointEvaluation(models=indModels, ensMods=ensModels)
    #These are the most suitable algorithms according to the joint evaluation
    #Check out the file outputs in the working folder!
    
    