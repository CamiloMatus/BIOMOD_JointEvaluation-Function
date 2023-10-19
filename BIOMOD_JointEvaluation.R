

BIOMOD_JointEvaluation <- function(models=x,ensMods=NULL,metrics=NULL,weights=NULL){ 

  if(file.exists("algoRankResults")==F){dir.create("algoRankResults")}else{"OK"}

  if(!is.null(metrics)&is.null(weights)){

    weights <- rep(1/length(metrics),length(metrics))
  }


  if(is.list(models)==T){isListMod <- "yes"}else{isListMod <- "no"}
  if(is.list(ensMods)==T){isListEnMod <- "yes"}else{isListEnMod <- "no"}



  if(isListMod!=isListEnMod&!is.null(ensMods)){message("Error: object of argument models or ensMods of class not supported")}else{




    norm10 <- function(x){(x-min(x))/(max(x)-min(x))}

    if(is.null(metrics)==TRUE&is.null(weights)==TRUE){metrics <- c("AUC","TSS","KAPPA","OF.AUC","OF.TSS","OF.KAPPA")
    
    weights <-  rep(0.1666667,6)}

    if(is.null(metrics)==TRUE&is.null(weights)==FALSE){metrics <- c("AUC","TSS","KAPPA","OF.AUC","OF.TSS","OF.KAPPA")
    }
    
    
    if(any(metrics%in% c("AUC","TSS","KAPPA","OF.AUC","OF.TSS","OF.KAPPA")==FALSE)==TRUE){

      er <- metrics[!metrics%in% c("AUC","TSS","KAPPA","OF.AUC","OF.TSS","OF.KAPPA")]
      message("Error")
      message(paste(" Wrong metric name: ",er,";"))}else{

        if(length(metrics)!=length(weights)){message(paste("Error, length of vector metrics", " (",length(metrics),") ", "does not fit length of vector weights",sep = ""))}else{
          if(class(weights)!="numeric"){message("Error, vector weights of non-numeric class")}else{
            if(is.null(ensMods)==FALSE){nEnsAlgo <- length(ensMods)}else{nEnsAlgo <- 0}
            if(sum(weights)>1.000001){warning("Sum of the weights greater than one")}
            if(sum(weights)<1){warning("Sum of the weights less than one")}


            sppsToWork <- 1:length(models)

            dat.out2 <- c()
            for (i in sppsToWork) {

              if(isListMod=="yes"){biomodIndModels.h <-   models[[i]]}else{
                biomodIndModels.h <-   models}


              myRespName.h <- biomodIndModels.h@sp.name

              resp.var <- get_formal_data(biomodIndModels.h,'resp.var')

              modPredIndi <- get_predictions(biomodIndModels.h, as.data.frame=TRUE,evaluation=F)

              calibData <-data.frame(get_calib_lines(biomodIndModels.h))[get_calib_lines(biomodIndModels.h) %>% data.frame %>%names() %>% str_detect("RUN")%in%T]


              if(nEnsAlgo>0){

                if(isListEnMod=="yes")  {
                  biomodIndModels.h2 <-   ensMods[[i]]
                  modPredEns <- get_predictions(biomodIndModels.h2, as.data.frame=TRUE,evaluation=F)
                  modPred <- cbind(modPredIndi,modPredEns)
                  namesModPred <- names(modPred)}
                else{ biomodIndModels.h2 <-   ensMods

                modPredEns <- get_predictions(biomodIndModels.h2, as.data.frame=TRUE,evaluation=F)
                modPred <- cbind(modPredIndi,modPredEns)
                namesModPred <- names(modPred)}}else{
                  modPred <- modPredIndi
                  namesModPred <- names(modPred)}

              namesModPred <- namesModPred[str_detect(namesModPred,"RUN")]


              if(nEnsAlgo>0){

                if(any(biomodIndModels.h2@sp.name==myRespName.h)==FALSE){
                  options(warn = 2)
                  warning("Error: species name of individual models is not the same as that of ensemble models")

                  }else{ options(warn = 1)}}


                algos <-c()
                if(TRUE%in%str_detect(namesModPred,"SRE")==TRUE){algos <-c(algos,"SRE") }
                if(TRUE%in%str_detect(namesModPred,"MARS")==TRUE){algos <-c(algos,"MARS") }
                if(TRUE%in%str_detect(namesModPred,"GLM")==TRUE){algos <-c(algos,"GLM") }
                if(TRUE%in%str_detect(namesModPred,"RF")==TRUE){algos <-c(algos,"RF") }
                if(TRUE%in%str_detect(namesModPred,"GBM")==TRUE){algos <-c(algos,"GBM") }
                if(TRUE%in%str_detect(namesModPred,"GAM")==TRUE){algos <-c(algos,"GAM") }
                if(TRUE%in%str_detect(namesModPred,"CTA")==TRUE){algos <-c(algos,"CTA") }
                if(TRUE%in%str_detect(namesModPred,"ANN")==TRUE){algos <-c(algos,"ANN") }
                if(TRUE%in%str_detect(namesModPred,"FDA")==TRUE){algos <-c(algos,"FDA") }
                if(TRUE%in%str_detect(namesModPred,"MAXENT.Phillips")==TRUE){algos <-c(algos,"MAXENT.Phillips") }
                if(TRUE%in%str_detect(namesModPred,"EMmeanByROC")==TRUE){algos <-c(algos,"EMmeanByROC") }
                if(TRUE%in%str_detect(namesModPred,"EMmeanByKAPPA")==TRUE){algos <-c(algos,"EMmeanByKAPPA") }
                if(TRUE%in%str_detect(namesModPred,"EMmeanByTSS")==TRUE){algos <-c(algos,"EMmeanByTSS") }
                if(TRUE%in%str_detect(namesModPred,"EMcaByROC")==TRUE){algos <-c(algos,"EMcaByROC") }
                if(TRUE%in%str_detect(namesModPred,"EMcaByKAPPA")==TRUE){algos <-c(algos,"EMcaByKAPPA") }
                if(TRUE%in%str_detect(namesModPred,"EMcaByTSS")==TRUE){algos <-c(algos,"EMcaByTSS") }
                if(TRUE%in%str_detect(namesModPred,"EMwmeanByROC")==TRUE){algos <-c(algos,"EMwmeanByROC") }
                if(TRUE%in%str_detect(namesModPred,"EMwmeanByKAPPA")==TRUE){algos <-c(algos,"EMwmeanByKAPPA") }
                if(TRUE%in%str_detect(namesModPred,"EMwmeanByTSS")==TRUE){algos <-c(algos,"EMwmeanByTSS") }

                dat.out1 <- c()

                for(w in 1:length(algos)){

                  dat.h <- modPred[,namesModPred[str_detect(namesModPred,algos[w]) ]]

                  nruns <- names(as.data.frame(get_calib_lines(biomodIndModels.h)))[get_calib_lines(biomodIndModels.h) %>% as.data.frame() %>% names() %>% str_detect("RUN")] %>% length()


                  for (e in 1:nruns) {

                    pred.h <- dat.h[, names(dat.h)[str_detect(names(dat.h),paste("RUN", e,"_",algos[w],sep=""))]]*0.001


                    if(length(pred.h)==0){
                      pred.h <- dat.h[, names(dat.h)[str_detect(names(dat.h),paste("RUN", e,"_","All",sep=""))]]*0.001
                    }


                    if(TRUE%in%is.na(pred.h)==T){warning(message(paste("Run",e, "not considered because it contains NA values")))


                    }else{message(paste("Computing predictive performance and overfitting for algorithm ",
                                        algos[w]," for RUN ",e, " (Especies: ",myRespName.h,")",sep = ""))


                      dat.h2 <- data.frame(1:length(pred.h),resp.var,pred.h)
                      th.op.h <- optimal.thresholds(dat.h2,threshold = 1000,na.rm = T,opt.methods="MaxSens+Spec")[2] %>% as.numeric()
                      dat.h2$presBin <- NA
                      dat.h2[dat.h2[!dat.h2[,3]%in%NA,3]>=th.op.h,"presBin"] <- 1
                      dat.h2[dat.h2[!dat.h2[,3]%in%NA,3]<th.op.h,"presBin"] <- 0
                      calibData.h <- calibData[str_detect(names(calibData),paste("RUN", e,"._",sep=""))]
                      dat.h2$calib <- calibData.h

                      dat.hTest <- dat.h2[dat.h2[,5]==FALSE,]
                      dat.hCalib <- dat.h2[dat.h2[,5]==TRUE,]

                      TPtest <- nrow(dat.hTest[dat.hTest[,2]==1&dat.hTest[,4]==1,])
                      FPtest <- nrow(dat.hTest[dat.hTest[,2]==0&dat.hTest[,4]==1,])
                      FNtest <- nrow(dat.hTest[dat.hTest[,2]==1&dat.hTest[,4]==0,])
                      TNtest <- nrow(dat.hTest[dat.hTest[,2]==0&dat.hTest[,4]==0,])

                      TPcalib <- nrow(dat.hCalib[dat.hCalib[,2]==1&dat.hCalib[,4]==1,])
                      FPcalib <- nrow(dat.hCalib[dat.hCalib[,2]==0&dat.hCalib[,4]==1,])
                      FNcalib <- nrow(dat.hCalib[dat.hCalib[,2]==1&dat.hCalib[,4]==0,])
                      TNcalib <- nrow(dat.hCalib[dat.hCalib[,2]==0&dat.hCalib[,4]==0,])

                      AUC <- auc(cbind((1:nrow(dat.hTest)),dat.hTest[,2],dat.hTest[,3]),na.rm = T)
                      
                      SEN <- TPtest/(FNtest+TPtest)*100
                      ESP <- TNtest/(FPtest+TNtest)*100
                      
                      TSS <- (SEN/100) + (ESP/100) -1
                      
                      Po <- (TPtest + TNtest)/nrow(dat.hTest)#the proportion of units in which the judges agreed
                      Pe <-((TPtest+FPtest)/nrow(dat.hTest))*((TPtest+FNtest)/nrow(dat.hTest))+((TNtest+FPtest)/nrow(dat.hTest))*((TNtest+FNtest)/nrow(dat.hTest))#the proportion of units for which agreement is expected by chance.
                      KAPPA <- (Po - Pe)/(1 - Pe)
                      

                      AUCcalib <- auc(cbind((1:nrow(dat.hCalib)),dat.hCalib[,2],dat.hCalib[,3]),na.rm = T)

                      PoCalib <- (TPcalib + TNcalib)/nrow(dat.hCalib)#the proportion of units in which the judges agreed
                      PeCalib <-((TPcalib+FPcalib)/nrow(dat.hCalib))*((TPcalib+FNcalib)/nrow(dat.hCalib))+((TNcalib+FPcalib)/nrow(dat.hCalib))*((TNcalib+FNcalib)/nrow(dat.hCalib))#the proportion of units for which agreement is expected by chance.
                      KAPPAcalib <- (PoCalib - PeCalib)/(1 - PeCalib)
                      
                      
                      SENcalib <- TPcalib/(FNcalib+TPcalib)*100
                      ESPcalib <- TNcalib/(FPcalib+TNcalib)*100

                      TSScalib <- (SENcalib/100) + (ESPcalib/100) -1
                      
                      OF.AUC <- AUCcalib[1]-AUC[1]
                      
                      OF.TSS <-TSScalib-TSS
                      OF.KAPPA <- KAPPAcalib-KAPPA
                      OF.SEN <- SENcalib-SEN
                      OF.SPE <- ESPcalib-ESP

                      dat.h3 <- data.frame(AUC[1],TSS,KAPPA,SEN,ESP,OF.AUC,OF.TSS,OF.KAPPA,OF.SEN,OF.SPE)

                      dat.h3$spp <- myRespName.h
                      dat.h3$run <- e
                      dat.h3$algo <- as.factor(algos[w])

                      names(dat.h3)[6] <- "OF.AUC"

                      dat.out1 <- rbind(dat.out1,dat.h3)
                    }
                  }
                }


                write.csv(dat.out1,file = paste("./algoRankResults/","modelsEvalFor_",myRespName.h,".csv",sep = ""),row.names = F)

                out.h <- dat.out1[,c(1:10,13)] %>% group_by(algo) %>% summarise(across(everything(), list(mean)))%>%data.frame()

                out.h2 <- data.frame(out.h$algo,sapply(out.h[,2:6],norm10),(-1*sapply(out.h[,7:11],norm10)+1))

                names(out.h2) <- str_remove_all(names(out.h2),"_1")

                names(out.h2)[1] <- "algo"

                out.h3 <- data.frame(out.h2$algo,data.frame(mapply(`*`,out.h2[,metrics],weights)) %>%
                                       mutate(sumWeightsScores = reduce(., `+`)))

                out.h3 <- out.h3[order(out.h3$sumWeightsScores,decreasing = T),]

                names(out.h3)[1] <- "algo"

                write.csv(out.h3,file = paste("./algoRankResults/","rankedEval_",myRespName.h,".csv",sep = ""),row.names = F)

                dat.out2.h <- data.frame(out.h3[,c("algo","sumWeightsScores")][1,],myRespName.h)
                names(dat.out2.h)[3] <- "spp"

                dat.out2 <- rbind(dat.out2,dat.out2.h)

              }


          }
          
            return(dat.out2)
 
          }
        }
      }
  }
