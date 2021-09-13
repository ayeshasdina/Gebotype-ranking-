##### for only fixed time point without TPO #####

rm(list=ls())
ls()
library(HMM)
library(stringr)
library(plyr)
setwd("/Users/dina/Documents/QTL mapping/operation_on_real_data/version_1/HMM_for_real_data_all/input_obs/")



checkEmission <- function(bw)
{
  count_num <- 0
  num_obs <- length(bw$hmm$Symbol)
  emission <- bw$hmm$emissionProbs
  state = bw$hmm$States
  # print(num_obs)
  for(i in c(1:num_obs)) ### going thru all the observation
  {
    if (sum(emission[,i]) == 0)
    {
      emission[state[1],i] = 0.00001 ######## introduced small value to avoid "0" as a emission probability
      emission[state[2],i] = 0.00001
      count_num <- count_num + 1
      
    }
  }
  
  temp <- which.max(emission[state[1],]) ## getting the max index
  emission[state[1],temp] <- emission[state[1],temp] - (0.00001*count_num) #### for unstress row
  temp <- which.max(emission[state[2],]) ## getting the max index
  emission[state[2],temp] <- emission[state[2],temp] - (0.00001*count_num) ### for stressed row
  
  bw$hmm$emissionProbs <- emission
  
  return(bw)
  
}

HMM<-function(input,hmm)
{
  
  # Sequence of observation
  observation = input
  
  # Baum-Welch
  bw = baumWelch(hmm,observation,100,0.01)
  bw1 <- checkEmission(bw)
  return (bw1)
}

Predict_nextObs<-function(genotype_rep,len,dict_S,temp_rep,hmm_model)
{
  
  for (rep in 1:len)
  {
    
    
    
    predicted_obs <- index[genotype_rep[rep],start_time_s1:(start_time_s1+1)] ### giving the 1st two obs from the original data
    match_obs <- c(1,1)
    
    for ( num_to_predict in 1:4) ### seq length to predict
    {
      first_obs <- index[genotype_rep[rep],start_time_s1:(start_time_s1+num_to_predict)] #### g1st two value of obs
      print(first_obs)
      max_likelihood = array()
      
      for (j in obs_array)
      {
        temp_obs <- as.character(c(first_obs,j))
        print(temp_obs)
        Next_obs<- exp(forward(hmm_model,temp_obs))
        
        max_likelihood <- rbind(max_likelihood,sum(Next_obs[,ncol(Next_obs)]))
        
      }
      ########################## Scoring 
      ori_obs <- index[genotype_rep[rep],(start_time_s1+num_to_predict)+1]
      #print(ori_obs)
      if((which.max(max_likelihood) - 1) == ori_obs)
      {
        match_obs <- c(match_obs,1)
      }
      else
      {
        match_obs <- c(match_obs,0)
      }
      
      
      #########################
      #print(max_likelihood)
      
      
      #print(which.max(max_likelihood) - 1)
      predicted_obs <- c(predicted_obs,which.max(max_likelihood) - 1) #### array has a NA value as a 1st element, return the index of the
      # print(as.character(predicted_obs))
    }
    print(match_obs)
    sVal <- sum(match_obs)/6 #### calculate the percentage
    print(sVal)
    # 
    #                
    #                 # print(as.character(predicted_obs))
    #                 temp_predicted_obs <- as.character(predicted_obs)
    #                ########### scoring predicted sequence ###########
    #                 match_obs <- vector()
    #                 for (ori in c(1:length(obs)))
    #                 {
    #                   for (pred in c(ori:ori))
    #                     if (obs[ori]== predicted_obs[pred])
    #                       {flag <- 1
    #                       match_obs <- c(match_obs,flag)}
    #                     else
    #                       {flag <- 0
    #                       match_obs <-c(match_obs,flag) }
    #                 }
    #                 
    #                 sVal <- sum(match_obs)/length(obs)
    #                 
    #                 # print(sVal)
    #                 temp_rep <- c(temp_rep,predicted_obs,sVal)
    temp_rep <- c(temp_rep," ",sVal)
    
    #                 # temp_rep <- c(temp_rep,sVal)
    #                 
    #               
    # 
  }
  # 
  max_value <- max(ncol(dict_S),length(temp_rep))
  # 
  ############ append on the dataframe ######
  if (ncol(dict_S) == max_value)
  {
    
    
    if (length(temp_rep)< max_value)
    {
      num_rep <- max_value - length(temp_rep)
      rep_value <- rep('NA',num_rep)
      temp_rep <- c(temp_rep,rep_value)
      
    }
    
    
    temp_rep <- setNames(temp_rep,c(1:max_value))
    
    colnames(dict_S) <- c(1:max_value)
    dict_S <- rbind(dict_S,temp_rep)
    
  }
  else if(ncol(dict_S) == 0)
  {
    dict_S <- rbind(dict_S,temp_rep)
    colnames(dict_S) <- c(1:max_value)
    
  }
  else
  {
    
    temp_rep <- setNames(temp_rep,c(1:max_value))
    num_rep <- max_value - ncol(dict_S)
    new_col <- c((ncol(dict_S) +1):max_value)
    
    dict_S1[new_col] <- NA
    colnames(dict_S) <- c(1:max_value)
    
    dict_S <- rbind(dict_S,temp_rep)
    
  }
  
  return(dict_S)
}



############# initializing variable ##########
temp_S1 <- vector()
temp_S2 <- vector()
temp_R1 <- vector()
temp_R2 <- vector()
temp_hmm <- vector()


dict_S1 <- data.frame()
dict_S2 <- data.frame()
dict_R1 <- data.frame()
dict_R2 <- data.frame()

inter_s1 <- data.frame(matrix(ncol = 12,nrow = 1))

inter_s2 <- vector()
inter_r1 <- vector()
inter_r2 <- vector()
i <-0

start_time_s1 <- 34
end_time_s1 <- 39

start_time_r1 <- 50 
end_time_r1 <- 55

start_time_s2 <- 60 
end_time_s2 <- 65

start_time_r2 <- 78 
end_time_r2 <- 85

obs_array <- c(1,2,3,4,5)  ## observation name array to calculate the maximum likelihood 


output_path <- "/Users/dina/Documents/QTL mapping/operation_on_real_data/version_1/HMM_for_real_data_all/output_R_transition/HMM_2/"


########## initialization of model ###########


hmm = initHMM(c("Unstressed","Stressed"),c("1","2","3","4","5"),
              transProbs=matrix(c(.5,.5,.5,.5),2),
              emissionProbs=matrix(c(.62,.01,.34,.01,.01,.07,.01,.90,.01,.01),5)) ###emissionProbs=matrix(c(.41,.01,.54,.4,.02,.06,.00001,.5,.03,.03),5))### emissionProbs=matrix(c(.34,.01,.47,.5,.02,.04,.00001,.31,.16999,.14),5)
index <- read.csv(file="New_seq_obsDiscritization.csv", header=TRUE, sep=",")
index_name <- index[,c(1)]
# print(index[,c(1)])



for (plant in index[,c(1)])##for (plant in "125B[set_1][cam_0]") ##for (plant in index[,c(1)])
{
  print(plant)
  i <- i+1 ##37 ##i+1
  # plant_num <- str_extract_all(plant, "^(\\d+)\\[set_.*|^[A-Z]+[0-9]+\\[|^[0-9]+[A-Z]+\\[") 
  plant_num <- str_match(plant, "^(\\w+)\\[") [2]
  print(plant_num)
  plant_rep_name <-c(plant_num)
  genotype_rep = grep(paste("^",plant_num,"\\[",sep = ""),index_name) ### getting replicate of same genotype 
  print(genotype_rep)
  
  len <- length(genotype_rep)
  temp_rep <- array()
  # # temp_rep <- index[i,start_time_s1:end_time_s1]
  # 
  ############# Stressed ###########

  ########### stressed 1 #######
  obs <- vector()


  df = index[i,start_time_s1:end_time_s1]
  input <- as.character(as.vector(t(df)))
  input <- input[!is.na(input)]

  obs <- c(obs,input)
  print (obs)





  hmm_result <- HMM(obs,hmm)

  transition <-(as.array(t(hmm_result$hmm$transProbs)))
  emission <-(as.array(t(hmm_result$hmm$emissionProbs)))
  temp <- c(plant)
  temp_hmm <- rbind(temp_hmm,c(plant,transition," ",emission)) #### intermediate result of HMM models
  temp_S1 <- rbind(temp_S1,temp)
  #print(hmm_result$hmm)

   ################################
  
  #dict_S1 <- Predict_nextObs(genotype_rep,len,dict_S1,temp_rep,hmm_result$hmm)  ## calling the function

  
  # # # index_name_predictedObs <- 0

      for (rep in 1:len)
      {



                  predicted_obs <- index[genotype_rep[rep],start_time_s1:(start_time_s1+1)] ### giving the 1st two obs from the original data
                  match_obs <- c(1,1)

                  for ( num_to_predict in 1:4) ### seq length to predict
                  {
                    first_obs <- index[genotype_rep[rep],start_time_s1:(start_time_s1+num_to_predict)] #### g1st two value of obs
                    #print(first_obs)
                    max_likelihood = array()

                    for (j in obs_array)
                    {
                      temp_obs <- as.character(c(first_obs,j))
                      #print(temp_obs)
                      Next_obs<- exp(forward(hmm_result$hmm, temp_obs))

                      max_likelihood <- rbind(max_likelihood,sum(Next_obs[,ncol(Next_obs)]))

                    }
                    ########################## Scoring
                    ori_obs <- index[genotype_rep[rep],(start_time_s1+num_to_predict)+1]
                    #print(ori_obs)
                    if((which.max(max_likelihood) - 1) == ori_obs)
                    {
                      match_obs <- c(match_obs,1)
                    }
                    else
                    {
                      match_obs <- c(match_obs,0)
                    }


                    #########################
                    #print(max_likelihood)


                    #print(which.max(max_likelihood) - 1)
                    predicted_obs <- c(predicted_obs,which.max(max_likelihood) - 1) #### array has a NA value as a 1st element, return the index of the
                    # print(as.character(predicted_obs))
                  }
                  #print(match_obs)
                  sVal <- sum(match_obs)/6 #### calculate the percentage
                  #print(sVal)
  #
  #
  #                 # print(as.character(predicted_obs))
  #                 temp_predicted_obs <- as.character(predicted_obs)
  #                ########### scoring predicted sequence ###########
  #                 match_obs <- vector()
  #                 for (ori in c(1:length(obs)))
  #                 {
  #                   for (pred in c(ori:ori))
  #                     if (obs[ori]== predicted_obs[pred])
  #                       {flag <- 1
  #                       match_obs <- c(match_obs,flag)}
  #                     else
  #                       {flag <- 0
  #                       match_obs <-c(match_obs,flag) }
  #                 }
  #
  #                 sVal <- sum(match_obs)/length(obs)
  #
  #                 # print(sVal)
  #                 temp_rep <- c(temp_rep,predicted_obs,sVal)
                    temp_rep <- c(temp_rep,sVal)
  #                 # temp_rep <- c(temp_rep,sVal)
  #
  #
  #
      }
  #
  max_value <- max(ncol(dict_S1),length(temp_rep))
  
  # if (plant_rep_name[i-1] == plant_num & plant_rep_name != NULL)
  # {
  #   
  # }
  #
  ############ append on the dataframe ######
  if (ncol(dict_S1) == max_value)
  {


     if (length(temp_rep)< max_value)
     {
       num_rep <- max_value - length(temp_rep)
       rep_value <- rep('NA',num_rep)
       temp_rep <- c(temp_rep,rep_value)

     }


    temp_rep <- setNames(temp_rep,c(1:max_value))

    colnames(dict_S1) <- c(1:max_value)
    dict_S1 <- rbind(dict_S1,temp_rep)

  }
  else if(ncol(dict_S1) == 0)
  {
    dict_S1 <- rbind(dict_S1,temp_rep)
    colnames(dict_S1) <- c(1:max_value)

  }
  else
  {

    temp_rep <- setNames(temp_rep,c(1:max_value))
    num_rep <- max_value - ncol(dict_S1)
    new_col <- c((ncol(dict_S1) +1):max_value)

    dict_S1[new_col] <- NA
    colnames(dict_S1) <- c(1:max_value)

    dict_S1 <- rbind(dict_S1,temp_rep)

  }




  
#######################
  
  
  
  
  # # # # #   ########### stressed 2 #######
  temp_rep <- array()
  obs <- vector()



  df <- index[i,60:65]

  input <- as.character(as.vector(t(df)))
  input <- input[!is.na(input)]

  obs <- c(obs,input)




  hmm_result <- HMM(obs,hmm)                    ##### Calling HMM Function . used a Cran Package and from that I used baum_weltch function to train my data
  transition <-(as.vector(t(hmm_result$hmm$transProbs)))
  temp <- c(plant)
  temp_S2 <- rbind(temp_S2,temp)

  ############


  for (rep in 1:len)
  {



    predicted_obs <- index[genotype_rep[rep],start_time_s2:(start_time_s2+1)] ### giving the 1st two obs from the original data
    match_obs <- c(1,1)

    for ( num_to_predict in 1:4) ### seq length to predict
    {
      first_obs <- index[genotype_rep[rep],start_time_s2:(start_time_s2+num_to_predict)] #### g1st two value of obs
      print(first_obs)
      max_likelihood = array()

      for (j in obs_array)
      {
        temp_obs <- as.character(c(first_obs,j))
        print(temp_obs)
        Next_obs<- exp(forward(hmm_result$hmm, temp_obs))

        max_likelihood <- rbind(max_likelihood,sum(Next_obs[,ncol(Next_obs)]))

      }
      ########################## Scoring
      ori_obs <- index[genotype_rep[rep],(start_time_s2+num_to_predict)+1]
      #print(ori_obs)
      if((which.max(max_likelihood) - 1) == ori_obs)
      {
        match_obs <- c(match_obs,1)
      }
      else
      {
        match_obs <- c(match_obs,0)
      }


      #########################
      #print(max_likelihood)


      #print(which.max(max_likelihood) - 1)
      predicted_obs <- c(predicted_obs,which.max(max_likelihood) - 1) #### array has a NA value as a 1st element, return the index of the
      # print(as.character(predicted_obs))
    }
    print(match_obs)
    sVal <- sum(match_obs)/6 #### calculate the percentage
    print(sVal)
    #
    #
    #                 # print(as.character(predicted_obs))
    #                 temp_predicted_obs <- as.character(predicted_obs)
    #                ########### scoring predicted sequence ###########
    #                 match_obs <- vector()
    #                 for (ori in c(1:length(obs)))
    #                 {
    #                   for (pred in c(ori:ori))
    #                     if (obs[ori]== predicted_obs[pred])
    #                       {flag <- 1
    #                       match_obs <- c(match_obs,flag)}
    #                     else
    #                       {flag <- 0
    #                       match_obs <-c(match_obs,flag) }
    #                 }
    #
    #                 sVal <- sum(match_obs)/length(obs)
    #
    #                 # print(sVal)
                    temp_rep <- c(temp_rep,sVal)
                      #temp_rep <- c(temp_rep," ",predicted_obs,sVal)
    #                 # temp_rep <- c(temp_rep,sVal)
    #
    #
    #
  }
  #
  max_value <- max(ncol(dict_S2),length(temp_rep))
  #
  ############ append on the dataframe ######
  if (ncol(dict_S2) == max_value)
  {


    if (length(temp_rep)< max_value)
    {
      num_rep <- max_value - length(temp_rep)
      rep_value <- rep('NA',num_rep)
      temp_rep <- c(temp_rep,rep_value)

    }


    temp_rep <- setNames(temp_rep,c(1:max_value))

    colnames(dict_S2) <- c(1:max_value)
    dict_S2 <- rbind(dict_S2,temp_rep)

  }
  else if(ncol(dict_S2) == 0)
  {
    dict_S2 <- rbind(dict_S2,temp_rep)
    colnames(dict_S2) <- c(1:max_value)

  }
  else
  {

    temp_rep <- setNames(temp_rep,c(1:max_value))
    num_rep <- max_value - ncol(dict_S2)
    new_col <- c((ncol(dict_S2) +1):max_value)

    dict_S2[new_col] <- NA
    colnames(dict_S2) <- c(1:max_value)

    dict_S2 <- rbind(dict_S2,temp_rep)

  }

  # ############
  # 
  # 
  # # #   
  # # # # #
  # # # # #   ############# Recover ###########
  # # # # #
  # # # # #   ########### recovered 1 #######
  obs <- vector()
  temp_rep <- array()
  #


  df <- index[i,50:55]

  input <- as.character(as.vector(t(df)))
  input <- input[!is.na(input)]

  obs <- c(obs,input)


  hmm_result <- HMM(obs,hmm)
  transition <-(as.vector(t(hmm_result$hmm$transProbs)))
  temp <- c(plant)
  temp_R1 <- rbind(temp_R1,temp)
  ################


  for (rep in 1:len)
  {



    predicted_obs <- index[genotype_rep[rep],start_time_r1:(start_time_r1+1)] ### giving the 1st two obs from the original data
    match_obs <- c(1,1)

    for ( num_to_predict in 1:4) ### seq length to predict
    {
      first_obs <- index[genotype_rep[rep],start_time_r1:(start_time_r1+num_to_predict)] #### g1st two value of obs
      print(first_obs)
      max_likelihood = array()

      for (j in obs_array)
      {
        temp_obs <- as.character(c(first_obs,j))
        print(temp_obs)
        Next_obs<- exp(forward(hmm_result$hmm, temp_obs))

        max_likelihood <- rbind(max_likelihood,sum(Next_obs[,ncol(Next_obs)]))

      }
      ########################## Scoring
      ori_obs <- index[genotype_rep[rep],(start_time_r1+num_to_predict)+1]
      #print(ori_obs)
      if((which.max(max_likelihood) - 1) == ori_obs)
      {
        match_obs <- c(match_obs,1)
      }
      else
      {
        match_obs <- c(match_obs,0)
      }


      #########################
      #print(max_likelihood)


      #print(which.max(max_likelihood) - 1)
      predicted_obs <- c(predicted_obs,which.max(max_likelihood) - 1) #### array has a NA value as a 1st element, return the index of the
      # print(as.character(predicted_obs))
    }
    print(match_obs)
    sVal <- sum(match_obs)/6 #### calculate the percentage
    print(sVal)
    #
    #
    #                 # print(as.character(predicted_obs))
    #                 temp_predicted_obs <- as.character(predicted_obs)
    #                ########### scoring predicted sequence ###########
    #                 match_obs <- vector()
    #                 for (ori in c(1:length(obs)))
    #                 {
    #                   for (pred in c(ori:ori))
    #                     if (obs[ori]== predicted_obs[pred])
    #                       {flag <- 1
    #                       match_obs <- c(match_obs,flag)}
    #                     else
    #                       {flag <- 0
    #                       match_obs <-c(match_obs,flag) }
    #                 }
    #
    #                 sVal <- sum(match_obs)/length(obs)
    #
    #                 # print(sVal)
    #                 temp_rep <- c(temp_rep,predicted_obs,sVal)
                      # temp_rep <- c(temp_rep," ",predicted_obs,sVal)
                      temp_rep <- c(temp_rep,sVal)
    #
    #
    #
  }
  #
  max_value <- max(ncol(dict_R1),length(temp_rep))
  #
  ############ append on the dataframe ######
  if (ncol(dict_R1) == max_value)
  {


    if (length(temp_rep)< max_value)
    {
      num_rep <- max_value - length(temp_rep)
      rep_value <- rep('NA',num_rep)
      temp_rep <- c(temp_rep,rep_value)

    }


    temp_rep <- setNames(temp_rep,c(1:max_value))

    colnames(dict_R1) <- c(1:max_value)
    dict_R1 <- rbind(dict_R1,temp_rep)

  }
  else if(ncol(dict_R1) == 0)
  {
    dict_R1 <- rbind(dict_R1,temp_rep)
    colnames(dict_R1) <- c(1:max_value)

  }
  else
  {

    temp_rep <- setNames(temp_rep,c(1:max_value))
    num_rep <- max_value - ncol(dict_R1)
    new_col <- c((ncol(dict_R1) +1):max_value)

    dict_R1[new_col] <- NA
    colnames(dict_R1) <- c(1:max_value)

    dict_R1 <- rbind(dict_R1,temp_rep)

  }
  # 
  # ###############
  # # 
  # # first_obs <- index[i,50:51]
  # # 
  # # for ( num_to_predict in 1:4)
  # # {
  # #   max_likelihood = array()
  # #   
  # #   for (j in obs_array)
  # #   {
  # #     temp_obs <- as.character(c(first_obs,j))
  # #     print(temp_obs)
  # #     Next_obs<- exp(forward(hmm_result$hmm, temp_obs))
  # #     
  # #     max_likelihood <- rbind(max_likelihood,sum(Next_obs[,ncol(Next_obs)]))
  # #     
  # #   }
  # #   
  # #   # print(max_likelihood)
  # #   first_obs <- c(first_obs,which.max(max_likelihood) - 1) #### array has a NA value as a 1st element
  # #   
  # # }
  # # # print(as.character(first_obs))
  # # # 
  # # # #print(temp_in)
  # # 
  # # temp_in <-c(temp,index[i,50:55],setNames(first_obs,c(1:6)))
  # # # temp_in <-c(temp,index[i,50:55])
  # # dict_R1 <- rbind(dict_R1,temp_in)
  # # # dict_R1 <- rbind(dict_R1,temp)
  # # 
  # # # #
  # # # # 
  # # # 
  # # # # #   ########### recovered 2 #######
  obs <- vector()
  temp_rep <- array()

  df <- index[i,78:85]
  input <- as.character(as.vector(t(df))) ### need some work in there

  input <- input[!is.na(input)]

  obs <- c(obs,input)

  print (obs)

  hmm_result <- HMM(obs,hmm)
  transition <-(as.vector(t(hmm_result$hmm$transProbs)))
  temp <- c(plant)
  temp_R2 <- rbind(temp_R2,temp)
  # ###########


  for (rep in 1:len)
  {



    predicted_obs <- index[genotype_rep[rep],start_time_r2:(start_time_r2+1)] ### giving the 1st two obs from the original data
    match_obs <- c(1,1)

    for ( num_to_predict in 1:6) ### seq length to predict for recover 2 which has 8 timepoints
    {
      first_obs <- index[genotype_rep[rep],start_time_r2:(start_time_r2+num_to_predict)] #### g1st two value of obs
      print(first_obs)
      max_likelihood = array()

      for (j in obs_array)
      {
        temp_obs <- as.character(c(first_obs,j))
        print(temp_obs)
        Next_obs<- exp(forward(hmm_result$hmm, temp_obs))

        max_likelihood <- rbind(max_likelihood,sum(Next_obs[,ncol(Next_obs)]))

      }
      ########################## Scoring
      ori_obs <- index[genotype_rep[rep],(start_time_r2+num_to_predict)+1]
      #print(ori_obs)
      if((which.max(max_likelihood) - 1) == ori_obs)
      {
        match_obs <- c(match_obs,1)
      }
      else
      {
        match_obs <- c(match_obs,0)
      }


      #########################
      #print(max_likelihood)


      #print(which.max(max_likelihood) - 1)
      predicted_obs <- c(predicted_obs,which.max(max_likelihood) - 1) #### array has a NA value as a 1st element, return the index of the
      # print(as.character(predicted_obs))
    }
    print(match_obs)
    sVal <- sum(match_obs)/8 #### calculate the percentage
    print(sVal)
    #
    #
    #                 # print(as.character(predicted_obs))
    #                 temp_predicted_obs <- as.character(predicted_obs)
    #                ########### scoring predicted sequence ###########
    #                 match_obs <- vector()
    #                 for (ori in c(1:length(obs)))
    #                 {
    #                   for (pred in c(ori:ori))
    #                     if (obs[ori]== predicted_obs[pred])
    #                       {flag <- 1
    #                       match_obs <- c(match_obs,flag)}
    #                     else
    #                       {flag <- 0
    #                       match_obs <-c(match_obs,flag) }
    #                 }
    #
    #                 sVal <- sum(match_obs)/length(obs)
    #
    #                 # print(sVal)
    #                 temp_rep <- c(temp_rep,predicted_obs,sVal)
                      # temp_rep <- c(temp_rep," ",predicted_obs,sVal)
                    temp_rep <- c(temp_rep,sVal)
    #
    #
    #
  }
  #
  max_value <- max(ncol(dict_R2),length(temp_rep))
  #
  ############ append on the dataframe ######
  if (ncol(dict_R2) == max_value)
  {


    if (length(temp_rep)< max_value)
    {
      num_rep <- max_value - length(temp_rep)
      rep_value <- rep('NA',num_rep)
      temp_rep <- c(temp_rep,rep_value)

    }


    temp_rep <- setNames(temp_rep,c(1:max_value))

    colnames(dict_R2) <- c(1:max_value)
    dict_R2 <- rbind(dict_R2,temp_rep)

  }
  else if(ncol(dict_R2) == 0)
  {
    dict_R2 <- rbind(dict_R2,temp_rep)
    colnames(dict_R2) <- c(1:max_value)

  }
  else
  {

    temp_rep <- setNames(temp_rep,c(1:max_value))
    num_rep <- max_value - ncol(dict_R2)
    new_col <- c((ncol(dict_R2) +1):max_value)

    dict_R2[new_col] <- NA
    colnames(dict_R2) <- c(1:max_value)

    dict_R2 <- rbind(dict_R2,temp_rep)

  }

  ##############
  #

}

######## binding two dataframe ####
dict_S1[1] <- NULL ### remove the 1st element of the dict_s1
colnames(dict_S1) <- c(1:ncol(dict_S1)) ### set the name again
temp_S1 <- as.data.frame(temp_S1) ### convert to dataframe
dict_S1 <-cbind(temp_S1,dict_S1) #### binding two dataframe

dict_S2[1] <- NULL ### remove the 1st element of the dict_s1
colnames(dict_S2) <- c(1:ncol(dict_S2)) ### set the name again
temp_S2 <- as.data.frame(temp_S2) ### convert to dataframe
dict_S2 <-cbind(temp_S2,dict_S2) #### binding two dataframe
# 
# 
dict_R1[1] <- NULL ### remove the 1st element of the dict_s1
colnames(dict_R1) <- c(1:ncol(dict_R1)) ### set the name again
temp_R1 <- as.data.frame(temp_R1) ### convert to dataframe
dict_R1 <-cbind(temp_R1,dict_R1) #### binding two dataframe
# 
# 
dict_R2[1] <- NULL ### remove the 1st element of the dict_s1
colnames(dict_R2) <- c(1:ncol(dict_R2)) ### set the name again
temp_R2 <- as.data.frame(temp_R2) ### convert to dataframe
dict_R2 <-cbind(temp_R2,dict_R2) #### binding two dataframe


# print(dict_S1)
# # #
# write.table(dict_S1,file=paste(output_path,"transition_s1.csv"),col.names = c("plantname","U2U","U2S","S2U","S2S"),sep=",",row.names = FALSE)
# write.table(dict_S2,file=paste(output_path,"transition_s2.csv"),col.names = c("plantname","U2U","U2S","S2U","S2S"),sep=",",row.names = FALSE)
# 
# write.table(dict_R1,file=paste(output_path,"transition_r1.csv"),col.names = c("plantname","U2U","U2S","S2U","S2S"),sep=",",row.names = FALSE)
# write.table(dict_R2,file=paste(output_path,"transition_r2.csv"),col.names = c("plantname","U2U","U2S","S2U","S2S"),sep=",",row.names = FALSE)
#temp_hmm <- set_names(temp_hmm,c("plantname","U2U","U2S","S2U","S2S","","U-1","U-2","U-3","U-4","U-5","S-1","S-2","S-3","S-4","S-5"))

write.csv(dict_S1,file=paste(output_path,"transition_s1_performance.csv"),sep=",",row.names = FALSE)
# write.csv(temp_hmm,file=paste(output_path,"hmm_S1.csv"),sep=",",row.names = FALSE)
write.csv(dict_S2,file=paste(output_path,"transition_s2_performance.csv"),sep=",",row.names = FALSE)
# # 
write.csv(dict_R1,file=paste(output_path,"transition_r1_performance.csv"),sep=",",row.names = FALSE)
write.csv(dict_R2,file=paste(output_path,"transition_r2_performance.csv"),sep=",",row.names = FALSE)

