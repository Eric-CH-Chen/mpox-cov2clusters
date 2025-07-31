#!/usr/bin/env -S Rscript 
#Importation from Nextflow JMC 2023
#
# Adapted from Jessica Caleta's covflow
#
# Note: the meta file needs to have 'StrainID' and 'Date' as the header for strain and date
#
#   On reusing old cluster names:
#   It requires "clusterFile", "pastTransProbs", "newClustering=FALSE" to be set
#
#

args = commandArgs(trailingOnly=TRUE)
input_tree <- args[1] #just placeholder in funct definition, but default if ordered not specified
input_meta <- args[2]
input_beta <- eval(parse(text=args[3]))  # so it evaluates the "c(1, 2, 3)"  string and creates the vector of beta.
input_probTh <- args[4]
input_fstring <- args[5]
input_label_string <- args[6]
input_returnProb <- as.logical(args[7]) # whether to return Prob calculation of cov2clust
test_run_tag <- args[8]

# input_probTh is likely to be coerced into a character value, so this checks and converts it back to numeric
input_probTh <- ifelse(class(input_probTh) == "character", eval(parse(text=input_probTh)), eval(parse(text=paste0("c(", input_probTh, ")"))) )

#---------------------
# Function `TransProbs` Replacement via data.table
#---------------------

#' Calculate pairwise logit probability between tree sequences
#'
#' @export
#' 
TransProbs <- function(inputData, contactData, probThreshold, dateThreshold, restrictClusters, beta) {
#TransProbs <- function(PatDist, date_df, contactData, probThreshold, dateThreshold, restrictClusters, beta) {
  # contactData and restrictClusters are ignored in this iteration, currently it will simply crash via stop
  if(restrictClusters | contactData ) {
    stop("Option restricClusters and contactData is not implemented in this version")
  }
  
  PatDist <- inputData$PatDist
  date_df <- inputData$dates
  date_df$Date <- as.Date(date_df$Date)
  #write.table(PatDist,paste0("_PatDist_",Sys.Date(),".txt"),quote = F,sep = "\t",row.names = F)
  # loop 1 - Date differences
  dt <- as.data.table(date_df)
  dt[, row_id := .I]  # add row number for sorting
  
  result <- dt[dt, on = .(row_id < row_id), allow.cartesian = TRUE, nomatch = 0]
  result <- result[, .(strain_1 = StrainID, strain_2 = i.StrainID, date_diff = as.numeric(abs(Date - i.Date)))]
    
  
  # loop 2
  dt_distance <- data.table(as.table(PatDist), 
                            stringsAsFactors = FALSE
  )
  # loop 3
  dt_com <- merge(dt_distance, result, all = FALSE, 
                  by.y = c("strain_1", "strain_2"), 
                  by.x = c("V1", "V2"))
  names(dt_com)[names(dt_com) == 'N'] <- 'PatD' 
  
  # this isn't needed if it is merge(result, dt_distance, all=false)
  names(dt_com)[names(dt_com) == 'V1'] <- 'strain_1' 
  names(dt_com)[names(dt_com) == 'V2'] <- 'strain_2' 
  
  # loop 4
  # cov2cluster logit model
  # trans = 1/(1 + exp(- (sum(beta*covariates))))
  
  dt_com <- dt_com[, Prob :=  ( 1 / (1 + exp(- (beta[1] + beta[2] * PatD + beta[3] * date_diff))) )]
  
  # loop 5
  # filter according to date threshold
#  dt_com <- dt_com[date_diff <= dateThreshold & Prob > pminProb]
  dt_com <- dt_com[, Prob:= ifelse(date_diff <= dateThreshold, Prob, 0)]
  dt_com$date_diff <- NULL
  dt_com$PatD <- NULL
  
  # loop 6
  # convert and return as data.frame object so it can fit back
  # 
  df_com <- as.data.frame(dt_com)
  colnames(df_com)<-c("host1","host2","Prob")
  #write.table(dt_com,paste0("_TransProbsMid_",Sys.Date(),".txt"),quote = F,sep = "\t",row.names = F)
  return(df_com)
  
}
#---------------------
# Function `numberClusters`
#---------------------

#' Hierarchical clustering of sequences based on accepted transmission links above probability threshold
#' EC: Unchanged from cov2cluster
#' @export

numberClusters<-function(acceptTrans){
  names<-unique(c(as.character(acceptTrans[,1]),as.character(acceptTrans[,2])))   # EC get list of unique names, flattens the 'acceptTrans' two columns, [,1] and [,2], into one column
  cluster_results<-matrix(NA,ncol = 2, nrow = length(names))                      # EC initialize matrix with 2-col and number of rows equal to length of unique names
  cluster_results[,1]<-names      # EC fill first column (V1) with unique names
  if (nrow(cluster_results)>0){   # EC This is a check to ensure cluster results(?) are tallied only if there are at least 1
    nums<-1:nrow(cluster_results) # EC get index of rows (ie. 1, 2, 3, 4)
    for (i in 1:nrow(cluster_results)){   # EC loop over all rows of the matrix cluster_results
      if (i%in%nums){
        clustnums<-which(cluster_results[,1] %in% unique(unlist(c(acceptTrans[which(acceptTrans[,1]==cluster_results[i,1] | 
                                                                                      acceptTrans[,2]==cluster_results[i,1]),1:2]))))
        assocnums<-which(cluster_results[,1] %in% unique(unlist(c(acceptTrans[which(acceptTrans[,1] %in% cluster_results[clustnums,1] | 
                                                                                      acceptTrans[,2] %in% cluster_results[clustnums,1]),1:2]))))
        allnums<-c(clustnums,assocnums)
        prevClust<-unique(cluster_results[allnums,2])
        prevClust<-as.numeric(prevClust[!is.na(prevClust)])   
        if (length(prevClust)>0){
          cluster_results[unique(c(which(cluster_results[,2]%in%prevClust),allnums)),2]<-min(prevClust)
        } else{
          cluster_results[allnums,2]<-min(nums) #
        }
        nums<-nums[!nums%in%clustnums]
      }
    }
  }
  return(cluster_results)
}

#---------------------
# Function `mergeOldRun`
#---------------------
#' Replace Prob from new run with previously ran data
mergeOldRun <- function(newProb, oldProb) {
  newProb <- as.data.table(newProb)
  oldProb <- as.data.table(oldProb)
  
  setkey(newProb, host1, host2)
  setkey(oldProb, host1, host2)
  merge_data <- merge(newProb, oldProb, by=c("host1", "host2"), all.x=TRUE, suffixes = c("", "old"))
  merge_data[, Prob := ifelse(!is.na(Probold), Probold, Prob)] 
  merge_data[, Probold := NULL]
  
  # Return in data.frame for consistency
  return(as.data.frame(merge_data))
}

#---------------------
# Function `cov2cluster`
#---------------------
#' updated to work with new function

#' @param treeName Name of the input phylogenetic tree, in newick format
#' @param metafile Name of the date file if in .csv format with columns "sequence ID" and "dates" (not required if json_dates = TRUE and restrictClusters = FALSE)
#' @param json_dates If TRUE - dates are supplied in json format file
#' @param json_file If json_dates = TRUE - name .json dates file
#' @param contactData If TRUE - file will be supplied with contact data 
#' @param contactFile File contacting contact, in .csv format (required if contactData = TRUE)
#' @param beta List of the beta coefficients beta0 - beta3 (Intercept, patristic distance, date difference, and contact data(if applicable))
#' @param returnTransProbs Return text file of pairwise logit probabilities (must = TRUE if newClustering = FALSE)
#' @param dateThreshold Integer for the hard upper limit of date difference (in days) to link pairs of sequences in a cluster (default = 40)
#' @param restrictClusters If TRUE - Cluster sequences only with a shared variable in the metafile column 'restrictCluster' (optional)
#' @param probThreshold Integer or list - Pairwise probability threshold for clustering sequences
#' @param newClustering If TRUE - perform clustering on all sequences for the first time, if FALSE - include past cluster designations of sequences in new clustering run
#' @param pastTransProbs If newClustering = FALSE, this will be the text file returned from returnTransProbs = TRUE in the past clustering run
#' @param clusterFile If newClustering = FALSE, this will be output cluster file returned from the past clustering run
#' @param clusternameIdent Character string to include in the cluster names (default = "clust")
#' @param outfile Prefix for the output file
#' @return Text file with three columns - Sequence ID from tree tip, cluster name, and past cluster name (if applicable)
#' @export

### EC example runs
# example run dir 
#   setwd("~/Projects/ncov/data")
# cov2clusters(treeName = "nwk/tree.nwk", outfile = "covid-test",json_dates = FALSE, metafile = "cluster/SARS-CoV-2_0.8_Test.perl-fakeday.csv",returnTransProbs=TRUE)
# cov2clusters(treeName = "nwk/tree.nwk", outfile = "g-covid-test-reuse",json_dates = FALSE, metafile = "cluster/SARS-CoV-2_0.8_Test.perl-fakeday.csv", newClustering = FALSE, clusterFile = "cluster/SARS-CoV-2_0.8_genomic-clusters_Run-1-2889_Tree_samples.tsv", pastTransProbs = "covid-test_TransProbs_2024-04-10.txt", clusternameIdent="BC")
# cov2clusters(treeName = "../../../runs/example_mpox/run_no_mask_dataset_tree.nwk", outfile = "no_mask_mpox_example",json_dates = FALSE, metafile = "../Raw_noMask_dataset_meta_date.csv",returnTransProbs=TRUE, clusternameIdent = 'noMask')


cov2clusters<-function(treeName="tree.nwk",metafile=NA,
                       json_dates=TRUE,json_file="branch_lengths.json",
                       contactData=FALSE,contactFile=NA,
                       beta=c(3,-19735.98,-0.075,-0.2),returnTransProbs=FALSE,
                       dateThreshold=40,restrictClusters=FALSE,
                       probThreshold=c(0.8,0.9),
                       newClustering=TRUE, pastTransProbs=NA,
                       clusterFile=NA,
                       clusternameIdent = "clust",
                       outfile="SARS-CoV-2",no.Cores=1){
  # Load packages
  library(ape)
  library(reshape2)
  library(rjson)
  library(stringi)
  library(data.table) # EC new
  options(stringsAsFactors = F)  ## EC option?
  
  #Data Setup
  dataInput<-list()
  
  #Patristic distance
  ## EC: Using ape, calculate patristic distance between all leaf and then sort
  tree<-read.tree(treeName)
  dataInput$PatDist<-as.matrix(cophenetic.phylo(tree))
  dataInput$PatDist<-dataInput$PatDist[order(row.names(dataInput$PatDist)),order(colnames(dataInput$PatDist))]
  #Load Metadata
  if (json_dates){
    json_data <- fromJSON(file=json_file)
    dates<-data.frame(names=colnames(dataInput$PatDist),date=NA)
    for (i in 1:nrow(dates)){
      dates$date[i]<-as.character(unlist(json_data$nodes[[which(names(json_data$nodes)==dates$names[i])]][3]))
    }
  } else {
    dates <- read.csv(metafile,check.names = F) ## EC load dates csv if JSON not present
  }
  dates<-dates[which(dates[,1] %in% colnames(dataInput$PatDist)),]
  dates<-dates[order(dates[,1]),]
  # remove patristic distance without dates
  dataInput$PatDist<-dataInput$PatDist[which(colnames(dataInput$PatDist) %in% dates[,1]),
                                       which(colnames(dataInput$PatDist) %in% dates[,1])]
  # EC commented out to preserve the date obj/ and just save the df directly
  #dataInput$dates<-round(lubridate::decimal_date(as.Date(dates[,2]))*365)
  dataInput$dates <- dates

  
  if (restrictClusters){
    variable <- read.csv(metafile,check.names = F)
    variable<-variable[which(variable[,1] %in% colnames(dataInput$PatDist)),]
    variable<-variable[order(variable[,1]),]
    dataInput$variable<-variable$restrictCluster
  }
  # Contact data
  if (contactData){
    contacts<-read.table(contactFile)
    contacts<-contacts[which(row.names(contacts) %in% colnames(dataInput$PatDist)),
                       which(colnames(contacts) %in% colnames(dataInput$PatDist))]
    contacts<-contacts[order(colnames(contacts)),order(colnames(contacts))]
    dataInput$contacts<-contacts
  }
  
  # Transmission matrix
  transmat<-TransProbs(dataInput,contactData,probThreshold,dateThreshold,restrictClusters,beta)
  # if (!newClustering){
  #  pastTransProbs<-read.table(pastTransProbs,header = T,check.names = F)
  #  newnames<-colnames(dataInput$PatDist)[-which(colnames(dataInput$PatDist) %in% 
  #                                                pastTransProbs$host1 |
  #                                                 colnames(dataInput$PatDist) %in% 
  #                                                 pastTransProbs$host2)]   # EC Get names of samples from dataInput that does not 
  #                                                                          #   exists in the previous ran pastTransProbs
  #  transmat<-transmat[which(transmat[,1] %in% newnames | 
  #                             transmat[,2] %in% newnames),]    # EC Filter transmat by keeping only samples that did not exist in
  #                                                              #   pastTransProbs via newnames
  #  if (nrow(transmat)>0){
  #    transmat<-rbind(transmat,pastTransProbs)  # EC combine new pairs and previous run.
  #  } else {
  #    transmat<-pastTransProbs  # EC in case no new data is added, just use the previous run.
  #  }
  #}
  # EC new way of writing the clust
  if(!newClustering) {
    pastTransProbs<-read.table(pastTransProbs,header = T,check.names = F)
    transmat <- mergeOldRun(transmat, pastTransProbs)
  }

  if (returnTransProbs){
    write.table(transmat,paste0(outfile,"_TransProbs_",Sys.Date(),".txt"),quote = F,sep = "\t",row.names = F)
  }
  
  # Cluster using defined thresholds 
  for (threshold in 1:length(probThreshold)){
    acceptTrans<-transmat[which(transmat[,3]>=probThreshold[threshold]),]
    
    # Number clusters
    cluster_results<-numberClusters(acceptTrans)
    # EC modified the below
    cluster_results<-data.frame(cluster_results,pastClustering=NA) #EC original
    #cluster_results<-data.frame(cluster_results)
    #cluster_results$pastClustering <- NA
    
    ## Rename clusters by past cluster names or by earliest case date and location if new clustering
    if (newClustering){
      clusters<-unique(cluster_results[,2])
      months<-c("JAN","FEB","MAR","APR","MAY","JUN","JUL","AUG","SEP","OCT","NOV","DEC")
      randstring<-stri_rand_strings(length(clusters)*2,3,"[b-df-hj-np-tv-z]")
      randstring<-unique(randstring) # make sure string is unique
      for (clust in 1:length(clusters)){
        mon<-month(min(dates[which(dates[,1] %in% 
                                     cluster_results[which(cluster_results[,2]==clusters[clust]),1]),2]))
        yea<-substr(year(min(dates[which(dates[,1] %in% 
                                           cluster_results[which(cluster_results[,2]==clusters[clust]),1]),2])),3,4)
        cluster_results[which(cluster_results[,2]==clusters[clust]),2]<-paste0(paste0(months[mon],yea),".",clusternameIdent,".",randstring[clust])
      }
    } else {
      pastClusters<-read.table(clusterFile,header = T,check.names = F)
      clusters<-unique(cluster_results[,2])
      clustersDF<-data.frame(ClusterName=rep(NA,length(clusters)),ClusterSize=NA,ClusterComposition=rep(NA,length(clusters)))
      months<-c("JAN","FEB","MAR","APR","MAY","JUN","JUL","AUG","SEP","OCT","NOV","DEC")
      randstring<-stri_rand_strings(length(clusters)*3,3,"[b-df-hj-np-tv-z]")
      randstring<-unique(randstring) # make sure string is unique
      pastClusterString<-stri_sub(unique(pastClusters$Cluster_no[!pastClusters$Cluster_no=="-1"]),-3)
      randstring<-randstring[!randstring%in%pastClusterString] # not in past cluster string
      for (clust in 1:length(clusters)){
        rownum<-which(cluster_results[,2]==clusters[clust])
        clustersDF$ClusterSize[clust]<-length(rownum)
        pastClusterNames<-pastClusters[which(pastClusters[,1] %in% cluster_results[rownum,1]),2]
        for (j in 1:length(rownum)){
          if (cluster_results[rownum[j],1] %in% pastClusters[,1]){
            cluster_results[rownum[j],3]<-pastClusters[which(pastClusters[,1] %in% cluster_results[rownum[j],1]),2]
          }}
        if (length(pastClusterNames)>0){
          pastClusterDF<-data.frame(table(pastClusterNames),prop.in.clust=0)
          for (i in 1:nrow(pastClusterDF)){
            if (pastClusterDF$pastClusterNames[i]=="-1"){
              pastClusterDF$prop.in.clust[i]<-NA
            } else {
              pastClusterDF$prop.in.clust[i]<-round(pastClusterDF$Freq[i]/length(which(pastClusters[,2] %in% pastClusterDF$pastClusterNames[i])),2)
            }}
          clustersDF$ClusterComposition[clust]<-paste0(paste0(pastClusterDF$pastClusterNames,",",pastClusterDF$Freq,",",pastClusterDF$prop.in.clust),":",collapse = "")
          pastClusterDF<-pastClusterDF[pastClusterDF$pastClusterNames!="-1",]
          if (nrow(pastClusterDF)>0 &&
              pastClusterDF[which.max(pastClusterDF$Freq),3]>=0.6){
            clustername<-as.character(pastClusterDF[which.max(pastClusterDF$Freq),1])
            cluster_results[rownum,2]<-clustername
            clustersDF$ClusterName[clust]<-clustername
          } else {
            mon<-month(min(dates[which(dates[,1] %in% 
                                         cluster_results[rownum,1]),2]))
            yea<-substr(year(min(dates[which(dates[,1] %in% 
                                               cluster_results[rownum,1]),2])),3,4)
            clustername<-paste0(paste0(months[mon],yea),".",clusternameIdent,".",randstring[clust])
            cluster_results[rownum,2]<-clustername
            clustersDF$ClusterName[clust]<-clustername
          }
        } else {
          mon<-month(min(dates[which(dates[,1] %in% 
                                       cluster_results[rownum,1]),2]))
          yea<-substr(year(min(dates[which(dates[,1] %in% 
                                             cluster_results[rownum,1]),2])),3,4)
          clustername<-paste0(paste0(months[mon],yea),".",clusternameIdent,".",randstring[clust])
          cluster_results[rownum,2]<-clustername
          clustersDF$ClusterName[clust]<-clustername
        }
      }
      write.csv(clustersDF,paste0(outfile,"_",probThreshold[threshold],"_ClustersSummary",Sys.Date(),".csv"),row.names = F)# cluster summary output
    }
    
    # Output file 
    cluster_results<-cluster_results[order(cluster_results[,2]),]
    row.names(cluster_results)<-NULL
    colnames(cluster_results)<-c("SampleID","Cluster.No","pastCluster.No")
    if (!newClustering){
      pastnoclust<-cbind(as.character(pastClusters[which(!as.character(pastClusters$SampleID) %in% 
                                                           cluster_results[,1]),1]),"-1",
                         as.character(pastClusters[which(!as.character(pastClusters$SampleID) %in% 
                                                           cluster_results[,1]),2]))
          ## EC:pastClusters coming from a previously made file, specified in the
          ##  input "clusterFile", but this is only for printing output.
      if (ncol(pastnoclust)==3){
        cluster_results<-rbind(as.matrix(cluster_results),pastnoclust)
      }
    }
    noclust<-cbind(as.character(dates[which(!as.character(dates[,1]) %in% cluster_results[,1]),1]),"-1",NA) # find 
    if (ncol(noclust)==3){
      cluster_results<-rbind(as.matrix(cluster_results),noclust)
    }
    # EC write the cluster_results into a tsv file.
    write.table(cluster_results,paste0(outfile,"_",probThreshold[threshold],"_GenomicClusters",Sys.Date(),".txt"),row.names = F,sep = "\t",quote = F)
  }
}


#---------------------
# Run cov2cluster
#---------------------
# Example to run independantly
#cov2clusters(treeName = input_tree_order,returnTransProbs=TRUE,pastTransProbs=input_past_trans_probs,
#             probThreshold = 0.8,newClustering = FALSE,
#             clusterFile = input_past_gen_clusts,clusternameIdent = "BC",no.Cores = 32)

# implmentation using metaData instead of json
cov2clusters(treeName = input_tree, 
      beta = input_beta,
			metafile = input_meta, 
      json_dates = FALSE, 
			outfile = input_fstring, 
      clusternameIdent = input_label_string, 
			returnTransProbs = input_returnProb, 
      probThreshold = input_probTh)


