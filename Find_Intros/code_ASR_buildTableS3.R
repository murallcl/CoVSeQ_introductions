


# Murall et al. 2021 Genome Medicine 

# Carmen Lia Murall and Sana Naderi
# McGill University, Genome Center, CovSeQ consortium
# September 2021





#----------------------------------- Part 1 : Find introductions and build table

################################# load libraries 
## Loading libraries
library(phylotools)
library(tidyverse)
library(tidytree)
library(ggtree)
library(mvSLOUCH)
library(phangorn)
library(phytools)
library(tictoc)
library(ggimage)
library(gridExtra)
library(treeio)
library(viridis)
library(RColorBrewer)
library(colorspace)
library(data.table)
library(RJSONIO)
library(matrixStats)
################################# 


################################# data loading and preprocessing starts here

#set working directory 
setwd("./path/folder") # ADD path to folder that contains the files loaded in this section

#load tree
tree = read.tree(file = "freeze2_globaltree.nwk")

#load processed tree in searchable table (generated with branch lengths table provided by Nextstrain output)
edge_table<-read.csv("edge_tabe_freeze2_globaltree.csv")

#load metadata 
metadata <-read.csv("metadataGlobal.csv", header = T, sep = ",") #n = 15,722 
metadataQC <- read.csv("metadataQC.csv", header = T, sep = ",")
as.Date(metadata$date) -> metadata$date
as.Date(metadataQC$date) -> metadataQC$date

#table for country/region matching
match <- read.csv("continents_list_regions.csv", header = T, sep = ",")

################################# 



##############################################################################################################
##################### 1. -- Find Intros: using 3 methods (ML, ACC, DEL) -- ###################################


## --- A. Reconstruction ------------------------------------------------------------------------------------

tic()

################################# MPR reconstruction 
# Whole-tree MPR reconstruction, fitch.mvsl function is used from "mvslouch" package.
# both 'deltran' and 'acctran' settinge will be applied.

dichotomous_tree = multi2di(tree, random = FALSE)   #format change of tree
dichotomous_tree$edge.length[dichotomous_tree$edge.length==0] = dichotomous_tree$edge.length[dichotomous_tree$edge.length==0] + 1e-8  #set very small values to zero (needed for ace to run properly)
i = str_detect(dichotomous_tree$tip.label , "Qc-")
QC_indices = which(i == "TRUE")
mpr_tipstates = rep(0, length(dichotomous_tree$tip.label))
mpr_tipstates[QC_indices] = 1
dichotomous_tree$tip.label -> names(mpr_tipstates)  #add names

mpr_recon_deltran = fitch.mvsl(dichotomous_tree, mpr_tipstates, deltran = TRUE, acctran = FALSE, root = 0)
mpr_recon_acctran = fitch.mvsl(dichotomous_tree, mpr_tipstates, deltran = FALSE, acctran = TRUE, root = 0)
################################# 


################################# ML reconstruction 
# labels vector format is same as the MPR analysis, the same vector was used
ml_recon = ace(mpr_tipstates, dichotomous_tree, type = "discrete", model="ER")
#################################


#### NOTE:
# the output of the mpr function corresponds to the edge attribute of the tree
# the output of the ace fucntion increments from the top of the tree (starting from the root)


comparison = comparePhylo(tree, dichotomous_tree) #compares the original polytomous tree to the dichotomous tree
comp = matrix(0, dim(comparison$NODES)[1], dim(comparison$NODES)[2])

for (i in 1:dim(comparison$NODES)[1]){ #loop for extracting the node number out of the compare$node table returned by comparePhylo()
  for (j in 1:dim(comparison$NODES)[2]){
    string = comparison$NODES[i,j]
    str = stringr::str_extract(string = comparison$NODES[i,j], pattern = "(?<=\\().*(?=\\))")
    comp[i,j] = as.numeric(str) #replacing the string number to numeric value, optional, might have to change later
  }
}

################################# cleaning up ace output 
ml_asr = ml_recon$lik.anc
ml_asr_polytomous = matrix(0, dim(comp)[1], 3)
for (i in 1:dim(comp)[1]){
  row = comp[i,2] - length(dichotomous_tree$tip.label)
  ml_asr_polytomous[i,1] = ml_asr[row,1]
  ml_asr_polytomous[i,2] = ml_asr[row,2]
  ml_asr_polytomous[i,3] = i + length(tree$tip.label)
}
ace_pies_dataframe = data.frame(zero = ml_asr_polytomous[,1] , one = ml_asr_polytomous[,2] , node = ml_asr_polytomous[,3])
# pies dataframe is a #internal_nodes by 3 matrix, the first two columns are the ASR estimates and the third column is the node number in the original tree

ace_dataframe = matrix(0, dim(comp)[1], 3)
for (i in 1:dim(ace_pies_dataframe)[1]){
  if(ace_pies_dataframe[i,1] >= ace_pies_dataframe[i,2]){
    ace_dataframe[i,1] = 1
  }

  if(ace_pies_dataframe[i,1] < ace_pies_dataframe[i,2]){
    ace_dataframe[i,2] = 1
  }
}
## the final output of the ace reconstruction. only internal nodes are included.
ace_dataframe = data.frame(zero = ace_dataframe[,1], one = ace_dataframe[,2], node = ace_pies_dataframe[,3])
################################# 

################################# cleaning up MPR output 
## Deltran
mpr_polytomous = matrix(0, dim(comp)[1], 2)
mpr_polytomous[1,1] = 0
mpr_polytomous[1,2] = as.numeric(comp[1,2])
for (i in 2:dim(comp)[1]){
  n_di = comp[i,2]
  in_edge = which(dichotomous_tree$edge[,2] == n_di)
  mpr_polytomous[i,1] = as.numeric(mpr_recon_deltran$branch_regimes[in_edge])
  mpr_polytomous[i,2] = as.numeric(comp[i,1])
}
mpr_del_dataframe = data.frame(zero = 1 - mpr_polytomous[,1] ,one = mpr_polytomous[,1], node = mpr_polytomous[,2])


## ACCtran
mpr_polytomous = matrix(0, dim(comp)[1], 2)
mpr_polytomous[1,1] = 0
mpr_polytomous[1,2] = as.numeric(comp[1,2])
for (i in 2:dim(comp)[1]){
  n_di = comp[i,2]
  in_edge = which(dichotomous_tree$edge[,2] == n_di)
  if(mpr_recon_acctran$branch_regimes[in_edge] == "ambiguous"){
    mpr_polytomous[i,1] = 0.5
    mpr_polytomous[i,2] = as.numeric(comp[i,1])
  }
  else{
    mpr_polytomous[i,1] = as.numeric(mpr_recon_acctran$branch_regimes[in_edge])
    mpr_polytomous[i,2] = as.numeric(comp[i,1])
  }
}
mpr_acc_dataframe = data.frame(zero = 1 -mpr_polytomous[,1] ,one = mpr_polytomous[,1], node = mpr_polytomous[,2])
################################# 

#build final table with all ASR results 
ones = data.frame(ace = ace_dataframe$one , mpr_acc = mpr_acc_dataframe$one, mpr_del= mpr_del_dataframe$one, node = ace_dataframe$node)
# ones is where all the reconstruction data is stored.
# "ones" only has the internal nodes information






## --- B. Find transitions ---------------------------------------------------------------------------------

#loop for each ASR method
methodlabel <-  c("ML","ACC","DEL")  ### 1 for ML, 2 for ACCTRAN, 3 for DELTRAN  

for(k in 1:3){
  method = k # determines which method is used
  
  #create a dataframe of all states (tips and internal nodes)
  tips = mpr_tipstates
  node = c(1:length(tree$tip.label), ones$node)
  states = data.frame(ace = c(tips, ones$ace),acctran = c(tips, ones$mpr_acc), deltran = c(tips, ones$mpr_del), node = node)
  # "states" also has all the reconstruction data stored in it, its difference with "ones" is that it has tips info as well
  
  
  #sort edge table by the second column (children)
  edge = tree$edge
  edge = data.frame(parent = edge[,1], child = edge[,2])
  edge = edge[order(edge$child),]
  
  
  # finding *all* transition events in the tree
  o = vector()
  for (i in 1:dim(edge)[1]){
    s = states[edge$child[i], method]
    if(s == 1){
      o = append(o,i)
    }
  }
  transitions_children = edge[o,]
  z = vector()
  for (i in 1:dim(transitions_children)[1]){
    s = states[transitions_children$parent[i], method]
    if(s == 0){
      z = append(z,i)
    }
  }
  transitions_parents = transitions_children[z,]
  #trans_nodes includes ALL events where e non-QC node has a QC child. this includes embedded transitions and singletons
  trans_nodes = sort(unique(transitions_parents$child))
  
  
  
  #function for finding the parents of a node up to the root
  # inputs: tree, node number outputs: a vector of all parents of the node up to the root
  parent_finder = function(tree, nodes){ #function starts here
    parents_list = list()
    for (node in nodes){
      edge = tree$edge
      parents = vector()
      flag = 0
      root = length(tree$tip.label) + 1
      
      parents = append(parents, node)
      while(flag == 0){
        if (node == root){
          flag = 1
          next
        }
        ind = which(edge[,2] == node)
        node = edge[ind, 1]
        parents = append(parents, node)
        if (node == root){
          flag = 1
        }
      }
      parents_list = append(parents_list, list(parents))
    }
    return(parents_list)
  } #function ends here
  
  parents_list = parent_finder(tree, trans_nodes)
  
  #function for finding transition parents of transition nodes
  # this function gets a list of vectors representing nodes ancestors, returns transition events among them
  
  
  checked_list = vector()
  transitions_obj = list()
  for (i in 1:length(parents_list)){ #loop starts here
    p = unlist(parents_list[i])
    if (p[1] %in% checked_list){
      next
    }
    
    parents = vector()
    for (j in 1:length(p)){
      if (p[j] %in% trans_nodes){
        parents = append(parents, p[j])
      }
    }
    checked_list = append(checked_list, parents)
    transitions_obj = append(transitions_obj, list(parents))
    
  }#loop ends here
  
  independent_transitions = vector()
  for (i in 1:length(transitions_obj)){ #loop starts here
    p = unlist(transitions_obj[i])
    independent_transitions = append(independent_transitions, p[length(p)])
  }#loop ends here
  independent_transitions = unique(independent_transitions)
  
  new_independent_transitions = independent_transitions
  new_trans_nodes = trans_nodes
  
  
  
  # finding the non-qc parents of "new transitions"
  new_independent_transitions_parents = vector()
  for (i in 1:length(new_independent_transitions)){
    ind = which(edge[,2] == new_independent_transitions[i])
    new_independent_transitions_parents = append(new_independent_transitions_parents, edge[ind,1])
  }
  
  
  ######## finding singletons:
  # here, every QC tip that has a nonQC parent is counted as a singleton even if some of them have common parents
  # then check every singleton's parents up to the root, and removed the ones that are embedded in new_independent_transitions clades
  # exc_singletons are all singletons not embedded in the clades, ALL OF THESE ARE TIPS
  # new_independent_transitions are transition nodes, ALL OF THESE ARE NODES
  ## if the code malfunctioned, remove everything below this line and it should work fine, though singletons should be counted again
  
  edge = tree$edge
  one_tips = which(mpr_tipstates == 1)
  their_parents = vector()
  for (i in 1:length(one_tips)){
    ind = which(edge[,2] == one_tips[i])
    their_parents = append(their_parents, edge[ind,1])
  }
  #their_parents = unique(their_parents)
  singletons = vector()
  
  for (i in 1:length(their_parents)){
    j = which(ones$node == their_parents[i])
    if (ones[j, method] == 0){
      singletons = append(singletons, one_tips[i])
    }
  }
  singletons = unique(singletons)
  
  
  #checking if the singleton is embedded in a clade already marked
  root = length(tree$tip.label) + 1
  embedded = vector() #this is a vector of all singletons that are embedded in the clades
  for (i in 1:length(singletons)){
    flag = singletons[i] #loop should start here
    while (flag != root){
      ch = which(tree$edge[,2]== flag)
      flag = edge[ch,1]
      if (flag %in% new_independent_transitions){
        embedded = append(embedded, singletons[i])
        break
      }
    }
  }
  
  e = match(embedded, singletons)
  exc_singletons = singletons[-e]
  
  
  #store:
  df <- data.frame(transition_nodes = new_independent_transitions, parents = new_independent_transitions_parents)
  write.csv(df, paste0("transitions_", methodlabel[k], ".csv"), row.names = F, quote = F)
  assign(paste0("dftransitions_", methodlabel[k]), df, .GlobalEnv) 
}





##############################################################################################################
##################### 2. -- Build Intros Table  -- ###########################################################


#Requires: 
# ones, dftransitions_ML, dftransitions_ACC, dftransitions_DEL
#to be loaded in the global environment


### ------- Loop for building 3 tables -----------------------------------------------------------------------

for(h in 1:3){
  
  #Load recon results
  dfT <- get(paste0("dftransitions_", methodlabel[h]))
  new_independent_transitions <- dfT$transition_nodes 
  new_independent_transitions_parents <- dfT$parents
  
  #Assign intros
  intros2 <- unique(dfT$parents) 
  nn <- length(intros2)
  
  #Assign node character (recon results)
  #recon results: ones dataframe has all 3 methods' results: 0 = nonQC, 1 = QC
  nodesQC <- ones$node[ones[,h]==1] #load specific method's recon results
  
  #-------- Pre-processing for tree exploration -----------------
  
  #make list of QC labeled tips and nodes, combine
  tipsQC <- edge_table$child[str_detect(edge_table$chi.name, "Qc-")] #list of QC tips
  #combine
  vecQC <- c(tipsQC, nodesQC)
  #make list of nonQC labeled tips and nodes, combine
  allnonQC <- edge_table$child[!str_detect(edge_table$chi.name, "Qc-")]
  dfallnonQC <- filter(edge_table, edge_table$child %in% allnonQC )
  tipsnonQC <- dfallnonQC$child[dfallnonQC$isTip==TRUE] #list of nonQC tips
  nodesnonQC <-  dfallnonQC$child[dfallnonQC$isTip==FALSE] #list of nonQC inner nodes
  vecnonQC2 <- c(tipsnonQC, nodesnonQC)
  #nums from tree
  Ntip <- length(tree$tip.label) #number of tips
  Nnode <- tree$Nnode #number of nodes
  tipnums <- seq(1, Ntip) #list tips
  nodenums <- seq((Ntip+1), (Ntip+Nnode)) #list of nodes
  #---------------------------------------------------------
  
  
  #-------- Build empty dataframe to fill : 
  
  #empty df to fill:
  df_Intros <- data.frame(node = rep(NA, length(intros2)),            #introduction node (nonQC node)
                          basalQC_lab = rep(NA, length(intros2)),     #basal QC leaf label
                          outgroup_lab = rep(NA, length(intros2)),    #nearest outgroup leaf's label
                          outgroup_country = rep(NA, length(intros2)),#nearest outgroup leaf's origin (i.e. where it was sampled, e.g. Ontario)
                          origin_country = rep(NA, length(intros2)),  #inferred country of origin (final decision)
                          origin_region = rep(NA, length(intros2)),   #inferred region of origin (based on origin_country)
                          node_date =rep(NA, length(intros2)))        #date of the node based on time tree from Nextstrain
  
  
  
  ## ------ Fill in intro nodes with basal QC leaves (1st 2 columns)
  
  #find intros with basal QC leaves
  introsbasalQC <-c()
  for(i in 1:length(intros2)){
    intro <- intros2[i]
    tipkids <- edge_table %>% filter(parent == intro & isTip == TRUE)
    if(any(str_detect(tipkids$chi.name, "Qc-"))){introsbasalQC <-c(introsbasalQC, intro)}
  }
  
  # 1. Nodes column: store nodes into df_intros
  for (i in 1:length(introsbasalQC)) {  #only do this with those with QC basal leaf
    introsbasalQC[i] -> df_Intros[i, 1] #store intros in column 1
  }
  
  #find label of basal QC leaf
  introsbasalQClab <-c()
  for(i in 1:length(introsbasalQC)){
    intro <- introsbasalQC[i]
    tipkids <- edge_table %>% filter(parent == intro & isTip == TRUE)
    qckids <- tipkids$chi.name[str_detect(tipkids$chi.name, "Qc-")]
    sub <- filter(tipkids, tipkids$chi.name %in% qckids)
    ifelse(length(qckids)==1, introsbasalQClab <-c(introsbasalQClab, qckids),
           introsbasalQClab <-c(introsbasalQClab, NA)) }
  
  #find label of ones with >1 basal QC leaf (choose one with shortest branch, if tie choose one with travel history other wise coin toss)
  multibasal <-introsbasalQC[is.na(introsbasalQClab)]
  lab<-c()
  for (i in 1:length(multibasal)){
    tipkids <-edge_table %>% filter(parent == multibasal[i] & isTip == TRUE)
    qckids <- tipkids$chi.name[str_detect(tipkids$chi.name, "Qc-")]
    sub <- filter(tipkids, tipkids$chi.name %in% qckids)
    min <-sub$chi.name[sub$branch_length == min(sub$branch_length)]
    if(length(min)==1){lab<-c(lab, min)}
    if(length(min)==1){next}
    if(length(min)>1){lab <-c(lab, min[sample(length(min),1)])}
    }
  
  #fill in NAs in labels list
  indv <- which(is.na(introsbasalQClab))
  for(i in 1:length(lab)){ # lab and # of NAs should be the same length, so should multibasal
    introsbasalQClab[indv[i]] <- lab[i]
  }
   
  # 2. Store QC basal labels in df_intros
  for (i in 1:length(introsbasalQClab)) {  #only do this with those with QC basal leaf
    introsbasalQClab[i] -> df_Intros[i, 2] #store in column 2
  }
  
  
  
  ## ------ Fill in with intros that have an inner QC node before the basal QC leaf (column 1 and 2)
  #(move up reconstruction one node)
  
  restintros <-setdiff(intros2, introsbasalQC)
  
  #store nodes into df_intros
  for (i in 1:length(restintros)) {  #only do this with those with QC basal leaf
    restintros[i] -> df_Intros[(i+length(introsbasalQC)), 1] #store intros
  }
  
  #find label of basal QC leaf
  restintrosQClab <-c()
  flaggedintros <-c()
  for(i in 1: length(restintros)){
    intro <- restintros[i]
    nodekids <- intersect(edge_table$child[edge_table$parent == intro], nodenums) #find children that are nodes
    leng <- length(restintrosQClab)
    for(j in 1: length(nodekids)){
      sub <- edge_table %>% filter(parent == nodekids[j]) #this method favours the first nodekids that has a QC tip instead of comparing across all nodekids
      if(!any(str_detect(sub$chi.name, "Qc-"))) { next}  #if none are QC next nodekids
      if(any(str_detect(sub$chi.name, "Qc-"))) {#if at least one is QC
        sub2 <-sub %>% filter(str_detect(chi.name, "Qc-"))
        min2 <-sub2$chi.name[sub2$branch_length == min(sub2$branch_length)]
        if(length(min2)==1){restintrosQClab<-c(restintrosQClab, min2)}
        if(length(min2)==1){ break} #exit nested for loop if value added to restintrosQClab
        if(length(min2)>1){restintrosQClab <-c(restintrosQClab, min2[sample(length(min2),1)])}
        break
      }
    }
    if(length(restintrosQClab) == leng){
      restintrosQClab <- c(restintrosQClab, NA)
      flaggedintros <- c(flaggedintros, intro)
    }
  }
  length(flaggedintros) 
  
  #--store
  for (i in 1:length(restintrosQClab)) {  #
    restintrosQClab[i] -> df_Intros[(i+length(introsbasalQClab)), 2] #finish filling in column 2
  }
  
  
  ## ------ Fill in outgroup and it's country
  
  #label of outgroup
  outgrouplab <-c()
  for(i in 1: length(df_Intros$node)){
    node <-df_Intros$node[i]
    kids <-edge_table$child[edge_table$parent == node]
    sub1 <-edge_table %>% filter(parent %in% kids & isTip == TRUE & !str_detect(chi.name, "Qc-"))
    #check all children for nonQC tip with shortest branch length
    if(dim(sub1)[1]!=0){ #if sub1 has entries
      outg <-sub1$chi.name[sub1$branch_length == min(sub1$branch_length)]
      if(length(outg)>0){paste(outg, collapse=";")-> outg}
      outgrouplab <-c(outgrouplab, outg)
    }
    #if doesn't exist, go back one node and check again for nonQC tip
    if(dim(sub1)[1]==0){ #if sub1 has columns but no observations
      node2 <-edge_table$parent[edge_table$child == node]
      sub2 <-edge_table %>% filter(parent %in% node & isTip == TRUE & !str_detect(chi.name, "Qc-"))
      ifelse( dim(sub2)[1]==0,
              outgrouplab <-c(outgrouplab, "unclear"),
              outgrouplab<-c(outgrouplab, paste(sub2$chi.name[sub2$branch_length == min(sub2$branch_length)], collapse = ";")))
    }
  }
  
  df_Intros$outgroup_lab <-outgrouplab #fill in
  
  #assign country
  outgrouprss <-c()
  for(i in 1:length(df_Intros$outgroup_lab)){
    lab <- df_Intros$outgroup_lab[i]
    if(str_detect(lab, ";")){ outgrouprss <-c(outgrouprss, "TBD")}
    if(lab=="unclear"){ outgrouprss <-c(outgrouprss, "unclear")}
    if( (!str_detect(lab, ";")) & (lab!="unclear")){
      outgrouprss <-c(outgrouprss,metadata$country[metadata$strain == lab])}
  }
  
  df_Intros$outgroup_country <-outgrouprss #fill in
  
  #NOTE: visually check those "TBD" (to be determined) cases, since some might be inferrable from the tree
  
  
  ## ------ Fill in origins:
  
  #intro node dates (tmrca's of introductions)
  datesvec <-c()
  for(i in 1:length(df_Intros$node)){ #screws up when there's NAs
    datesvec <-c(datesvec, edge_table$date[edge_table$child == df_Intros$node[i]] )}
  
  df_Intros$node_date <-datesvec
  
  
  #designate the origin (and the region) of the introduction:
  
  #origin_country (without travel history this can only be assigned by the outgroup)
  df_Intros$origin_country <- df_Intros$outgroup_country #thus these columns are identical

  #origin_region
  df_Intros$origin_region <- rep(NA, length(df_Intros$node))
  for(i in 1:length(df_Intros$origin_region)){
    ifelse(any(match$Country.or.Area==df_Intros$origin_country[i]),
           match$Region.3[match$Country.or.Area==df_Intros$origin_country[i]]-> df_Intros$origin_region[i],
           next)
  }
 
  
  
  ## ------ Characterize QC transmission lineages
  
  # 1. size, i.e. number of QC sequences in the clade 
  linsize <- c()
  for(i in 1:length(df_Intros$node)){ #
    nodei <- df_Intros$node[i]
    cladei <- extract.clade(tree, node = nodei)
    #don't fill those with basal polytomy:
    ifelse(sum(edge_table$isTip[edge_table$parent== nodei])>2, #if has >2 tips, then it's a polytomy
           linsize <- c(linsize, NA),
           linsize <- c(linsize, length(cladei$tip.label[str_detect(cladei$tip.label, "Qc-")])))
  }
  #linsize <- c(linsize, NA)
  df_Intros$lin_size <- linsize #QC lineage size
  
  # 2. is node a polytomy? y/n
  polyt <-c()
  for(i in 1:length(df_Intros$node)){ #
    nodei <- df_Intros$node[i]
    cladei <- extract.clade(tree, node = nodei)
    ifelse(sum(edge_table$isTip[edge_table$parent== nodei])>2, #if has >2 tips, then it's a polytomy
           polyt <- c(polyt, "y"), polyt <- c(polyt, "n"))
  }
  #polyt <-c(polyt, NA)
  df_Intros$polytomy <-polyt
  
  
  # 2.b add lineage sizes of polytomies
  polyies <- df_Intros$node[df_Intros$polytomy == "y"] #list of intro node with a polytomy as a MRCA
  listlineages <- list() #nested list of QC seqs in each polytomy intro
  linsizepoly <- c() #assigned lin_size for each polytomy node based on  listlineages result
  
  for(i in 1:length(polyies)){ #
    embeddedqcseqs <- c() #clear for every node
    nodei <- polyies[i]
    cladei <- extract.clade(tree, node = nodei)
    nodesincladei <- edge_table$child[edge_table$chi.name %in% cladei$node.label]
    qcseqincladei <- cladei$tip.label[str_detect(cladei$tip.label, "Qc-")]
    embedded <- df_Intros$node[df_Intros$node %in% nodesincladei] #embedded intro nodes
    embedded <- embedded[embedded != nodei] #remove nodei from list
    if(length(embedded) > 0){
      for(j in 1:length(embedded)){
        cladej <- extract.clade(tree, node = embedded[j])
        add <-  cladej$tip.label[str_detect(cladej$tip.label, "Qc-")]
        embeddedqcseqs <- c(embeddedqcseqs, add)  }}
    if(!is.null(embeddedqcseqs)){qcseqincladei <- setdiff(qcseqincladei, embeddedqcseqs)}
    linsizepoly <- c(linsizepoly, length(qcseqincladei)) #store in lin size list
    listlineages[[as.character(nodei)]] <- qcseqincladei #store qc seqs in nested list
  }
  for(i in 1:length(polyies)){
    df_Intros$lin_size[df_Intros$node == polyies[i]] <- linsizepoly[i] #fill in NAs
  }
  
 
  
  ## ------ Fill in pango lineages and dates
  
  listlin <-c()
  for(i in 1:length(df_Intros$basalQC_lab)){
    ifelse(is.na(df_Intros$basalQC_lab[i]),
           listlin <- c(listlin, NA),
           listlin <- c(listlin, metadata$cladePang[metadata$strain == df_Intros$basalQC_lab[i]]))
  }
  df_Intros$pangolins <- listlin
  
  #add dominant pango-lineage of each intro clade to df_Intros
  listlin2 <-c()
  for(i in 1:length(df_Intros$node)){
    testclade <- extract.clade(tree, df_Intros$node[i])
    lin <-names(which.max(table(metadata$cladePang[metadata$strain %in% testclade$tip.label[str_detect(testclade$tip.label, "Qc-")]])))
    listlin2 <-c(listlin2, lin)
  }
  df_Intros$pangolinmax <- listlin2
  
  startdate<-as.Date(c())
  for(i in 1:length(df_Intros$node)){
    testclade <- extract.clade(tree, df_Intros$node[i])
    dat<-min(metadata$date[metadata$strain %in% testclade$tip.label[str_detect(testclade$tip.label, "Qc-")]])
    startdate <-c(startdate, dat)
  }
  df_Intros$start_date <-startdate
  
  enddate<-as.Date(c())
  for(i in 1:length(df_Intros$node)){
    testclade <- extract.clade(tree, df_Intros$node[i])
    dat<-max(metadata$date[metadata$strain %in% testclade$tip.label[str_detect(testclade$tip.label, "Qc-")]])
    enddate <-c(enddate, dat)
  }
  df_Intros$end_date <-enddate
  
  #store:
  assign(paste0("dfIntros", methodlabel[h]), df_Intros, .GlobalEnv)
  
}

# ----------- end of loop



## ---- Build one final results table ---------------------------------------------

dfIntrosML$ASRmethod <- rep("ML",length(dfIntrosML$node))
dfIntrosACC$ASRmethod <- rep("ACC",length(dfIntrosACC$node))
dfIntrosDEL$ASRmethod <- rep("DEL",length(dfIntrosDEL$node))
dfAll <- rbind(dfIntrosML, dfIntrosACC, dfIntrosDEL)

#harmonize 'unclear' category
dfAll$outgroup_country <- str_replace_all(dfAll$outgroup_country, "TBD", "unclear")
dfAll$origin_country <- str_replace_all(dfAll$origin_country, "TBD", "unclear")
dfAll$origin_region <- str_replace_na(dfAll$origin_region, "unclear")



write.csv(dfAll, "Table_S3_introductions.csv", quote = F)

toc()