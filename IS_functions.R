

# Some functions that can be used in the context of diet analyses focusing on individual specialisation.
#
# Literature related to the functions are:
#
# Araújo, M. S., D. I. Bolnick, G. Machado, A. A. Giaretta, and S. F. dos Reis. 2007. Using ??13C stable isotopes to quantify individual-level diet variation. Oecologia 152:643-654.
# Bolnick, D. I., R. Svanbäck, M. S. Araújo, and L. Persson. 2007. Comparative support for the niche variation hypothesis that more generalized populations also are more heterogeneous. Proceedings of the National Academy of Sciences 104:10075-10079.
# Bolnick, D. I., L. H. Yang, J. A. Fordyce, J. M. Davis, and R. Svanbäck. 2002. Measuring individual-level resource specialization. Ecology 83:2936-2941.
# Lawlor, L. R. 1980. Overlap, Similarity, and Competition Coefficients. Ecology 61:245-251.
# Pianka, E. R. 1973. The Structure of Lizard Communities. Annual Review of Ecology and Systematics 4:53-74.
# Vander Zanden, M., and J. Rasmussen. 2001. Variation in d15N and d13C trophic fractionation: Implications for aquatic food web studies. Limnology and Oceanography 46:2061-2066.
# Winemiller, K. O., and E. R. Pianka. 1990. Organization in Natural Assemblages of Desert Lizards and Tropical Fishes. Ecological Monographs 60:27-55.


# The cum_prey function provide cummulative distributions of prey items when randomly sampling
# individuals within a population
# 
# arguments used are:
# data = diet matrix
# N = number of replicated cummulative distribution
#
# It returns a data.frame with N columns representing the numbe rof prey items 
# following successive random sampling(without replacement) from the observed diet matrix

cum_prey<-function(data,N){
  
  res<-data.frame(matrix(NA,ncol=N,nrow=nrow(data)))
  for(k in 1:N){
    
      # random sample stomach content
  n<-sample(seq(from=1,to=nrow(data),by=1))
  
  # intermediate dataset receiving prey names
  inter<-data.frame(matrix(NA,ncol=ncol(data),nrow=nrow(data)))
  for (i in 1:nrow(data)){
    prey<-which(data[n[i],]>0)
    name_preys<-names(data[n[i],])[prey]
    L<-length(name_preys)
    inter[i,1:L]<-name_preys
  }
  
  # Receive the cumulative prey items
  np<-data.frame(matrix(NA,ncol=ncol(data),nrow=nrow(data)))
  
  # Initialize the first prey item in np
  L2<-length(inter[1,-which(is.na(data.frame(inter[1,])))])
  np[1,1:L2]<-as.character(inter[1,-which(is.na(data.frame(inter[1,])))])
  
  # Loop to add progressively new items
  for (j in 1:(nrow(data)-1)){
    
    # Match between prey items in successive stomachs
   del<- which(is.na(data.frame(np[j,])))
   
  if(length(del)>0){L3<-match(inter[j+1,-which(is.na(data.frame(inter[j+1,])))],np[j,-del])
  
   # Conditions to fill "np" considering prey item for successive stomachs
    
    # different preys
    if(length(na.omit(L3))==0){
      L4<-length(data.frame(np[j,-which(is.na(data.frame(np[j,])))],inter[j+1,-which(is.na(data.frame(inter[j+1,])))]))
      q<-as.matrix(c(np[j,-which(is.na(data.frame(np[j,])))],inter[j+1,-which(is.na(data.frame(inter[j+1,])))]))
      np[j+1,1:L4]<-q
    }
    
    # different and identical preys 
    if(length(L3)>length(na.omit(L3))){
      L4<-length(data.frame(np[j,-which(is.na(data.frame(np[j,])))],inter[j+1,-which(is.na(data.frame(inter[j+1,])))][which(is.na(L3))]))
      q<-as.matrix(c(np[j,-which(is.na(data.frame(np[j,])))],inter[j+1,-which(is.na(data.frame(inter[j+1,])))][which(is.na(L3))]))
      np[j+1,1:L4]<-q
    }
    
    # identical preys  
    if(length(L3)==length(na.omit(L3))){
      np[j+1,]<-np[j,]
    }
  
  }
     
  if(length(del)==0){np[j+1,]<-np[j,]}
}

  
res[,k]=data.frame(apply(np,1,function (x) length(na.omit(x))))
  
    }

  return(res=res)
}





# The IS function compute the following different specialization indices
#
# The arguments used are :
# data = diet matrix, 
# method = "sum" or "average" to compute population niche breath qj
# 
# It returns:
# 1) indexes of diet diversity for each consumer: :
#   - PS  (Bolnick 2002) overlap between one indiv and its pop
#   - lambda (and p-values) (Bolnick 2002) How differ an individual from its pop (chi2)
#   - W (Bolnick 2007) # Corrected values

# 2) indexes of diet specilisation for the population :
#   - WIC, TNW, WIC/TNW  (Bolnick 2002)
#   - D (Levins, 1968)
#   - mean population PS = IS (Bolnick 2002)
#   - mean population V (Bolnick 2007)
#   - mean population W (Bolnick 2002)


IS<-function(data,method){
  
  pij<-data/rowSums(data)#proportion des donnees pour chaque consommateur
  pip<-rowSums(data)/sum(data)#toutes les proies pour un conso dans toute la pop bolnich 2002 pi point = pop 
  
  #qj = propo de ressource pour toute la pop = nicghe breadth
  if(method=="sum"){qj<-colSums(data)/sum(data)}#proportion moyenne
  if(method=="average"){qj<-colMeans(pij)}# mieu =proportion moyenne de proie pour toute la pop evite les effets de ce qui ont beaucoup de proie = biais 
 
  #alphaij<- sweep(data, 2, colSums(data), `/`)#= alpha pour tous les conso on sait quelle prop ol mange par rapport à la proportion mangée par la pop
 
  df<-apply(data,1,function(x) length(which(x>0)))
  
  inter<-matrix(NA,ncol=5,nrow=nrow(data))
  for (z in 1:nrow(data)){
    inter[z,1]<- pip[z]*(-sum(pij[z,]*log(pij[z,]),na.rm=T))
    inter[z,2]<-1-0.5*sum(abs(pij[z,]-qj))#ps
    inter[z,3]<-prod((qj/pij[z,])^data[z,]) #lamda
    inter[z,4]<-(prod((qj/pij[z,])^data[z,]))^(1/sum(data[z,]))#W = corrigé par raaport au nombre d'item alimentaire
    inter[z,5]<-pchisq(-2*log(prod((qj/pij[z,])^data[z,])), df=df[z], lower.tail = FALSE) #pvalue = landa suit qui2 donc on a prob ave df = au nmbre de proies diff de celui modelisé au hasard avec model null
  }
  
  #inter2<-matrix(NA,ncol=1,nrow=ncol(data))
  #for (r in 1:ncol(data)){
  #  inter2[r,1]<- qj[r]*(-sum(alphaij[,r]*log(alphaij[,r]),na.rm=T))
  #}
  
  #BICs=sum(pip*log(pip))-sum(inter2)
  
  res<-data.frame(WICs=sum(inter[,1]),TNWs=-sum(qj*log(qj)),
                  WICs_TNWs=sum(inter[,1])/-sum(qj*log(qj)),
                  D_Levins=1/sum(qj^2),IS=sum(inter[,2])/nrow(data),
                  V=1-(sum(inter[,2])/nrow(data)),W=mean(inter[,4]))
  
  return(list(res,data.frame(PS=inter[,2],lambda=inter[,3],p_lambda=inter[,5],W=inter[,4])))
}


# The make_null function simulates random diet for consumers based on multinomial sampling 
# according to Bolnick (2002). The probability that an item falls in a prey category is equal
# to the proportions of the prey category in the population diet (qj).
#
# The arguments used are : 
# data = diet matrix
# n = number of random diets to be simulated
# N = number of items to be distributed among prey categories (random draws)
# method = "sum" or "average" to obtain population diet breath (qj)
#
# It returns a data.frame corresponding to a simulated population of n consumers

make_null<-function(data,n,N,method){
  
  pij<-data/rowSums(data);pij<-na.omit(pij)# proportion des proies
  
  if(method=="sum"){qj<-colSums(data)/sum(data)} # donne beaucoup de poids à des indiv rare avec bcp d'item alim spé
  if(method=="average"){qj<-colMeans(pij)}#lisse les indiv rare 
  
  null<-t(rmultinom(n, N,qj))
  return(data.frame(null))
}


# the random_IS function computes random estimates of individual specialization (V, WIC_TNW and W)
# for a population based on the make_null function 
#
# The arguments used are : 
# data = diet matrix
# replicate = the number of random populations to be simulated 
# method = "sum" or "average" to obtain population diet breath (qj)
#
# It returns:
# 1) the observed individual specialization estimates
# 2) the probability that the observed individual specialization estimates differed from the random expectation
# 3) a data.frame with the individual specialization estimates for the random populations

random_IS<-function(data,replicate,method){
  
  ind_sp<-IS(data,method)
  
  obs_V<-ind_sp[[1]]$V
  obs_WICs_TNWs<-ind_sp[[1]]$WICs_TNWs
  obs_W<-ind_sp[[1]]$W
  obs<-data.frame(var=c("V","WICs_TNWs","W"),obs=rbind(obs_V,obs_WICs_TNWs,obs_W))
  
  items<-rowSums(data)
  
  res<-data.frame(matrix(NA,ncol=3,nrow=replicate));colnames(res)<-c("V","WICs_TNWs","W")
  for (i in 1:replicate){
    
    surr<-data.frame(matrix(NA,ncol=ncol(data),nrow=nrow(data)))
    for (j in 1:nrow(data)){
      surr[j,]<-make_null(data=data,n=1,N=items[j],method)# n=1 on prend pour chaque indiv une simulation de bouffe
    }
    
    if(length(which(colSums(surr)==0))>0){surr<-surr[,-which(colSums(surr)==0)]}# vire col = 0 car proba très faible 
    
    sim_ind_sp<-IS(surr,method)
    
    res[i,1]<- sim_ind_sp[[1]]$V
    res[i,2]<-sim_ind_sp[[1]]$WICs_TNWs
    res[i,3]<-sim_ind_sp[[1]]$W
    
  }
  
  p_V<-length(which(res[,1]>obs_V))/replicate# lequels simulés sont sup au v obs fait avec le mod nulle = proportion = sgnificatif de spe indiv
  p_WICs_TNWs<-length(which(res[,2]<obs_WICs_TNWs))/replicate  # spé = proche de 0 donc on cherche les inferieurs
  p_W<-length(which(res[,3]<obs_W))/replicate
  
  return(list(obs=obs, data.frame(p_V=p_V,p_WICs_TNWs=p_WICs_TNWs,p_W=p_W),data.frame(res)))
}




# The pairs_overlap function compute proportional similarity (Bolnick 2002) 
#
# The argument used is a diet matrix
#
# It returns:
# 1) mean population index
# 2) mean individual index
# 3) values for each pair of consumers
# 4) a matrix of all consumer pairs combinaisons

pairs_overlap<-function(data){
  
pij<-data/rowSums(data)

m=seq(from=1,to=nrow(pij),by=1)
uniq.pairs <- unique(as.data.frame(t(apply(expand.grid(m, m), 1, sort))))#test que les paire uniques 
uniq.pairs<-uniq.pairs[-which(uniq.pairs[,1]==uniq.pairs[,2]),]

overlap<-data.frame(matrix(NA,ncol=1,nrow=nrow(uniq.pairs)));names(overlap)<-"PS" #pour les couples
for (k in 1:nrow(overlap)){
  
  couple<-as.matrix(uniq.pairs[k,]) #quelle paire est unique ? enlève les indiv 1 avec indiv1
  overlap[k,1]<-sum(pmin(pij[couple[1],],pij[couple[2],])) #pmin prend la plus petite valeur 
  
  }
overlap<-data.frame(uniq.pairs,overlap)

  attach(overlap)

mean_ind_overlap<-matrix(NA,ncol=1,nrow=nrow(data));names(mean_ind_overlap)<-c("PS")# moyenne de l'overlap de indiv & à nombre de data lequel dans la col1 et col2 lequel est = de 1 à la dernière valeur overlap moyen pour un inviv pa rapp aux autres
for (l in 1:nrow(data)){
  mean_ind_overlap[l,1]<-mean(overlap[which(overlap[,1]==l | overlap[,2]==l),3])
 
}

mean_overlap<-mean(overlap[,3])#overlap moyen de toutes les paires uniques dans la pop

sim_mat_p<-data.frame(matrix(NA,ncol=nrow(data),nrow=nrow(data)))# mat d'overlap si on a le 1-2 on a pas le 2-1
for (i in 1:nrow(data)){
  for (j in 1:nrow(data)){
  
  sim_mat_p[i,j]<-sum(pmin(pij[i,],pij[j,]))
 
   }
}

return(list(mean_overlap=mean_overlap,mean_ind_overlap=mean_ind_overlap,
            overlap=overlap,sim_mat_p=sim_mat_p))
}

# The pianka_overlap function computes Pianka index (1973) for different resource utilisation metrics 
# (Lawlor 1980, Winemiller and Pianka 1990).
#
# The argument used are:
# data = diet matrix
# resource metric must be either:
#           - "proportion" = diet matrix converted to proportional resource utilization( pij)
#           - "electivity" = diet matrix converted to electivity resource utilisation (eij) 
#           - "G" = diet matrix converted to gij = sqrt(pij*eij) 
#
# It returns:
# 1) mean population index
# 2) mean individual index
# 3) Index for each pair of consumers
# 4) a matrix of diet overlap with all pair combinaisons

pianka_overlap<-function(data, resource_metric){
  
  if(resource_metric=="proportion"){metric<-data/rowSums(data)}
  if(resource_metric=="electivity"){metric<-sweep(data/rowSums(data), 2, colSums(data), `/`)}
  if(resource_metric=="G"){metric<-sqrt(data/rowSums(data)*sweep(data/rowSums(data), 2, colSums(data), `/`))}
  
  m=seq(from=1,to=nrow(metric),by=1)
  uniq.pairs <- unique(as.data.frame(t(apply(expand.grid(m, m), 1, sort))))
  uniq.pairs<-uniq.pairs[-which(uniq.pairs[,1]==uniq.pairs[,2]),]
  
  overlap_pairs<-data.frame(matrix(NA,ncol=1,nrow=nrow(uniq.pairs)))
  for (k in 1:nrow(overlap)){
    
    couple<-as.matrix(uniq.pairs[k,])
    
    if(sum(metric[couple[1],]*metric[couple[2],])==0){overlap_pairs[k,1]<-0}
    else{
overlap_pairs[k,1]<-sum(metric[couple[1],]*metric[couple[2],])/sqrt(sum(metric[couple[1],]^2)*sum(metric[couple[2],]^2))
        }
  }
  overlap_pairs<-data.frame(uniq.pairs,overlap_pairs);colnames(overlap_pairs)<-c("C1","C2","overlap")
  
  attach(overlap_pairs)
  
  mean_ind_overlap<-matrix(NA,ncol=1,nrow=nrow(data))
  for (l in 1:nrow(data)){
    
    mean_ind_overlap[l,1]<-mean(overlap_pairs[which(overlap_pairs[,1]==l | overlap_pairs[,2]==l),3])
  }
  
  mean_overlap<-mean(mean_ind_overlap)
  
  overlap_mat<-data.frame(matrix(NA,ncol=nrow(data),nrow=nrow(data)))
  for (i in 1:nrow(data)){
    for (j in 1:nrow(data)){
      
      if(sum(metric[i,]*metric[j,])==0){overlap_mat[i,j]<-0}
      else{
      overlap_mat[i,j]<-sum(metric[i,]*metric[j,])/sqrt(sum(metric[i,]^2)*sum(metric[j,]^2))
      }
    }
  }
  
  return(list(mean_overlap=mean_overlap,mean_ind_overlap=mean_ind_overlap,
              overlap_pairs=overlap_pairs,overlap_mat=overlap_mat))
}



# The function iso_convert is derived from Araùjo (2007) and converts the diet of a consumer
# to its expected stable isotope value considering individual variations for trophic fractionation 
# according to Vander Zanden and Rasmussen (2001) set at 1+-1 for d13C and 3.4+-1 for d15N
# 
# The arguments used are:
# data = diet matrix 
# iso = the stable isotopes of interest either "d13C" or "d15N"
# iso_values = the stable isotope values of each preys
# mass = the mass of each prey item
#
# It returns a vector of inferred stable isotope values for each consumer of the diet matrix

iso_convert<-function(data,mass,iso,iso_mean,iso_sd){
  
  pij<-data/sum(data)
  
  if(is.na(iso_sd)[1]){
    
    if(iso=="d13C"){  iso_cons=sum(((pij*mass)/sum(pij*mass))*iso_mean)+rnorm(1,1,0)}
    if(iso=="d15N"){  iso_cons=sum(((pij*mass)/sum(pij*mass))*iso_mean)+rnorm(1,3.4,0)}
    
  }else{
    
    iso_values<-matrix(NA,ncol=2,nrow=length(data))
    for (i in 1:length(data)){
      iso_values[i,] <-c(rnorm(1,iso_mean[i],iso_sd[i]),rnorm(1,iso_mean[i],iso_sd[i]))
    }
    
    if(iso=="d13C"){  iso_cons=sum(((pij*mass)/sum(pij*mass))*iso_values)+rnorm(1,1,1)}
    if(iso=="d15N"){  iso_cons=sum(((pij*mass)/sum(pij*mass))*iso_values)+rnorm(1,3.4,1)}
  }
  
  return(data.frame(iso_cons=iso_cons))
}


# The random_iso function uses the make_null function to simulate random populations 
# and then uses the convert_iso function to estimate isotope values for each consumer of the simulated populations.
# It finally estimates the isotope variances of the simulated populations. 
# The observed isotope variance is confronted to those obtained from simulated random populations to obtain a
# probability than the observed variance differed from random expectation.
#
# The arguments used are:
# data= diet matrix
# replicate = the number of random populations to be simulated
# method = "sum" or "average" to obtain population diet breath (qj)
# iso = the stable isotopes of interest either "d13C" or "d15N"
# iso_values = the averaged stable isotope values of each prey
# iso_sd = the sd of stable isotope values of each prey
# mass = the mass of each single prey item
# var_iso_obs = the observed isotope variance
#
# it returns:
# p_val = the probability that observed isotope variance differed from random expectation
# sim_var_iso = the isotope variance for the simulated random populations
# sim_iso = the inferred isotope values for each consumers within the simulated random populations


random_iso<-function(data,replicate,method,mass,iso,iso_values,var_iso_obs){
 
  sim_iso<-data.frame(matrix(NA,ncol=2,nrow=replicate*nrow(data)));colnames(sim_iso)<-c("replicate","sim_iso")
  sim_var_iso<-data.frame(matrix(NA,ncol=1,nrow=replicate));colnames(sim_var_iso)<-c("var_iso")
  
  items<-rowSums(data)
  
   for (i in 1:replicate){
   
    surr_iso<-data.frame(matrix(NA,ncol=1,nrow=nrow(data)))
    for (j in 1:nrow(data)){
      
      surr_iso[j,]<-iso_convert(data=make_null(data=data,n=1,N=items[j],method),mass=mass,iso=iso,iso_values=iso_values)
    }
 
    sim_var_iso[i,1]<- var(surr_iso[,1])
  
    sim_iso[(1+(i-1)*nrow(data)):(nrow(data)+(i-1)*nrow(data)),1]<-rep(i,nrow(data))
    sim_iso[(1+(i-1)*nrow(data)):(nrow(data)+(i-1)*nrow(data)),2]<-surr_iso
  
   }
  
  p_val<-length(which(sim_var_iso[,1]>var_iso_obs))/replicate
  
  return(list(p_val=p_val,sim_var_iso=sim_var_iso, sim_iso=sim_iso))
}



#  The make_null_R1 function is a randomisation procedure that randomly  reshuffle rows for each column
make_null_R1<-function(data){
  
  rand<-matrix(rep(sample(seq(from=1,to=nrow(data),by=1),replace=F),ncol(data)),ncol=ncol(data),byrow=T)
  
  random_mat<-matrix(NA,ncol=ncol(data),nrow=nrow(data))
  for (i in 1:ncol(data)){
    random_mat[,i]<-data[rand[,i],i]
  }
  
  return(random_mat=random_mat)
}









