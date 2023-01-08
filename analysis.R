#### Main codes for Paired fecal and blood metabolomics associate gut microbiota and cardiometabolic diseases (for questions, please contact Kui Deng, dengkui@westlake.edu.cn or dengkui_stat@163.com)
#1.Correlations between paired fecal and blood metabolites
#2.Paired comparison between fecal and blood metabolome in their associations with gut microbiota and microbial pathways
#3.Associations of well-predicted metabolites with cardiometabolic diseases

## All codes were executed using R, except for GREML analysis for genetic correlation, which was executed by GCTA tool

###1.Correlations between paired fecal and blood metabolites
##Phenotypic correlations 
corr_observation <- c()
count <- 1
for (var_temp in Metabolite_fecal){
  
  cat("##########",count,"###########\n")
  data_fit_temp <- data.frame(Met_fecal=data_all_metabolite_fecal_new[,var_temp],Met_serum=data_all_metabolite_serum_new[,var_temp],data_all_phenotype_new,stringsAsFactors = F)
  
  fit_temp_serum <- glm(Met_serum~age_follow+sex+BMI_follow,data=data_fit_temp)
  fit_temp_fecal <- glm(Met_fecal~age_follow+sex+BMI_follow,data=data_fit_temp)
  
  serum_residue <- fit_temp_serum$residuals
  fecal_residue <- fit_temp_fecal$residuals
  
  corr_temp <- cor.test(serum_residue,fecal_residue,method="spearman")
  corr_temp_new <- c(r=corr_temp$estimate,p=corr_temp$p.value)
  corr_temp_new <- c(Metabolite=var_temp,corr_temp_new)
  corr_observation <- rbind(corr_observation,corr_temp_new)
  count <- count+1
}

corr_observation <- as.data.frame(corr_observation)
corr_observation[,-1] <- apply(corr_observation[,-1],2,as.numeric)
corr_observation$FDR <- p.adjust(corr_observation$p,method="fdr")

##Genetic correlation (using bivariate GREML analysis by GCTA tool)
for i in $(seq 1 132);do

cat fecal_serum_merge/$i.txt | sed 's/\"//g' >temp 

gcta64 --reml-bivar --reml-bivar-nocove --grm pca_grm_ld_c --pheno temp --reml-bivar-lrt-rg 0 --out $i

done


###2.Paired comparison between fecal and blood metabolome in their associations with gut microbiota and microbial pathways
##Using random forest model with five-fold cross-validation
#Function
randomForest_MODEL<-function(X,
                             Y,
                             kfold=5,
                             seed=1234
){
  
  set.seed(seed)
  kfold10 <- sample(rep(1:kfold,ceiling(nrow(X)/kfold))[1:nrow(X)])
  
  prob <- matrix(NA,nrow(X),1)
  for (i in (1:kfold)){
    cat(paste("##########fold",i,"##############\n"))
    library(randomForest)
    set.seed(seed)
    rfvimSecond <- randomForest(X[kfold10!=i,],Y[kfold10!=i])
    
    prob[kfold10==i,1]  <- predict(rfvimSecond , X[kfold10==i,,drop=F], type="response")
  }
  
  cor_temp <- cor.test(Y,prob,method="spearman")
  
  r_value <- cor_temp$estimate
  p_value <- cor_temp$p.value
  
  result<-list(r_value=r_value,p_value=p_value,pred=prob)
  return(result)
}

#Example: for the association between species and metabolites
Model_fit <- randomForest_MODEL(X=data_all_species_new,Y=data_Metabolite_temp,kfold=5,seed=1234)
r_value_temp <- Model_fit$r_value
p_value_temp <- Model_fit$p_value


##Using LightGBM model with five-fold cross-validation
#Function
lightgbm_model_GM <- function(X,
                              Y,
                              kfold=5,
                              seed=1234,
                              categorical_feature=NULL){
  set.seed(seed)
  library(lightgbm)
  kfold10 <- sample(rep(1:kfold,ceiling(nrow(X)/kfold))[1:nrow(X)])
  
  prob <- matrix(NA,nrow(X),1)
  for (i in (1:kfold)){
    cat(paste("##########fold",i,"##############\n"))
    
    
    lgb.train = lgb.Dataset(data=as.matrix(X[kfold10!=i,]), label=Y[kfold10!=i],free_raw_data=F)
    
    lgb.grid = list(objective = "regression",
                    metric = "l2",
                    learning_rate=0.005,
                    feature_fraction=0.2,
                    min_data_in_leaf=15,
                    n_estimators=2000,
                    bagging_fraction=0.8,
                    bagging_freq=1
    )
    
    dtest <- lgb.Dataset(data=as.matrix(X[kfold10==i,]), label = Y[kfold10==i],free_raw_data=F)
    valids <- list(test = dtest)
    
    set.seed(seed)
    
    lgb.model = lgb.train(params = lgb.grid, data = lgb.train,record=F,valids=valids,categorical_feature=categorical_feature,verbose=-1)
    
    prob[kfold10==i,1]  <- predict(lgb.model,data=as.matrix(X[kfold10==i,]))
    
  }
  
  cor_temp <- cor.test(Y,prob,method="spearman")
  
  r_value <- cor_temp$estimate
  p_value <- cor_temp$p.value
  
  
  result<-list(r_value=r_value,p_value=p_value,pred=prob)
  return(result)
}

#Example: for the association between species and metabolites
Model_fit <- lightgbm_model_GM(X=data_all_species_new,Y=data_Metabolite_temp,kfold=5,seed=1234)
r_value_temp <- Model_fit$r_value
p_value_temp <- Model_fit$p_value

##Differences between the associations of gut microbiota/pathways with paired fecal and blood metabolites 
library(cocor)
Corr_diff_overlap_metabolites <- c()
for (i in 1:nrow(data_temp_A)){
  cat("#########",i,"###########\n")
  
  corr_temp <- cocor.dep.groups.overlap(r.jk=data_temp_A$r[i],r.jh=data_temp_B$r[i],
                                        r.kh=data_temp_C$r.rho[i],n=1007,test="hittner2003")
  
  
  result_temp <- data.frame(Metabolite=data_temp_A$sample[i],r_fecal=data_temp_A$r[i],
                            r_serum=data_temp_B$r[i],heterogeneity=corr_temp@hittner2003$p.value)
  Corr_diff_overlap_metabolites <- rbind(Corr_diff_overlap_metabolites,result_temp)
}


Corr_diff_overlap_metabolites$heterogeneity_FDR <- p.adjust(Corr_diff_overlap_metabolites$heterogeneity,method="fdr")


##Validation in the independent validation cohort
set.seed(1234)
library(randomForest)
RF_model <- randomForest(data_all_species_new,data_Metabolite_GNHS_temp)

prob  <- predict(RF_model , data_FH_species_new, type="response")


cor_temp <- cor.test(data_Metabolite_FH_temp,prob,method="spearman")

r_value_temp <- cor_temp$estimate
p_value_temp <- cor_temp$p.value


###3.Associations of well-predicted metabolites with cardiometabolic diseases
#using the multivariable logistic model
fit_temp <- glm(factor(phenotype)~Metabolite+age_follow+sex+BMI_follow+smoke_follow+alc_follow+factor(edu3)+factor(income4)+MET+energy,data=data_fit_temp,family=binomial(link="logit"))

