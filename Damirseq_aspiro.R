# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")

setwd("Box/gandreoletti/ASPIRO trial_03.19.2021 (1)/")
setwd("~/Dropbox/Audentes_gaia/")
<!-- BiocManager::install("DaMiRseq") -->

## ----style-knitr, eval=TRUE, echo=FALSE, results="asis"--------------------
BiocStyle::latex(relative.path = TRUE)

## ----chu_1-----------------------------------------------------------------
library(DaMiRseq)
## only for example:
# rawdata.path <- system.file(package = "DaMiRseq","extdata")
# setwd(rawdata.path)
filecounts <- read.csv(file = "all_cohort_counts_all_54Samples.csv",row.names = 1,header = T, check.names=FALSE)
filecounts <- filecounts[ -c(1:10) ]

filecovariates <- read.table("all_54_samples_noControls.txt",row.names = 1, header = 1)
# filecovariates$sampleID <- rownames(filecovariates)
count_data <- round(filecounts)
#  covariate_data <- read.delim(filecovariates)

ncol(count_data) == nrow(filecovariates)
count_data <- count_data[,order(colnames(count_data))]
covariate_data <- filecovariates[order(rownames(filecovariates)),]
colnames(count_data) == rownames(covariate_data)

SE<-DaMiR.makeSE(count_data, covariate_data)

## ----chu_2-----------------------------------------------------------------
# data(SE)
assay(SE)[1:5, c(1:5, 21:25)]
colData(SE)

## ----chu_4-----------------------------------------------------------------
data_norm <- DaMiR.normalization(SE, minCounts=10, fSample=0.7,
                                 hyper = "no")


## ----chu_5-----------------------------------------------------------------
data_norm <- DaMiR.normalization(SE, minCounts=10, fSample=0.7,
                                 hyper = "yes", th.cv=3)
print(data_norm)
assay(data_norm)[c(1:5), c(1:5, 21:25)]

## ----chu_6, eval=FALSE-----------------------------------------------------
#  # Time Difference, using VST or rlog for normalization:
#  #
data_norm <- DaMiR.normalization(SE, minCounts=10, fSample=0.7, th.cv=3)
#  # VST: about 80 seconds
#  #
#  #data_norm <- DaMiR.normalization(SE, minCounts=10, fSample=0.7, th.cv=3,
#  #                                 type="rlog")
#  # rlog: about 8890 seconds (i.e. 2 hours and 28 minutes!)

## ----chu_7-----------------------------------------------------------------
data_filt <- DaMiR.sampleFilt(data_norm, th.corr=0.8)
dim(data_filt)

typeof (data_filt)
## ----chu_8, dev="pdf"------------------------------------------------------
sv <- DaMiR.SV(data_filt)

## ----chu_9, dev="pdf"------------------------------------------------------
DaMiR.corrplot(sv, colData(data_filt), sig.level = 0.01)

## ----chu_10, dev="pdf"-----------------------------------------------------
data_adjust<-DaMiR.SVadjust(data_filt, sv, n.sv=9)
assay(data_adjust[c(1:5), c(1:5, 21:25)])

## ----chu_11, dev="pdf"-----------------------------------------------------
# After gene filtering and normalization
DaMiR.Allplot(data_filt, colData(data_filt))

## ----chu_12, dev="pdf"-----------------------------------------------------
# After sample filtering and sv adjusting
DaMiR.Allplot(data_adjust, colData(data_adjust))

## ----chu_export, dev="pdf", eval=FALSE-------------------------------------
#  outputfile <- "DataNormalized.txt"
#  write.table(data_norm, file = outputfile_norm, quote = FALSE, sep = "\t")

## ----chu_13----------------------------------------------------------------
set.seed(12345)
data_clean<-DaMiR.transpose(assay(data_adjust))
df<-colData(data_adjust)
data_reduced <- DaMiR.FSelect(data_clean, df, th.corr=0.4)

## ----chu_14, dev="pdf"-----------------------------------------------------
data_reduced <- DaMiR.FReduct(data_reduced$data)
DaMiR.MDSplot(data_reduced, df)

## ----chu_15, dev="pdf"-----------------------------------------------------
# Rank genes by importance:
df.importance <- DaMiR.FSort(data_reduced, df)
head(df.importance)

## ----chu_16, dev="pdf"-----------------------------------------------------
# Select Best Predictors:
selected_features <- DaMiR.FBest(data_reduced, ranking=df.importance,
                                 n.pred = 5)
selected_features$predictors
# Dendrogram and heatmap:
DaMiR.Clustplot(selected_features$data, df)

## ----chu_17, dev="pdf"-----------------------------------------------------

Classification_res <- DaMiR.EnsembleLearning(selected_features$data,
                                             classes=df$class, fSample.tr = 0.5,
                                             fSample.tr.w = 0.5, iter = 30)

## ----chu_17bis1, dev="pdf"-------------------------------------------------
# Dataset for prediction
set.seed(10101)
nSampl_cl1 <- 5
nSampl_cl2 <- 5

## May create unbalanced Learning and Test sets
# idx_test <- sample(1:ncol(data_adjust), 10)

# Create balanced Learning and Test sets
idx_test_cl1<-sample(1:(ncol(data_adjust)/2), nSampl_cl1)
idx_test_cl2<-sample(1:(ncol(data_adjust)/2), nSampl_cl2) + ncol(data_adjust)/2
idx_test <- c(idx_test_cl1, idx_test_cl2)

Test_set <- data_adjust[, idx_test, drop=FALSE]
Learning_set <- data_adjust[, -idx_test, drop=FALSE]


## ----chu_17bis12, dev="pdf"------------------------------------------------

# Training and Test into a 'nfold' Cross Validation
nfold <- 3
cv_sample <- c(rep(seq_len(nfold), each=ncol(Learning_set)/(2*nfold)),
               rep(seq_len(nfold), each=ncol(Learning_set)/(2*nfold)))

# Variables initialization
cv_models <- list()
cv_predictors <- list()
res_df <- data.frame(matrix(nrow = nfold, ncol = 7))
colnames(res_df) <- c("Accuracy",
                      "N.predictors",
                      "MCC",
                      "sensitivity",
                      "Specificty",
                      "PPV",
                      "NPV")

## ----chu_17bis2, dev="pdf", results="hide"---------------------------------
for (cv_fold in seq_len(nfold)){
  
  # Create Training and Validation Sets
  idx_cv <- which(cv_sample != cv_fold)
  TR_set <- Learning_set[,idx_cv, drop=FALSE]
  Val_set <- Learning_set[,-idx_cv, drop=FALSE]
  
  #### Feature selection
  data_reduced <- DaMiR.FSelect(t(assay(TR_set)),
                                as.data.frame(colData(TR_set)),
                                th.corr=0.4)
  data_reduced <- DaMiR.FReduct(data_reduced$data,th.corr = 0.9)
  df_importance <- DaMiR.FSort(data_reduced,
                               as.data.frame(colData(TR_set)))
  selected_features <- DaMiR.FBest(data_reduced,
                                   ranking=df_importance,
                                   autoselect = "yes")
  # update datasets
  TR_set <- TR_set[selected_features$predictors,, drop=FALSE]
  Val_set <- Val_set[selected_features$predictors,drop=FALSE]
  
  ### Model building
  ensl_model <- DaMiR.EnsL_Train(TR_set,
                                 cl_type = c("RF", "LR"))
  # Store all trained models
  cv_models[[cv_fold]] <- ensl_model
  
  ### Model testing
  res_Val <- DaMiR.EnsL_Test(Val_set,
                             EnsL_model = ensl_model)
  
  # Store all ML results
  res_df[cv_fold,1] <- res_Val$accuracy[1] # Accuracy
  res_df[cv_fold,2] <- length(res_Val$predictors) # N. of predictors
  res_df[cv_fold,3] <- res_Val$MCC[1]
  res_df[cv_fold,4] <- res_Val$sensitivity[1]
  res_df[cv_fold,5] <- res_Val$Specificty[1]
  res_df[cv_fold,6] <- res_Val$PPV[1]
  res_df[cv_fold,7] <- res_Val$NPV[1]
  
  cv_predictors[[cv_fold]] <- res_Val$predictors
}


## ----chu_17bis3, dev="pdf"-------------------------------------------------
# Model Selection
res_df[,1:5]

idx_best_model <- DaMiR.ModelSelect(res_df,
                                    type.sel = "mode",
                                    npred.sel = "min")

## ----chu_17bis4, dev="pdf"-------------------------------------------------
# Prediction on the the independent test set
res_predict <- DaMiR.EnsL_Predict(Test_set,
                                  bestModel = cv_models[[idx_best_model]])

# Predictors
cv_predictors[[idx_best_model]]

# Prediction assessment for Ensemble learning
id_classifier <- 1 # Ensemble Learning
table(colData(Test_set)$class, res_predict[,id_classifier])

# Prediction assessment for Logistic regression
id_classifier <- 3 # Logistic regression
table(colData(Test_set)$class, res_predict[,id_classifier])


## ----chu_17ter1, dev="pdf"-------------------------------------------------
data(SE)
# create Independent test set and Learning set (raw counts)
idx_test <- c(18,19,39,40)
Ind_Test_set <- SE[, idx_test, drop=FALSE]
Learning_set <- SE[, -idx_test, drop=FALSE]

# DaMiRseq pipeline on Learning Set
data_norm <- DaMiR.normalization(Learning_set,
                                 minCounts=10,
                                 fSample=0.7,
                                 hyper = "yes",
                                 th.cv=3)
sv <- DaMiR.SV(data_norm)
data_adjust <- DaMiR.SVadjust(data_norm, sv, n.sv=4)


## ----chu_17ter2, dev="pdf"-------------------------------------------------
# remove not expressed genes from Learning_set and Ind_Test_set
expr_LearningSet <- Learning_set[rownames(data_norm)]
expr_Ind_Test_set <- Ind_Test_set[rownames(data_norm)]

# Independent test set Normalization
norm_ind_ts <- DaMiR.iTSnorm(expr_LearningSet,
                             expr_Ind_Test_set,
                             normtype = "vst",
                             method = "precise")

# Independent test set batch Adjusting
adj_norm_ind_ts <- DaMiR.iTSadjust(data_adjust, norm_ind_ts)


## ----chu_17ter3, dev="pdf"-------------------------------------------------
# Prediction on independent test set

prediction <- DaMiR.EnsL_Predict(t(adj_norm_ind_ts),
                                 bestModel = cv_models[[idx_best_model]])
prediction

# confusion matrix for the Ensemble Learner
table(Ind_Test_set@colData$class, prediction[,1])

## ----chu_18, dev="pdf"-----------------------------------------------------
## Feature Selection
set.seed(12345)
data_clean_2<-DaMiR.transpose(assay(data_filt))
df_2<-colData(data_filt)

data_reduced_2 <- DaMiR.FSelect(data_clean_2, df_2, th.corr=0.4)
data_reduced_2 <- DaMiR.FReduct(data_reduced_2$data)
df.importance_2 <- DaMiR.FSort(data_reduced_2, df_2)
head(df.importance_2)

selected_features_2 <- DaMiR.FBest(data_reduced_2, ranking=df.importance_2,
                                   n.pred=5)
selected_features_2$predictors

## Classification

Classification_res_2 <- DaMiR.EnsembleLearning(selected_features_2$data,
                                               classes=df_2$class,
                                               fSample.tr = 0.5,
                                               fSample.tr.w = 0.5,
                                               iter = 30)


## ----sessInfo, results="asis", echo=FALSE----------------------------------
toLatex(sessionInfo())
