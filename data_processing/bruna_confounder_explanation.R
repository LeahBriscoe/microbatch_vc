library(ggcorrplot)
library(mvtnorm)
library(glmnet)

# Define parameters
N=1000 # number of samples
n=1000 # number of OTUs  
b0=5 # effect of confounder batch and phenotype
b1=1 # effect of confounder batch on microbiome
b2=5  # effect of non-confounder batch on microbiome
b3=10 # effect of non-confounder batch on microbiome

# Simulate batch
batch1=rbinom(n = N,size = 1,prob = .5)
batch2=rbinom(n = N,size = 1,prob = .5)
batch3=rbinom(n = N,size = 1,prob = .5)

# Simulate phenotype with batch1 as confounder
phenotype=b0*batch1+rnorm(n = N)

# Correlation between batches and phenotype
# batch 1 is a potential confounder
ggcorrplot(corr = cor(cbind(batch1,batch2,batch3),phenotype), lab = T)

# Simulate microbiome with batch effect in the absence of phentype effect
microbiome=cbind(b1*batch1 + rmvnorm(n = N,mean = rep(0,round(n/3)),sigma = diag(round(n/3))),
                 b2*batch2 + rmvnorm(n = N,mean = rep(0,round(n/3)),sigma = diag(round(n/3))),
                 b3*batch3 + rmvnorm(n = N,mean = rep(0,round(n/3)),sigma = diag(round(n/3))))

# Principal component analysis
pcs=prcomp(microbiome,center = T, scale. = T)$x[,1:10]

# Correlation between batches and phenotype with microbiome PCs
# Induced correlation between microbiome and phenotype
ggcorrplot(corr = cor(pcs, cbind(batch1,phenotype,batch2,batch3)), lab = T)

# Predict phenotype from microbiome
microbiome_train=microbiome[1:round((2*N/3)),]
microbiome_test=microbiome[-c(1:round((2*N/3))),]

phenotype_train=phenotype[1:round((2*N/3))]
phenotype_test=phenotype[-c(1:round((2*N/3)))]

fit_lasso_cv = cv.glmnet(x = microbiome_train, y = phenotype_train, alpha = .5)
my_prediction_test=predict(fit_lasso_cv, newx = microbiome_test, type = "response", s = fit_lasso_cv$lambda.min)[,1]

cor(phenotype_test,my_prediction_test)^2 # erroneously high prediction performance

# Regress PC 1 and 2 (not capturing confounder batch) and predict phenotype from microbiome
microbiome_res=residuals(lm(microbiome~pcs[,1:2]))
microbiome_train=microbiome_res[1:round((2*N/3)),]
microbiome_test=microbiome_res[-c(1:round((2*N/3))),]

fit_lasso_cv = cv.glmnet(x = microbiome_train, y = phenotype_train, alpha = .5)
my_prediction_test_res=predict(fit_lasso_cv, newx = microbiome_test, type = "response", s = fit_lasso_cv$lambda.min)[,1]
cor(phenotype_test,my_prediction_test_res)^2 # erroneously (often) higher prediction performance


# Regress PC 3 which captures confounder batch and predict phenotype from microbiome
microbiome_res=residuals(lm(microbiome~pcs[,3]))
microbiome_train=microbiome_res[1:round((2*N/3)),]
microbiome_test=microbiome_res[-c(1:round((2*N/3))),]

fit_lasso_cv = cv.glmnet(x = microbiome_train, y = phenotype_train, alpha = .5)
my_prediction_test_res_no_confounder=predict(fit_lasso_cv, newx = microbiome_test, type = "response", s = fit_lasso_cv$lambda.min)[,1]

cor(phenotype_test,my_prediction_test_res_no_confounder)^2 # low or zero prediction performance

