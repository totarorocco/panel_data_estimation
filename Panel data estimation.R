# Project by Luca Colaci (3176608), Rocco Totaro (3159435) & Giovanni Parlangeli (3152699)

#### Setup ####

# Installing and loading the libraries used
libraries = c('plm', 'lmtest', 'haven', 'gplots', 'ggplot2', 'jtools', 'rstudioapi', 'broom', 'broom.mixed')
for(i in libraries){
  print(paste('Installing', i))
  if(system.file(package = i) == ""){
    install.packages(i)
  }
  library(i, character.only = TRUE)
}
rm(i)

#Loading the dataset
setwd(dirname(getActiveDocumentContext()$path))
data = read_dta("prodfn_data_4_5.dta")


#### Exploratory analysis ####
pdata = pdata.frame(data, index = c('ivar', 'tvar'))

# Checking for NA values
na_count = sum(is.na(pdata))
if (na_count > 0) {
  print("There are NA values in the dataset.")
  } else {
  print("There are no NA values in the dataset.")
}

# Checking for firms with different numbers of timesteps
firm_counts = table(pdata$ivar)
distinct_timesteps = unique(table(pdata$tvar))

if (length(distinct_timesteps) > 1) {
  print("Different firms have a different number of timesteps.")
  } else {
  print("All firms have the same number of timesteps.")
}

# Checking for timesteps with different numbers of firms
timestep_counts = table(pdata$tvar)
distinct_firms = unique(table(pdata$ivar))

if (length(distinct_firms) > 1) {
  print("Different timesteps have a different number of firms.")
  } else {
  print("All timesteps have the same number of firms.")
}

rm(distinct_firms, distinct_timesteps, firm_counts, na_count, timestep_counts)

# Plotting heterogeneity
plotmeans(y ~ ivar, 
          data = pdata, 
          main = 'Heterogeneity across firms',
          xlab = 'Firms',
          ylab = 'Average y',
          col = '#3f7e44',
          bars = FALSE,
          connect = FALSE,
          )


plotmeans(y ~ tvar, 
          data = pdata, 
          main = 'Heterogeneity across timesteps',
          xlab = 'Timesteps',
          ylab = 'Average y',
          col = '#3f7e44',
          barwidth = 1.5,
          barcol = '#192954',
          ccol = '#d2df64',
)
  


#### Fixed Effect model - Individual effect ####

fe_formula = y ~ k + l
fe = plm(fe_formula, data = pdata, model = 'within', effect = 'individual')
summary(fe)

# BP test for homoscedasticity
# H0: Homoscedasticity is present in data
# H1: Heteroscedasticity is present in data
bptest(fe)
# The test fails to reject homoscedasticity

# Woolridge test for serial correlation
# H0: Serial correlation is not present in the error term
# H1: Serial correlation is present in the error term
pwartest(fe)
# The test fails to reject no serial correlation
# In FE models we use a more specific test, different from the one we applied elsewhere


#### Fixed Effect model - Time effect ####

fe2 = plm(fe_formula, data = pdata, model = 'within', effect = 'time')
summary(fe2)

# BP test for homoscedasticity
bptest(fe2)
# The test fails to reject homoscedasticity

# Woolridge test for serial correlation
pwartest(fe2)
# The test reject no serial correlation


#### Random Effect model - Individual effect ####

re_formula = y ~ k + l
re = plm(fe_formula, data = pdata, model = 'random', effect = 'individual')
summary(re)

# BP test for homoscedasticity
bptest(re)
# The test fails to reject homoscedasticity

# Woolridge test for serial correlation
pbgtest(re)
# The test reject no serial correlation, so we will need to adjust var-cov matrices

vcov_re = vcovHC(re, type = 'HC1')


#### Random Effect model - Time effect ####

re2 = plm(re_formula, data = pdata, model = 'random', effect = 'time')
summary(re2)

# BP test for homoscedasticity
bptest(re2)
# The test fails to reject homoscedasticity

# Woolridge test for serial correlation
pbgtest(re2)
# The test reject fails to reject no serial correlation


#### Pooled model ####

pols_formula = y ~ l + k
pols <- plm(pols_formula, data = data, model = 'pooling')
summary(pols)

# BP test for homoscedasticity
bptest(pols)
# The test fails to reject homoscedasticity

# Woolridge test for serial correlation
pbgtest(pols)
# The test reject no serial correlation

vcov_pols = vcovHC(pols, type = 'HC1')


#### Between model ####

be_formula = y ~ l + k
be = plm(be_formula, data = pdata, model = 'between')
summary(be)

# BP test for homoscedasticity
bptest(be)
# The test fails to reject homoscedasticity

# Woolridge test for serial correlation
pbgtest(be)
# The test fails to reject no serial correlation


#### Twoway Fixed Effect model ####

fe3 = plm(fe_formula, data = pdata, model = 'within', effect = 'twoway')
summary(fe3)

# BP test for homoscedasticity
bptest(fe3)
# The test fails to reject homoscedasticity

# Woolridge test for serial correlation
pbgtest(fe3)
# The test reject no serial correlation


#### Twoway Random Effect model ####

re3 = plm(re_formula, data = pdata, model = 'random', effect = 'twoway')
summary(re3)

# BP test for homoscedasticity
bptest(re3)
# The test fails to reject homoscedasticity

# Woolridge test for serial correlation
pbgtest(re3)
# The test fails to reject no serial correlation


#### Specificity tests ####

# Testing POLS vs FE models
pFtest(fe, pols)
pFtest(fe2, pols)
pFtest(fe3, pols)

# The test reject no significant effects only in the time and two-way models,
# therefore discarding the individual FE model in favor of POLS

# Testing RE vs FE models
# Set 1: Individual RE vs FE (we use a robust variance-covariance matrix)
phtest(fe, re, cov.fe = vcov_re)
phtest(fe2, re, cov.fe = vcov_re)
phtest(fe3, re, cov.fe = vcov_re)

# In all cases the RE model is preferred over the FE

# Set 2: Time RE vs FE
phtest(fe, re2)
phtest(fe2, re2)
phtest(fe3, re2)

# In all cases the RE model is preferred over the FE

# Set 3: Two-way RE vs FE
phtest(fe, re3)
phtest(fe2, re3)
phtest(fe3, re3)

# In all cases the RE model is preferred over the FE

# Testing RE vs Between model
phtest(be, re)
phtest(be, re2)
# The Hausman test is not applicable in the Between - two-way RE case

# All tests fail to reject the null, therefore pointing towards the RE model

plot_summs(fe, re2, be, pols, plot.distributions = TRUE, model.names = c('Individual FE', 'Time RE', 'Between', 'POLS'), robust = 'HC1')

#### Conclusion & Final visualizations ####
# Chosen model: Time Random Effect model

summary(re2)

# Plotting estimated coefficients and their distributions
plot_summs(re2, plot.distributions = TRUE, model.names = c('Time RE'), robust = 'HC1')



