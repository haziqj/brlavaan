##==============================================================================
## Project:     Resampling based bias correction for small sample SEM
##
## Script:      Illustration on how to obtain bias-corrected estimates
##==============================================================================

# Load package
library(lavaan)

# Set seed
set.seed(123)

# Define population model 
pop.model <- ' fx =~ 1*x1 + 0.8*x2 + 0.6*x3 '

# Generate some data
n <- 200
data <- simulateData(model = pop.model, sample.nobs = n)

# Fit model and inspect estimates
fit <- cfa(data = data, model = ' fx =~ x1 + x2 + x3 ')

#-------------------------------------------------------------------------------
# Jackknife bias correction
#-------------------------------------------------------------------------------

# Define function that returns a list of n resamples, each of size n-1 
resample.jackknife <- function(data){
  df.list <- vector(mode = "list", length = nrow(data))
  for(i in 1:nrow(data)) { df.list[[i]] <- data[-i, ] }
  return(df.list)
}

# Create a list of all n jackknife resamples
jack.resamples <- resample.jackknife(data)

# Fit the model to all dataframes in the list
jack.estimates <- semList(model = ' fx =~ x1 + x2 + x3 ',
                          dataList = jack.resamples,
                          store.slots = "partable")


# Inspect coefficients obtained from the first five resamples
coef(jack.estimates)[ , 1:5]

# Compute bias-corrected estimates (equation 10)
jack.corrected <- n*coef(fit) - (n-1)*rowMeans(coef(jack.estimates))

#-------------------------------------------------------------------------------
# Bootstrap bias correction
#-------------------------------------------------------------------------------

# Define function that returns a list of 500 resamples, each of size n
resample.bootstrap <- function(data){
  df.list <- vector(mode = "list", length = 500)
  for(i in 1:500) { df.list[[i]] <- data[sample.int(nrow(data), 
                                                    size = nrow(data),
                                                    replace = TRUE), ] }
  return(df.list)
}

# Create a list of all 500 Bootstrap resamples
boot.resamples <- resample.bootstrap(data)

# Fit the model to all dataframes in the list
boot.estimates <- semList(model = ' fx =~ x1 + x2 + x3 ',
                          dataList = boot.resamples,
                          store.slots = "partable")

# Inspect coefficients obtained from the first five resamples
coef(boot.estimates)[ , 1:5]

# Compute bias-corrected estimates (equation 12)
boot.corrected <- 2*coef(fit) - rowMeans(coef(boot.estimates))

#-------------------------------------------------------------------------------
# Compare estimates
#-------------------------------------------------------------------------------

rbind("ML" = coef(fit),
      "Jackknife" = jack.corrected,
      "Bootstrap" = boot.corrected)
