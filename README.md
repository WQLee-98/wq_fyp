# Functional consequences of global avian frugivore loss

## Scripts
### Geographic analysis
`retrieveThreats.R`
* Retrieve threat magnitudes for every avian frugivore, calculates missing magnitude scores, and
creates a dataframe with the maximum magnitude scores per threat for every species

`defaunationAnalysis.R`
* simulates defaunation,calculates functional dispersion before and after defaunation, calculates z-scores

`defaunationPlot.R`
* plots maps and graphs

`defaunationDataExploration.R`
* exploring relationship between factors from dataset

`fdFunc1.R`
* Functional dispersion functions used for 'defaunationAnalysis'

`nullModel.R`
* Null model function to calculate 999 randomised functional dispersion values for each cell

### Network analysis
`networkAnalysis.R`
* Calculates bird morphological uniqueness, interaction uniqueness, fits a linear mixed effects model, and plots graphs

`fdFunc2`
* functional dispersion function used for 'networkAnalysis'
