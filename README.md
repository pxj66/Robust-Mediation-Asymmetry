# RMSkew: Robust Mediation Analysis with Skewed (Asymmetric) Data
An R package implements a new robust linear regression method for conducting mediation analysis to handle skewed data.
## Installation
install.packages("RMSkew_0.0.0.9000.tar.gz", 
                 repos = NULL, 
                 type = "source")
## Main Functions
- `rlmskew` (a robust linear model for skewed data) implements a robust regression method with asymmetric and symmetric Huber, Tukey losses.
- `extract_mediation_formula` can extract the mediation model from an R formula similar to this:   `y~x+m(med1, medi2)+covariates(cov1,cov2)`.
- `fit_skew_mediation` fits the mediation model specified by the **mediation formula**.
- `test_skew_mediation` test the indirect, direct, total effects.
