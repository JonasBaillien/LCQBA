# LCQBA
R-code accompanying "Flexible asymmetric multivariate distributions based on two-piece univariate distributions"

We recommend installing the latest version of R to make sure all functions run properly. This code will later be made available in an R-package for quick and easy use.

For theoretical result on which this code is based, we refer to
  - Baillien, J et al. (2021). Flexible asymmetric multivariate distributions based on two-piece univariate distributions. Under review.
  - Gijbels, I. et al. (2019). On quantile-based asymmetric family of distributions: Properties and inference. International Statistical Review, 87(3):471–504.

Required libraries:
  "DAAG"
  "sn"
  "mrfDepth"
  "fields"
  "doParallel"
  "combinat"
  "nloptr"
  "optimx"
  "utils"
  "mvnTest"
  "mixtools"
  "LaplacesDemon"
  "QBAsyDist"
  "expm"
  "alabama"
  
The repository consists of 3 main parts:
  - Functions used in either simulation studies or practical examples. These will be made available either in a new R-package, or as an extension to the "QBAsyDist" package.
  - Simulation studies. Results from simulation studies need to be stored locally. Directories to store and bundle these outputs need to be changed to their local counterpart. Run times might vary depending on the chosen model.
  - Examples and plots. For examples, the data needs to be loaded in. Same holds for plots of simulations. Bundles output is loaded in from a local directory.

Comments are added to the functions, please read these to ensure proper use.
