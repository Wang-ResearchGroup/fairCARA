# fairCARA
To install this package in R, run the following commands:

```R
library(devtools) 
install_github("Wang-ResearchGroup/fairCARA")
```

This package implements fair covariate-adjusted response-adaptive (CARA) design. Our design:
- Improves data use efficiency;
- Achieves an “envy-free” treatment assignment guarantee;
- Ensures the welfare of participating subpopulations. 

Our package has the following features:
- Simulate covariate-adjusted response-adaptive experiments;
- Implement the proposed fair CARA design;
- Comparison with conventional designs.

In the following example, we show how to simulate adaptive experiments using our design. 

## Example usage:

```R

# load packages
library(fairCARA) 
library(nloptr)

# user-specified subgroup proportions
true_p <- c(0.5,0.5)

# specify mean and standard deviation in treatment and control arm for adaptive experiments simulation
mu1.t <- 4  
mu1.c <- 2  
sd1.t <- 2.5  
sd1.c <- 1.5 

# specify parameters for the second subgroup
mu2.t <- 1 
mu2.c <- 4    
sd2.t <- 1.2
sd2.c <- 3.5 

# true average treatment effects
tau_true <- c(2,-3)
tau_true_all <- 2*0.5+(-3)*0.5

# initialize the experiment with 100 subjects 
m <- 100

# sequentially enroll 200 subjects
n <- 200

# specify the maximum difference of subgroup treatment assignment probabilities
envy_c <- 0.1

# simulate adaptive experiments and obtain outputs
results <- fairCARA(m,n,envy_c)
                  
```

### Reference
- Wei, W., Ma, X., and Wang, J. (2023). Fair Adaptive Experiments.
