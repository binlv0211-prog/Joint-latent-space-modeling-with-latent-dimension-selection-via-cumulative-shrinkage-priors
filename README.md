### **Description**

The repository contains code associated with the following article (i.e., the code supporting the numerical findings)  and describes how to reproduce the findings.
> Bin, lv. and Siliang, Zhang. (2025). Joint Latent Space Modeling with Latent Dimension Selection via Cumulative Shrinkage Priors.
> 
### **Core Functionality**

The primary functions of the code in this repository are:

*   **Data Simulation and Evaluation:** A dedicated script (`data_process.R`) for generating simulated datasets and for the calculation of evaluation criteria.
*   **Gibbs Sampling Algorithms:** A comprehensive suite of scripts to perform Gibbs sampling for various data types, including:
    *   Network data only (`network_only.R`)
    *   Continuous node variable data only (`Y_only_normal.R`)
    *   Binary node variable data only (`Y_only_binary.R`)
    *   Network with continuous node variable data (`network_normal.R`)
    *   Network with binary node variable data (`network_binary.R`)
*   **Cross-Validation and Prior Variations:** Scripts that incorporate normal priors and functions for cross-validation in the context of:
    *   Binary node variables in network data (`network_normal_kl.R`)
    *   Continuous node variables in network data (`network_binary_kl.R`)
*   **Missing Data Imputation:** Specialized functions for imputing missing binary node variables using:
    *   Network latent space modeling (`network_miss.R`, `network_binary_miss.R`)
    *   The Multiple Imputation by Chained Equations (MICE) algorithm (`mice.R`)

### **Replication Scripts**

The repository includes high-level scripts to reproduce the main findings of the paper:

*   **Simulation 1 (`simulation1.R`):** Investigates the impact of sample size on model performance.
*   **Simulation 2 (`simulation2.R`):** Explores the effects of network density in a six-node simulation setup.
*   **Real Data Analysis 1 (`real_data1.R`):** Applies the proposed methodology to the French financial elites dataset.
*   **Real Data Analysis 2 (`real_data2.R`):** Demonstrates the method's application on a Facebook social network dataset.

### **Code Format**

All scripts are written in the R programming language. To execute the code, users will need a functioning R environment with the necessary packages installed.
