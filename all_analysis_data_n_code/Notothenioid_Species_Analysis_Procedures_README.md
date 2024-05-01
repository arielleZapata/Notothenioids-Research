## Overview
This document outlines the methods and procedures utilized for a comprehensive analysis of Notothenioid species using various R packages and data files. The project encompasses a wide array of analyses, including geometric morphometric analysis, phylogenetic data processing, depth data integration, and morphospacial visualization.

### Software and Libraries
- **R Studio**
- **R Libraries**:
  - `ggplot2`: For creating complex and customized plots.
  - `ggdist`: Enhances ggplot2 for distributional representations.
  - `tidyquant`: Provides themes and color palettes for financial analysis, used here for aesthetic enhancements.
  - `phytools`: For analyzing and visualizing phylogenetic data.
  - `geomorph`: For geometric morphometrics analysis.
  - `dplyr`: For data manipulation.
  - `geiger`: For phylogenetic analysis.
  - `splancs`: For spatial and cluster analysis.
  - `gridExtra`: For arranging multiple grid-based plots on a page.

### Data Files
- Phylogenetic Tree (`notothenioid_timetree.tre`): Represents evolutionary relationships among Notothenioid species.
- Landmark Data (`all_tps_file.tps`): Contains geometric morphometrics data.
- Species Names Update List (`updated_names.csv`): Lists original and updated species names.
- Depth Data (`CombinedAntarcticData.csv`, `Polarstern_2022_Fish_Sampled.csv`): Includes depth information for fish specimens.
- Name Replacement for Depth Data (`FIXEDdepthNames.csv`): For reconciling species names.
- PCA Data (`pca.csv`): Contains PCA scores for each specimen.
- Various RDS files for depth ranges (e.g., `depth.0TO100_OUTPUTS.RData`).

## Procedures

### Part 1: Initial Data Analysis
- **Environment Setup**: Load libraries (`geomorph`, `phytools`) and set working directory.
- **Data Preparation**: Read and clean the phylogenetic tree and landmark data.
- **Geometric Morphometric Analysis (GPA)**: Conduct GPA and visualize the consensus configuration.
- **Phylogenetic Data Processing**: Convert and rename species in the phylogenetic tree.
- **Principal Component Analysis (PCA)**: Perform PCA and visualize results.
- **Phylomorphospace Analysis**: Plot data combining morphological and phylogenetic information.
- **Output Generation**: Export PCA components, trees, and plots.

### Part 2: Depth Data Analysis
- **Load Libraries**: `dplyr`, `phytools`, `geiger`, `splancs`.
- **Import Functions and Data**: Source custom functions and read depth data.
- **Data Name Reconciliation**: Replace species names using mappings.
- **Phylogenetic and Morphometric Integration**: Use the pruned tree with depth data for analysis.
- **Output Generation**: Save filtered data and visualizations of convex hull areas.

### Part 3: Advanced Morphospacial Visualization
- **Setup Environment**: Load `ggplot2`, `phytools`, `gridExtra`.
- **Import Data**: Read the pruned phylogenetic tree and PCA data.
- **Create Morphospaces**: Utilize custom functions for different depth ranges.
- **Generating Facet Plots**: Convert to ggplot and arrange plots in facets.
- **Output Generation**: Export to PDF all the facet plots.

### Part 4: Visualization of Statistical Distributions
- **Setup Environment**: Load `ggplot2`, `ggdist`, `tidyquant`.
- **Data Preparation**: Import and group depth data.
- **Raincloud Plot Creation**: Visualize the distribution of lengths across depth groups.
- **Output Generation**: Export the raincloud plot to PDF.

### Part 5: Analysis of Interquartile Ranges (IQR)
- **Setup Environment**: Load `ggplot2`.
- **Data Import and Setup**: Import depth data subsets.
- **Calculate and Visualize IQR**: For 100s and 200s depth ranges.
- **Output Generation**: Export IQR plots to PDF.

### Part 6: Variance Analysis
- **Setup Environment**: Load `ggplot2`, `dplyr`.
- **Data Preparation**: Read depth data subsets.
- **Variance Calculation and Visualization**: For 100s and 200s depth ranges.
- **Output Generation**: Visualizations exported to PDF.

### Part 7: Integrated Phylogenetic and Morphometric Analysis
- **Setup Environment**: Load `rfishbase`, `phytools`.
- **Data Integration and Cleaning**: Match PCA to fish data and update the phylogenetic tree.
- **Phylogenetic and Morphometric Analysis**: Generate and visualize phylomorphospace plots.
- **Output Generation**: Save integrated data and visualizations.

### Part 8: Depth and Length Visualization
- **Visualization of Species by Depth**: Create density ridge plots using `ggplot2` and `ggridges`.
- **Data Manipulation for Length Analysis**: Prepare data for visualization.
- **Output Generation**: Export plots to PDF files.
