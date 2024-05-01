# Notothenioid Research Analysis Repository

This repository contains all the analysis data and code for studying the evolutionary patterns, morphological variations, and depth adaptations of Notothenioid fish species. The repository is structured into multiple folders, each dedicated to a specific aspect of the research process.

## Repository Structure

- **all_analysis_data_n_code** - Contains scripts, data processing procedures, and analyses related to phylogenetic and morphometric studies.
- **all_data_Notothenioids** - Hosts the initial datasets including phylogenetic trees, depth data, and morphometric data of Notothenioids.
- **figures_n_designs** - Comprises code and outputs for generating figures and design elements used in the research publications.
- **.gitignore** - Specifies intentionally untracked files to ignore.

## Contents and Usage

### 1. all_analysis_data_n_code

#### Materials
- **Software and Libraries:**
  - R Studio
  - R Libraries:
    - ggplot2 (3.4.4)
    - ggdist (3.3.1)
    - tidyquant (1.0.7)
    - phytools (2.0-3)
    - geomorph (4.0.6)
    - dplyr (1.1.4)
    - geiger (2.0.11)
    - splancs (2.01-44)
    - gridExtra (2.3)

#### Data Files
- Phylogenetic Tree (`notothenioid_timetree.tre`)
- Landmark Data (`all_tps_file.tps`)
- Species Names Update List (`updated_names.csv`)
- Depth Data (`CombinedAntarcticData.csv`, `Polarstern_2022_Fish_Sampled.csv`)
- Name Replacement for Depth Data (`FIXEDdepthNames.csv`)
- PCA Data (`pca.csv`)
- RDS Files for Depth Ranges (e.g., `depth.0TO100_OUTPUTS.RData`, `depth.101TO200_OUTPUTS.RData`)

### Procedures
Detailed steps from data setup, analysis to output generation are outlined within each script in the `all_analysis_data_n_code` folder. These include:
- Environment setup
- Data preparation and cleaning
- Advanced morphometric and phylogenetic analyses
- Visualization and output generation

### 2. all_data_Notothenioids
Initial raw and processed data are maintained here for reference and reproducibility.

### 3. figures_n_designs
Adobe illustrator files for figure adjustments and alignments specific to publication needs are found here.

## How to Use This Repository
1. Clone the repository to your local machine.
2. Ensure that R Studio is installed along with all required libraries.
3. Follow the procedures outlined in the scripts within the `all_analysis_data_n_code` folder to perform the analyses.
4. Generated outputs including plots, processed data files, and PDFs are saved within designated output folders as specified in the scripts.

## Contributing
Contributions to this research are welcome. Please fork the repository and submit pull requests with suggested changes. Ensure that any additions are well-documented and tested with the existing dataset.


## Contact
For any further queries regarding the research or data usage, please contact the repository maintainers.
