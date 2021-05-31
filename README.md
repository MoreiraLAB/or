# Opioid Receptors
This is the repository for Opiod Receptor Partners analysis code. By calling following this instructions, the user should be able to run the full pipeline for the opioid receptors profile. In order to make sure everything runs smoothly, the user should follow the upcoming steps:
1. Clone the repository
2. In the **gpcr_variables.py** script, change the **DEFAULT_FOLDER** variable into the users' location of the repository
3. Create a conda environment: `conda create --name opioid_receptors R=3.6 python=3.9`
4. Activate the conda environment: `conda activate opioid_receptors`
5. Add some Python packages: `pip install pandas selenium numpy toolz bs4`
6. Prepare for webscrapping:
	- Confirm your installed Chrome version
	- Download the corresponding chromedriver executable
	- Move the executable to the folder
7. Add some R packages: `RScript  -e "install.packages(c('circlize','stringr','ggplot2','bio3d','tidyverse','ggrepel','ggsci','cowplot','svglite'), repos='https://cran.rstudio.com/')"`
8. Run the pipeline: `python call.py`

# Folder structure
The **pdb** files containing the dimeric structures should be in the root folder. Subdirectories should be:
- results
- processed_results
- images
- summary
- templates
