# Using R 
1. Load the script `01-download-ereefs-data.R` into RStudio. 
2. Set the working directory to the folder corresponding to the `01-download-ereefs-data.R` script. Use the Files tab to navigate to the script folder, the using the `More / Set As Working Directory`. If you know what the directory path is you can also use the `setwd` command directly in the Console such as:
```R
setwd("~/2025/ereefs/GBR_AIMS_eReefs-climatology_2025")
```
3. Run the download script
```R
source("01-download-ereefs-data.R")
```

# Setting up Python
1. Create the Conda environment. 
    ```bash
    cd {path to the GBR_AIMS_ereefs-climatology_2025 dataset} 
    conda env create -f environment.yml
    ```
2. Activate the environment
    ```bash
    conda activate ereefs-climate
    ```