# Ab1toAIRR

Ab1toAIRR is a wrapper bash script designed to automatically extract sequences from Eurofins plate-based single B cell Sanger sequencing results "platename_SCF_SEQ_ABI.zip" file and perform QC and preliminary VDJ annotation.

## How to use it

Two options are possible:

### **Bash wrapper** (easiest)

To run any number of plate (but also works for one), open a terminal and run the following codes:

-   if running IgG plates [default primer option]:

``` bash
Ab1toAIRR.sh full_filepath_to_ABI.zip_file(1)/platename_SCF_SEQ_ABI(1).zip full_filepath_to_ABI.zip_file(2)/platename_SCF_SEQ_ABI(2).zip full_filepath_to_ABI.zip_file(3)/platename_SCF_SEQ_ABI(3).zip
```

-   to input primers just add it anywhere in the line, order is not important here:

``` bash
Ab1toAIRR.sh IgM full_filepath_to_ABI.zip_file(1)/platename_SCF_SEQ_ABI(1).zip full_filepath_to_ABI.zip_file(2)/platename_SCF_SEQ_ABI(2).zip full_filepath_to_ABI.zip_file(3)/platename_SCF_SEQ_ABI(3).zip
```

-   and you can also choose saving option the same way:

``` bash
Ab1toAIRR.sh IgL png full_filepath_to_ABI.zip_file(1)/platename_SCF_SEQ_ABI(1).zip full_filepath_to_ABI.zip_file(2)/platename_SCF_SEQ_ABI(2).zip full_filepath_to_ABI.zip_file(3)/platename_SCF_SEQ_ABI(3).zip
```

1.  to rapidly copy paste "full_filepath_to_ABI.zip_file/platename_SCF_SEQ_ABI.zip", drag and drop to terminal (doesn't work in Rstudio...) or select file in finder and Command + Option/Alt + C -\> Command + V in Terminal (or Right click and Paste)

2.  **additional info can be automatically added if a "platename_info.xlsx" is present in the same folder** as the "platename_SCF_SEQ_ABI.zip" file. The "platename_info.xlsx" file **must contain a "well_id" collumn** corresponding to the Eurofins sequencing wells (!!!). All other collumn present will be added to the final recap file. You can **use the updated "template_scSangerBCR_info_v1.3.xlsx" file**, key columns (like well_id) and colnames are protected to avoid issues (pwd=u1151, if needed), you can change all other columns and add as many as you want (see Readme inside the excel workbook).

3.  **options for primers are : IgG, IgM, IgK, IgL and MixL (for mixed IgK/IgL plates)**; if no info is provided will default to IgG. Choosing the right primers adapts length cut-offs for QC.

4.  **options for primers are : none, png, html**; by default, will generate "png" QC plots. "html" is more interactive but size is around 25MB per file!

### **R function** (more flexible)

Same options as the bash wrapper with the added possibility to return the final aligned VDJ data table in an AIRR format (when multiple plates, a list of data tables is returned). Also possible to choose which steps to perform (QC, igblast, update_c_call, SHM, full_seq_aa) or to point directly to the wanted database if not installed as recommended below. Simply requires a vector of full file paths as input.

``` r
my_Ab1_file_list <- Ab1toAIRR(files=c("path_to_file1_SCF_SEQ_ABI.zip", "path_to_file2_SCF_SEQ_ABI.zip", "path_to_file3_SCF_SEQ_ABI.zip"), primers = "IgG", save = "png", return_df = FALSE)
```

or

``` r
my_plate <- my_Ab1_file_list <- Ab1toAIRR(files=c("path_to_file1_SCF_SEQ_ABI.zip", "path_to_file2_SCF_SEQ_ABI.zip", "path_to_file3_SCF_SEQ_ABI.zip"), primers = "IgG", save = "png", return_df = TRUE)
```

you can also skip QC if already performed:

``` r
my_plate <- my_Ab1_file_list <- Ab1toAIRR(files=c("path_to_file1_SCF_SEQ_ABI.zip", "path_to_file2_SCF_SEQ_ABI.zip", "path_to_file3_SCF_SEQ_ABI.zip"), QC = FALSE, return_df = TRUE)
```

or skip both QC and further igblast alignment step to simply updated info linked to the excel template:

``` r
my_plate <- my_Ab1_file_list <- Ab1toAIRR(files=c("path_to_file1_SCF_SEQ_ABI.zip", "path_to_file2_SCF_SEQ_ABI.zip", "path_to_file3_SCF_SEQ_ABI.zip"), update_info = TRUE, return_df = TRUE)
```

## What it will do:

1.  unzip zip file in a temporary folder

2.  extract and QC sequences from all .ab1 files, using Eurofins well_ids as unique identifiers (name of .ab1 files should be "your_name-A01.ab1"; "A01" will be extracted in that case); output folder will be "/platename/" and individual well QC plots are exported in "/platename/QC_files/" folder

3.  run AssignGenes.py igblast (changeo) %\> MakeDb.py (changeo) %\> blastn for C region %\> createGermlines (dowser) %\> observedMutations (shazam)

4.  reformat output; **add additional infos if a "zip_file_info.xlsx" file is provides in the same folder as the initial zip file**; add VH somatic mutations information and reconstruct full sequences with and without IgG1 CH1 domain for use in AlphaFold. [Note: could be amended to use identified C chain CH1 domain (to be done in later versions)]

A final "platename_full_recap.xlsx" file is exported with all sequences info (one sheet each for sequences passing all QC steps and sequences failing at individual Blast/QC steps)

**!createGermlines and observedMutations() functions are run before clonal analysis here and should be run again after.**

## Guide to "platename_full_recap.xlsx" file:

In addition to the classical collumns outputted by igblast and Immcantation packages (AIRR format) and the collumns imported from the user provided "platename_info.xlsx", the following collumns will be created:

-   sequence: **trimmed reverse-complemented sequence** with N for bases with Phred score \<10 and two potential base calls [avoid artificially counting them as somatic mutations]; trimming is performed using the Mott's trimming algorithm as implemented in the sangeranalyseR package (<https://sangeranalyser.readthedocs.io/en/latest/>), see also the QC_plot_png file for each well;

-   primary_sequence: trimmed reverse-complemented sequence with primary base calls.

-   raw_sequence: untrimmed sequence as provided by Eurofins;

-   primer: primer info provided by user;

-   sequence_length: length of the trimmed sequence;

-   pct_under_30QC_in_trimmed: percent of bases in the trimmed sequence that have a Phred confidence score below 30;

-   low10QC_alternate_calls: list of base in the trimmed reverse-complemented sequence with a Phred confidence score below 10(\<10% of chance of having the correct base call). format = primary base call/position in trimmed reverse-complemented sequence/secondary base call;

-   low30QC_alternate_calls: list of base in the trimmed reverse-complemented sequence with a Phred confidence score between 10 and 30 (between 10 and 0.01% of chance of having the correct base call). format = primary base call/position in trimmed reverse-complemented sequence/secondary base call.

-   c_call and c_call_igblast; c_call correspond to the output of our homemade Blastn pipeline to call c gene (with no more selection based on alignment length, instead multiple calls are return, i.e. IGHG1\|IGHG2\|IGHG3); c_call_igblast correspond to the output of igblast.

-   missing_v_bp: number of bases missing at the beginning of the aligned V;

-   missing_j_bp: number of bases missing at the end of the aligned V;

-   full_sequence: reconstructed full sequence in nucleotide format, with missing beginning of V and end of J reverted to germline;

-   full_sequence_aa: reconstructed full sequence in amino-acid format, with missing beginning of V and end of J reverted to germline;

-   full_sequence_fab_aa: reconstructed full sequence in amino-acid format with added CH1 domain from IgG1, with missing beginning of V and end of J reverted to germline;

-   comments: comments on the full sequence reconstruction.

Finally, **rows for sequences with partial length issues are highlighted in light orange** (aligned V \< 270 \| missing too much base pairs in V for reconstruction \| full length \< tight cutoffs for each type of VDJ sequences). And **rows for sequences with partial quality issues are highlighted in light orange with text in red** (\>10% of base calls with Phred score \< 30). These should be checked manually.

## How to install:

1.  Install RERB first, including standalone versions of **blast** and **igblast**, **Immcantation ChangeO**, **Alakazam**, **Shazam** and **Dowser** packages (<https://changeo.readthedocs.io/en/stable/examples/igblast.html>) and blastable IMGT database.

2.  To use the bash script: modify line 59 of the **Ab1toAIRR.sh** script to adapt it to the path to your RERB folder (default: devtools::load_all("\~/R_packages/RERB")). Then copy the script to the /usr/local/bin folder (or any other folder in your \$PATH, just change the first line of the following script accordingly) and to make bash script executable, run in Terminal:

``` bash
cd /usr/local/bin
chmod 755 Ab1toAIRR.sh
```

6.  If not already installed upon RERB installation, install the following R packages:

``` r
install.packages("optparse") 
install.packages("tidyverse") 
install.packages("ggrepel") 
install.packages("openxlsx")

if (!require("BiocManager", quietly = TRUE)) 
install.packages("BiocManager") 
BiocManager::install("Biostrings") 
BiocManager::install("sangeranalyseR") 
BiocManager::install("sangerseqR")  
```

5.  "png" export of images generated through the plotly package in sangeranalyseR requires the installation of the reticulate package in R and the kaleido/plotly packages in python as follow :

    *[Note] if this isn't enough to ensure reticulate uses the correct python env (where plotly is installed), you can use the reticulate_py_env argument when using the Ab1toAIRR function in R or add RETICULATE_PYTHON="/Users/yourname/Library/r-miniconda-arm64/envs/r-reticulate/bin/python" to your .Renviron file.*

``` r
install.packages('reticulate') 
reticulate::install_miniconda() 
reticulate::conda_install('r-reticulate', 'python-kaleido') 
reticulate::conda_install('r-reticulate', 'plotly', channel = 'plotly') 
reticulate::use_miniconda('r-reticulate')
```

6.  finally, "html" export of images generated through the plotly package in sangeranalyseR requires the htmlwidgets package to be installed:

``` r
install.packages("htmlwidgets")
```
