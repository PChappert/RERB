## Version 0.9.7 - devel: Oct, 2025

-   new COVID dataset included
-   new set of functions to interact with AlphaFold3 server for antigen-antibodies interaction predictions
-   new CircosClonotypes function
-   updated reconstructFullVDJ to deal with Ns added at the createGermlines step of Dowser and remaining gaps in the J segment.
-   updated plotVGenePairing to had possibility to look at the family, gene or allele levels
-   plotCDR3logo is now rename CDR3logos and SingleCDR3logo and allows plotting of CDR3 sequenecs from different clone (and difference sizes) using ClustalW alignments.
-   updated Ab1toAIRR to solve issue with reticulate connection
-   fixed a few issues when merging Sanger plates with different metadata
-   updated bcr_info/tcr_info in addAIRRMetadata
-   updated runBlastnC() to filter out all alignment <25 pb

## Version 0.9.6: June 29, 2025

-   Initial release on github
-   added return_plot options to DonutPlotClonotypes
-   new plotCDR3logo, plotVGenePairing and CircosClonotypes (beta) plotting functions
-   All functions have been updated to be TCR compatible (!still beta version for the TCR part)
-   flagBdoublets and flagNonBdoublets renamed as homotypivVDJdoublets et heterotypicVDJdoublets
-   DonutPlotClonotypes3D and HexbinPlotClonotypes3D are renamed DonutPlotClonotypes and HexbinPlotClonotypes and underlying DonutPlotClonotypes and HexbinPlotClonotypes function renames as SingleDonutPlotClonotypes and SingleHexbinPlotClonotypes.
-   plot_multi_histogram is renamed plotHistogram
-   run_Igblast is renamed runAssignGenes
-   Added a beta version of runIgblastn to directly export an AIRR format table from igblast add IMGT gap to V sequences.
-   Added option to switch from runAssignGenes [default] to runIgblastn in most wrapper functions
-   runAssignGenes and runIgblastn now work with TCR and mouse/rat/rabbit/rhesus-monkey datasets
-   Added importFJGates function to import results from FlowJO gating on index-sorting data
-   Added Ab1toAIRR to directly QC and align .Ab1 zipped files from plate-based Sanger sequencing (based on former Ab1toAIRR.sh script now deprecated). Provides added flexibility regarding which steps need to be run a second time to simply update alignment or imported well_id-associated info.
-   Fully updated logs throughout using a new time_and_log fonction
-   Fixed bug in clone frequency calculation in AddAIRRMetadata

## Version 0.9.5: May 6, 2025

-   Base functions: scImportVDJ, scFindBCRClones, reformatVDJinput, resolveMultiContigs, reconstructFullVDJ, importSangerVDJ, AddAIRRMetadata
-   First ploting function: DonutPlotClonotypes3D, HexbinPlotClonotypes3D, plot_multi_histogram
