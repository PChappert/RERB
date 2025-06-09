## Version 0.9.6: June 9, 2025

-   Initial release as a github R package
-   run_Igblast is renamed runAssignGenes
-   Added a beta version of runIgblastn to directly export an AIRR format table from igblast add IMGT gap to V sequences.
-   Added option to switch from runAssignGenes [default] to runIgblastn in most wrapper functions
-   runAssignGenes and runIgblastn now work with TCR and mouse/rat/rabbit/rhesus-monkey datasets
-   Added importFJGates function to import results from FlowJO gating on index-sorting data
-   Added Ab1toAIRR to directly QC and align .Ab1 zipped files from plate-based Sanger sequencing (based on former Ab1toAIRR.sh script now deprecated)
-   Fixed bug in clone frequency calculation in AddAIRRMetadata

## Version 0.9.5: May 6, 2025

-   Base functions: scImportVDJ, scFindBCRClones, reformatVDJinput, resolveMultiContigs, importSangerVDJ, AddAIRRMetadata
-   First ploting function: DonutPlotClonotypes3D, HexbinPlotClonotypes3D, plot_multi_histogram
