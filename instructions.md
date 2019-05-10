## Welcome to MultiEditR version 1.0!

We hope that our program can be of use for your edit detection and quantification needs. Please note that this is a beta version of the application, and as such aspects of it are still under development. However, if you are noticing any errors when using the application please feel free to reach out Mitch at klues009@umn.edu for troubleshooting and input on the application.

## What this is

The purpose of this software is to detect multiple nucleic acid editing events from either RNA or DNA samples. You provide a sample .ab1 of the region of interest, an control unedited .fasta sequence or .ab1 file of the region, as well as the motif and edit you are interested in studying. MultiEditR will isolate your motif(s) of interest in the sample and report whether or not significant editing is detected as well as a measurement of the percent editing in the motif of interest. 

### Citation

If our software helps you -- please cite us!

Mitchell Kluesner, Annette Arnold, Taga Lerner, Rafail Nikolaos Tasakis, Sandra Wüst, Marco Binder, Branden Moriarity, and Riccardo Pecori **MultiEditR:  An easy validation method for detecting and quantifying RNA editing from Sanger sequencing** *Biorxiv*

## Running the prediction software

To run this software, please provide:

*   A sanger sequencing file of your region of interest
*   A control fasta file containing the sequence without any edits or a control sanger sequencing without any edits
*	A motif of interest (e.g. YAR, TCA, YURACN, or GGCGCGGGCCGCTCGCTCTA)
*	A WT base of interest (e.g. A, C, G, or T)
*	Edits of interest, which if theres multiple separated by "|" (e.g., T, T|G, T|G|A)

You can also provide a custom P-value cutoff for calling base-editing and Phred score, however we recommend using the default parameters

If you want to see some output of the software, just check the “Load example data” checkbox.

### Special thanks to / sources

Hill JT, Demarest BL, Bisgrove BW, Su YC, Smith M, Yost HJ. (2014) **Poly Peak Parser: Method and software for identification of unknown indels using Sanger Sequencing of PCR products.** *Developmental Dynamics.* [PMID: 25160973](http://www.ncbi.nlm.nih.gov/pubmed/25160973)

Brinkman EK, Chen T, Amendola M, van Steensel B. (2014) **Easy quantitative assessment of genome editing by sequence trace decomposition.** *Nucleic Acids Research.* doi: 10.1093/nar/gku936

Derek Nedveck

