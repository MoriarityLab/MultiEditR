# MultiEditR v1.0.7
 
The purpose of this software is to detect and quantify base editing events in RNA or DNA. Simply, provide an editing motif or guide RNA protospacer sequence and a .ab1 Sanger sequencing file of an edited region of interest (~300-700bp). MultiEditR will then generate distributions of the expected level of noise based on your sample. From these expectations, EditR assigns a probability that each background peak in your target regions are from noise -- in other a words a *P*-value for it being edited as opposed to noise.

## How to use MultiEditR

1. Upload an .ab1 sample file for a region of interest
2. Upload an .ab1 or .fasta control file that contains the region of interest without any edits
	+ Make sure the sample file and control file are in the same orientation
2. Enter your motif of interest, gRNA protospacer sequence, or any sequence of interest.
	+ If your gRNA is antisense to the ctrl .ab1 file, make sure to reverse complement it first
3. Enter the WT base suspected to be edited
4. Enter the edited base, if multiple editing outcomes are possible separate each base with "|"
	+  For example, C > T or G can be specified as T|G
4. On the top of the page, click the "Analysis" tab
5. Check the "Raw sample signal by position" to see if your sample was appropriately trimmed to a region of interest
	+ If the sample appears to be **over-trimmed, increase the phred cutoff value** by an order of magnitude in the left panel
	+ If the samples appears to be **under-trimmed, reduce the phred cutoff value** an order of magnitude in the left panel
6. Look at the various plots under the "Analysis" tab
	+ The "trimmed sample noise by position" plot shows where motifs or gRNA sequences are located n the sample
		- Colored dots represent significant edits, while translucent dots are non-significant editing events
	+ The "barplot of editing data" shows the distribution of edits vs. non-significant noise
	+ Subsequent summary tables provide descriptive statistics of the sample
7. On the top of the page, click the "Chromatograms" tab to view the chromatogram at a specific motif instance
	+ To look close at a specific motif, hover over the bottom "Edits and their significance" plot to get the "motif_id" of a specific instance
	+ Enter that motif id into the left panel under "Enter desired motif instance to examine" to update the chromatograms of the sample and control files
	+ To change the region included in the chromatogram, change the "Enter desired padding to chromatogram" input as needed. This parameter dictates how many bases surrounding the motif of interest are included in the chromatograms
8. To view the raw data of the analysis, click the "Raw Data" tab
	+ This data can be downloaded by clicking on the "Download Raw Data" in the left side panel

## Citation

If our software helps you -- please cite us!

Mitchell Kluesner, Annette Arnold, Taga Lerner, Rafail Nikolaos Tasakis, Sandra Wüst, Marco Binder, Branden S. Moriarity§, Riccardo Pecori§  [MultiEditR: An easy validation method for detecting and quantifying RNA editing from Sanger sequencing](https://www.biorxiv.org/content/10.1101/633685v1). *bioRχiv*. *2019*.

§ These corresponding authors contributed equally.

### Special thanks to

Derek Nedveck, M.S. for his work and mentorship on the original EditR application

Hill JT, Demarest BL, Bisgrove BW, Su YC, Smith M, Yost HJ. (2014) **Poly Peak Parser: Method and software for identification of unknown indels using Sanger Sequencing of PCR products.** *Developmental Dynamics.* [PMID: 25160973](http://www.ncbi.nlm.nih.gov/pubmed/25160973)

Brinkman EK, Chen T, Amendola M, van Steensel B. (2014) **Easy quantitative assessment of genome editing by sequence trace decomposition.** *Nucleic Acids Research.* doi: 10.1093/nar/gku936

## Troubleshooting

* The motif of interest or gRNA sequence is not matching
	+ Firstly, try reverse complementing the sequence
	+ The .ab1 file may be overly noisy, check the "Analysis" tab first, if so you may need to re-sequence your sample
	+ If there are SNPs within your region of interest, noisy sequencing or other artifacts, replace ambiguous bases with "N"s in your sequence of interest
		- MultiEditR can use any IUPAC nucleotides for the gRNA input
		- Additionally, the gRNA sequence length can be extended indefinitely to subset a continous region of interest (e.g. a 30-50 bp region)
		- However, note that noisy sequencing will yield unreliable results
* When uploading your .ab1 file the error `<!DOCTYPE HTML><html_lang="en"><head><title>Error</title><style type="text/css">` is given
	+ Open the .ab1 file in a .ab1 file viewer such as SnapGene Viewer
	+ Manually delete the 5' and 3' end of the file to create a region of relatively clean sequencing on either side
	+ Save this edited file as a new .ab1 file, then re-upload to MultiEditR
* Your sequencing is too noisy to detect an edit
	+ Try gel extracting your PCR product
		- Before sending your PCR product for sequencing load 25-50 uL of PCR reaction on a 1% agarose gel in TAE buffer
		- Run the gel and excise your band, and then purify using a gel extraction kit (Qiagen QIAquick works well)
		- Elute your sample in water and send to a trusted Sanger sequencing provider (ACGTinc, single pass works well)
