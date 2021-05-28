## Troubleshooting

* The motif of interest or gRNA sequence is not matching
	+ Firstly, make sure that the sequence is in the correct 5'>3' orientation
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
* Error not resolving? Please contact us at mori0164[at]umn.edu
