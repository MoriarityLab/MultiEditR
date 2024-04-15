## Frequently-Asked Questions

> What phred cutoff should I use?

We find 0.001 is usually sufficient, but for lower quality sequences, 0.01 may be better.

> In the results, what is "target_base"

The index of the base potentially edited, relative to the position of the start of the motif.

> Why isn't the motif direction auto-detected?

We considered this, but decided it was better to require the user to know the orientation of the motif relative to the sample. Right now, we will auto-revcom the control file, however.

> My control sequence was a Sanger file. Why didn't it's chromatogram show up?

If the control sequence needed to be rev-com, we do not show the control sanger sequence.

> When I upload the next file, the previous one I uploaded disappears!

You must upload all sample and control files in one step. This is why we recommend putting them in a single directory.

> I want to edit the plots made. Can I get R code to do that?

This shiny app is based on multiEditR, which you can install from github and use from the command line: the package is here https://github.com/MoriarityLab/multiEditR.pckg and instructions to install are on the main page. 

> I'm getting a weird message about indels.

Sometimes during alignment of the control to the sample, indels are detected, especially if the base quality is low. Sometimes these can be fixed by increasing the phred_cutoff, but the message may suggest that there are actually indels between your sample / control.

> Your old version of the editR app didn't require a control sequence, but now you do. I don't have a control sequence, what do I do?

You can make a control fasta file manually and upload that. A fasta file looks like this:

```

> name_of_control_sequence
ATGCCGTTAACTACTAGGACACAGGATTCANNGATCCATA

```

> I can't figure out how to get my parameters.xlsx and my sequence files to work, what should I do?

Try the example dataset which is a download at the bottom of the instructions page. That should work, so once that gets working, use it as a guide. 
