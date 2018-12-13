###bamLiftOver - Lift BAM alignments from one build to another (not just the coordinates!)\

bamLiftOver requires the same input as the standard liftover tool, i.e. a chain file from the current coordinate system to the desired one and of course a bam file containing the alignments to lift over. \
It performs coordinate conversion, CIGAR string updates and fragment size corrections for all read pairs along with updating the bam file header according to the supplied chain file.\

**command:** `bamLiftOver alignmentsInOldBuild.bam(input) chain-file.chain(input) alignmentsInNewBuild.bam(output)`\

Optionally, you can use mergeHaps to merge two BAM files containing the same reads aligned to 2 different references(haplotypes for example) after lifting them back to the reference space.\
mergeHaps selects the best alignment(based on alignment score, AS tag) for each read/read-pair and writes those to primary.bam along with a new tag XH:1/2 indicating which file the alignment came from (0 if they are equal).\
Other files generated by mergeHaps are:\
hap1.specific.bam -> Contains all alignments from file1 without an equivalent alignment in file2.\
hap2.specific.bam -> Contains all alignments from file2 without an equivalent alignment in file1.\
equalSecondary.bam -> Contains all alignments equal for both files without being the best alignments.\
A third parameter can also be supplied to mergeHaps in the form of a BAM file containing reads aligned directly to the build the others were lifted to.\
mergeHaps then generates summary statistics regarding between the litfted over alignments and the direct alignment (also based on alignment score).\

**command:** `mergeHaps 1.liftover.bam(input) 2.liftover.bam(input) [refAligned.bam](input)`\

Compilation commands:\
`make bamLiftOver`\
`make mergeHaps`\
Requires htslib, available `[here](https://github.com/samtools/htslib)` and SeqAnHTS(included as a submodule, git cloning must be done recursively)\
