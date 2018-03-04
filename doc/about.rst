About SAVNet
============

Genomic variants causing splicing aberations are an important class of functional variants and playing important roles in pathogenesis in many diseases including cancer.
To identify splicing associated variants, using matched transcriptome sequencing data would help to see the consequences of these variants in terms of transcription.
However, even when we have matched transcriptome sequencing data, there are number of challenges.
1. Transcriptome sequencing reads are usually short (50 ~ 150bp). So we can only see the part of splicing changes. Also, transcription is noisy biological process, and there are many variations among tissues as well as individuals. 
2. Most of somatic variants are sporadic (usually just one or so in the cohort at nucleotides levels). So splicing QTL (sQTL) approaches, which have been successfully used to identify splicing causing common variants in population studies, cannot be used in their current forms.
3. One genomic variant can generate multiple transcipt isoforms. Also, the consequences of several distinct genomic variants can converge to common splicing changes.

To overcome these challenges we have developed a novel statistical approach, SAVNet, which can sensitively and accurately identify splicing associated variants from genome and transcriptome data. 
The basic idea of extracting SAVs is to identify the pair of genomic variants and splicing changes occuring specific to individuals with these variants. However, as we noted above, naive approach would fail because of the noise of transcription and complexity of genomic variants and splicing changes.

The keys to success of SAVNet is that:
1. We have greatly restricted to the class of "association" between genomic variants and splicing changes; splicing donor/acceptor disruption and splicing donor/acceptor creation. By doing this. the false discovery ratio of SAVNet is below 5% for most cohort.
2. To solve the problem of redundant relationships, we have modeled the relationships as network, where after creating network with maximal links considering all the possible links, we prune links and extract the most reliable network that can explain the mutation status and supporting read number of splicing changes most effectively. Here, we utilized the Bayesian Model Averaging framework to perform statistical test in this complex settings.

SAVNet is developed mainly for cancer genome and transcriptome sequencing data. However, we believe that this can be also sed for identifying germline "rare" variants that cause splicing changes.
Also, we are currently improving the framework so that it can extract genomic structural variations causing splicing changes and so on.


Class of association

Here, we use the term splicing donor site as regions around the end of exons, and the start of exons. For the donor site, we set last 3bp in the last exonic regions and first 6bp in the intronic regions. For acceptor sites, we set the last 6bp in the acceptor sites and the first 1bp of exons. These values are determined by checking false discovery ratios when performing SAVNet. Also, we call the first 2bp and the last 2bp of introns (GT-AG) as "cannonical splice sites".


Splicing donor/acceptor disruption

When the genomic variants are affecting the annotated splicing donor/acceptor sites and splicing aberations are overlaping this somatic mutations (allowing some amount of mergins, 50bp by default). We nominate these relationship as splicing site disruption. 

Splicing donor/acceptor creation

When unnanotated edges of Alternative 5' or 3'SS are overlapping with the genomic variants, we nominate these relationship as splicing site creations


Classification of splicing changes

We consider "Exon skipping", "Alternative 5' splice site (SS)", "Alternative 3'SS" and "Intron retention" for the splicing aberattions caused by genomic variants. The first three types of aberrations are collected by split-alignment (therefore considered as introns). Using control panels (set of non-disease transcriptome data), we identify "normal" splicing aberations and remove them to improve the accuracy.
 

Exon skipping: Both the edges of introns are annotated exon-intron junctions. However, there is no annotated intron supporting in the transcriptome database (we use RefSeq).

Alternative 5'SS: One of the edge of introns corresponds to the annotated splicing acceptor site. When unnanotated edge is located in intronic region (>= 30bp from the annotated splicing donor site), we identify them as Intronic alternative 5'SS. 
Alternative 3'SS: One of the edge of introns corresponds to the annotated splicing donor site.

For identifying intron retention, we count the number of short reads covering both exon and intron regions (>= 8bp each in the default setting) at each annotated splicing donor/acceptor site.



