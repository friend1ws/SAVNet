About
=====

Introduction
------------

Splicing associated variant (SAV), genomic variant driving splicing alterations,
is an important class of functional variants.
SAVs often play crucial roles in pathogenesis in many diseases including cancer.
To identify SAVs, using matched transcriptome sequencing data
(collected from the patient having the genomic variant of interest) would help to see
the consequence of each genomic variant in terms of transcription.
However, even when transcriptome sequencing data is available,
there are number of challenges to sensitive and accurate identification of SAVs.

1. Diversity of transcription

  There are huge number of variations in many transcripts
  among tissues as well as individuals.
  In addition, transcription is a noisy process,
  and a number of somewhat immature forms of transcripts can be observed in transcriptome sequencing data.
  Furthremore, since the lengths of short reads of transcriptome sequencing data is short (50bp ~ 150bp).
  a number of short reads can be mis-aligned, increasing the noise.

2. Rarity of somatic variants

  Most of somatic variants are sporadic
  (usually just one or so in the cohort at nucleotides levels).
  So splicing QTL (sQTL) approaches,
  which have been successfully used to identify common SAVs in population studies,
  show low statistical power, and cannot be used in their current forms.

3. Redundant relationships between genomic variants and splicing changes

  One genomic variant can generate multiple transcipt isoforms.
  Also, the consequences of several distinct genomic variants can converge to common splicing changes.

To overcome these challenges,
We have developed a novel statistical approach, SAVNet,
which can sensitively and accurately identify splicing associated variants
from genome and transcriptome data.
The basic idea of extracting SAVs is to identify the pair of genomic variants and splicing changes occurring specific to individuals with these variants.
However, as we noted above, naive approach would fail.

The keys to success of SAVNet is that:

- Restriction of the class of association

  We have greatly restricted to the class of **association**
  between genomic variants and splicing changes;
  splicing donor/acceptor disruption and splicing donor/acceptor creation (which will be explained below).
  By doing this. the false discovery ratios of SAVNet are below 5% for most cohort.

- Network modeling

  To solve the problem of redundant relationships,
  we have modeled the relationships as network,
  where after creating network with maximal links considering all the possible links,
  we prune links and extract the most reliable network that can explain
  the mutation status and supporting read number of splicing changes most effectively.
  Here, we utilized the Bayesian Model Averaging framework to perform statistical test
  in this complex settings.

SAVNet is developed mainly for cancer genome and transcriptome sequencing data.
However, we believe that this can be also sed for identifying germline **rare** variants
that cause splicing changes.
Also, we are currently improving the framework so that it can extract genomic structural variations
causing splicing changes and so on.


Classification of splicing changes
----------------------------------
We consider *Exon skipping*, *Alternative 5' splice site (SS)*, *Alternative 3'SS* and *Intron retention*
for the splicing aberrations caused by genomic variants.
The first three types of aberrations are collected by split-alignment (therefore considered as introns).
Using control panels (set of non-disease transcriptome data) is highly recommended
to remove *normal* splicing aberrations to improve the accuracy.

:*Exon skipping*:
  Both the edges of introns are annotated exon-intron junctions.
  However, there is no annotated intron supporting in the transcriptome database (we use RefSeq).

:*Alternative 5'SS*:
  One of the edge of introns corresponds to the annotated splicing acceptor site.
  When unannotated edge is located in intronic region (>= 30bp from the annotated splicing donor site),
  we identify them as Intronic alternative 5'SS.

:*Alternative 3'SS*:
  One of the edge of introns corresponds to the annotated splicing donor site.

For identifying intron retention,
we count the number of short reads covering both exon and intron regions
(>= 8bp each, in the default setting)
at each annotated splicing donor/acceptor site.


Classification of SAVs
----------------------

Here, we use the term *splicing donor sites* and *splicing acceptor sites* as regions
around initial and final edges of introns, respectively.
More specifically, in the default setting of SAVNet,
*splicing donor sites* include the last 3bp in the exons and the first 6bp in the introns,
whereas *splicing acceptor* sites cover the last 6bp in the introns and the first 1bp in the exons.
These values are determined by evaluating false discovery ratios of SAVNet
at each relative position in splicing sites when performing SAVNet.
Also, we call the first and the last 2bp of introns (GT-AG) as *canonical splice sites*.


:Splicing donor/acceptor disruption:

  When the genomic variants are affecting the annotated splicing donor/acceptor sites
  and splicing aberrations are overlapping this somatic mutations
  (allowing some amount of margins, 50bp by default),
  we nominate these relationships as *splice site disruption*.

:Splicing donor/acceptor creation:

  When unannotated edges of Alternative 5' or 3'SS are overlapping with the genomic variants,
  we nominate these relationships as *splicing site creations*.
