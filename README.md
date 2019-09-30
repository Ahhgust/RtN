# Removing Numts


[Numts](https://en.wikipedia.org/wiki/NUMT) (pronounced "new mights") are mitochondrial DNA segments that have been 
transposed into the nuclear genome. If you work with mitochondrial genome data 
they're a problem as there's no *good* way to identify them. <br><br>

Numts have different properties than their mitochondrial counterparts; the nuclear genome has a different (generally lower) rate of molecular evolution 
and a different mode of inheritance (not just matrilineal). In general if you know that your data are single 
source then Numts don't tend to be a problem; most of the time the majorty basecall is correct.
However, if you're trying to call [heteroplasmies](https://en.wikipedia.org/wiki/Heteroplasmy) or you 
think you may have a mixture of DNAs it would be good to remove such sequences from your analysis. <br><br>

One way to remove Numts is to use an alignment-based approach-- map your reads to, say, the GRCh38 reference genome, as it contains (some) Numts, reads that look more similar to the Numts in this reference will be mapped there instead of to the mitochondrial genome. This has two major short-comings. GRCh38 only has some Numts; not all. In fact, some Numt insertions are polymorphic, which means that some individuals will have them, others will not, and further, that the set of all Numts will always be incomplete (as you haven't sampled all individuals in the population). The second shortcoming to reference-based approaches that when you map reads you're mapping an allele, an allele that may be fairly different from the reference genome(s). E.g., what if your read has 4 differences to the rCRS, and 3 differences to some Numt? What does that mean? And does the answer change if we know that the read is an exact hit to some known mitochondrial sequence, just not the rCRS (hint, it does)? <br><br>

The premise behind RTN is simple: only keep reads that map *well* to some known mitochondrial sequence. It uses annotated genomes from [HmtDB](https://www.hmtdb.uniba.it/). If a read is similar to some known sequence it is kept, otherwise it's mapping quality is set to 0. RTN also maps reads to a database of annotated Numt alleles. The database includes alleles from: [Dayama et al](https://doi.org/10.1093/nar/gku1038), [Calabrese et al](https://doi.org/10.1186/1471-2105-13-S4-S15) and [Smart et al](https://doi.org/10.1016/j.fsigen.2019.102146). The reads are then decorated with two tags: ZH, which gives the minimum distance of the read to some annotated human sequence, and ZN, which gives the same distance but to the Numt alleles. The minimum distance is defined in two ways: either ignoring indels (just using the matched bases in the CIGAR string), or it can include the number of indels as well-- I would guess that the former is better for Ion sequencing, while the latter would be better for Illumina.<br><br>

## Quick start





