# Remove the Numts!


[Numts](https://en.wikipedia.org/wiki/NUMT) (pronounced "new mights") are mitochondrial DNA segments that have been 
transposed into the nuclear genome. If you work with mitochondrial genome data 
they're a problem as there's no *good* way to identify them. <br><br>

Numts have different properties than their mitochondrial counterparts; the nuclear genome has a different (generally lower) rate of molecular evolution 
and a different mode of inheritance (not just matrilineal). In general if you know that your data are single 
source then Numts don't tend to be a problem; most of the time the majorty basecall is correct.
However, if you're trying to call [heteroplasmies](https://en.wikipedia.org/wiki/Heteroplasmy) or you 
think you may have a mixture of DNAs it would be good to remove such sequences from your analysis. <br><br>

One way to remove Numts is to use an alignment-based approach-- map your reads to, say, the GRCh38 reference genome, as it contains (some) Numts, reads that look more similar to the Numts in this reference will be mapped there instead of to the mitochondrial genome. This has two major short-comings. GRCh38 only has some Numts; not all. In fact, some Numt insertions are polymorphic, which means that some individuals will have them, others will not, and further, that the set of all Numts will always be incomplete (as you haven't sampled all individuals in the population). The second shortcoming to reference-based approaches that when you map reads you're mapping an allele, an allele that may be fairly different from the reference genome(s). E.g., what if your read has 4 differences to the rCRS, and 3 differences to some Numt? What does that mean? And does the answer change if we know that the read is an exact hit to some known mitochondrial sequence, just not the rCRS (hint, it does)? <br><br>

The premise behind RtN! is simple: only keep reads that map *well* to some known mitochondrial sequence. The simplicitly of this approach is that Numts don't have to be annotated; if a read is very far from all known mitochondrial sequences it very likely harbors a lot of sequencing errors (better to ignore it), or it's an off-target alignment (which we don't want). RTN uses annotated genomes from [HmtDB](https://www.hmtdb.uniba.it/). If a read is similar to some known sequence it is kept, otherwise it's mapping quality is set to 0. RTN also maps reads to a database of annotated Numt alleles. The database includes alleles from: [Dayama et al](https://doi.org/10.1093/nar/gku1038), [Calabrese et al](https://doi.org/10.1186/1471-2105-13-S4-S15) and [Smart et al](https://doi.org/10.1016/j.fsigen.2019.102146). The reads are then decorated with two tags: Z**H**, which gives the minimum distance of the read to some annotated **H**uman sequence, and Z**N**, which gives the same distance but to the **N**umt alleles. The minimum distance is defined in two ways: either ignoring indels (just using the matched bases in the CIGAR string), or it can include the number of indels as well-- I would guess that the former is better for Ion sequencing, while the latter would be better for Illumina.<br><br>

As we describe in our paper, most reads are ambiguous and exactly match both a human sequence and a Numt sequence.
The principle of RtN! is to keep ambiguous reads and flag them for manual review. 
There are processing steps that can be done to make RtN! more aggressive in these cases, 
but to be clear, some types of Numts are intractable. 

<br><br>
## Quick start

Clone this repo (NO recursive downloads SeqLib, htslib, bwa, fermi-kit, ... )
```
git clone https://github.com/Ahhgust/RtN.git
cd RtN
```

Change directories.
Uncompress and index *humans.fa.bz2* a la:

```
bunzip2 humans.fa.bz2 && bwa index humans.fa
```
(This will take a while)
<br>
The Numt databases need no decompression/indexing
<br>
Select the appropriate static binary <br>
[*nix](https://github.com/Ahhgust/RtN/tree/master/Nix_binary)
<br> or <br>
[Windows](https://github.com/Ahhgust/RtN/tree/master/WSL_binary)
<br>
And run it (help):
```
./rtn -h
```
<br>

Recommendations for single-end reads (likely Ion sequencing)

```
./rtn -h humans.fa -n Calabrese_Dayama_Smart_Numts.fa -b yourBam.bam
```

In words:<br>
-h: Compute near-neighbor distances to the human sequences (and remove reads that align poorly to all known human haplotypes)<br>
-n: And compute the near-neighbor distance to an (always incomplete) database of Numt alleles<br>
-b: And do so using the bam files specified...<br>


<br>

Recommendations for paired-end reads (likely Illumina sequencing)

```
./rtn -p  -h humans.fa -n Calabrese_Dayama_Smart_Numts.fa -b yourBam.bam
```

<br>
-p: The same as above, but run in paired-end mode<br>


## Single-source samples
If you *know* that your data are single source there's an additional trick that can be used.
Rather than aligning the reads to HmtDB it suffices to align the reads to just the single source haplotype.
We call this running a bootstrap (as in the "pull yourself up by" as opposed to the resampling technique).
To do a bootstrap, you can either use a fasta sequence (for the given individual) and then index it with bwa
(aka,
```bwa index yourFasta.fa```
)
<br>
As an added convenience we've added a simply python script (empop2fa.py)
that reads in an Empop-style encoded haplotype, and produces an indexed fasta.

Use it as:
```
fgrep -v '#' empopFile | python3 empop2fa.py rCRS.fa
```

By Empop-style file we mean a file that describes sequence differences to the rCRS sequence.
The file is presumed to be tab-separated, with the variants starting at column 4 onwards.
Click [here](https://raw.githubusercontent.com/gmitirol/empophub/master/AUT273_spec.emp) for an example.
<br>

RtN! can then be run similar to the above way.
In practice we have had excellent luck using a more stringent likelihood threshold in this case (1e-3),
```
./rtn -L 1e-3 -p -h mySingleSource.fa -n Calabrese_Dayama_Smart_Numts.fa -b yourBam.bam
```




# How to compile!
 (somewhat unusual) Requirements:
 c++ json library (install as necessary)
 
```
sudo apt-get install libjsoncpp-dev
```

## Compilation instructions (with dynamic libraries; necessary for curl support)
```
git clone --recursive https://github.com/Ahhgust/RtN.git
cd RtN/SeqLib/htslib # Make htslib
autoconf
autoheader
./configure --enable-libcurl
make
cd ../bwa        # Make BWA
make
cd ../fermi-lite # and Fermi-lite
make
cd ..            # and now the SeqLib library
./configure LDFLAGS='-lcurl -lcrypto'
make
cd ..        
make             # and make RtN
```


## Compilation instructions (with static libraries)
```
git clone --recursive https://github.com/Ahhgust/RtN.git
cd RtN/SeqLib/htslib # Make htslib
autoconf
autoheader
./configure --disable-libcurl
make
cd ../bwa        # Make BWA
make
cd ../fermi-lite # and Fermi-lite
make
cd ..            # and now the SeqLib library
./configure
make
cd ..        
make static             # and make RtN
```


