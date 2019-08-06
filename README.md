# Removing Numts


Numts (pronounced "new mights") are mitochondrial DNA segments that have been <br>
transposed into the nuclear genome. If you work with mitochondrial genome data <br>
they're a problem as there's no *good* way to identify them. <br><br><br>

Numts have different properties than their mitochondrial counterparts; they have different rate of molecular evolution <br>
and a different mode of inheritance. In general if you know that your data are single <br>
source then Numts don't tend to be a problem; most of the time the majorty basecall is correct.<br> However, if you're trying to call heteroplasmies or think you may have <br> a mixture of DNAs it would be good to remove such sequences from your analysis. <br><br>

One way to remove Numts is to use an alignment-based approach-- map your reads to, say, the GRCh38 reference genome, as it contains (some) Numts, reads that look more similar to the Numts in this reference will be mapped there instead of to the mitochondrial genome. <br>




