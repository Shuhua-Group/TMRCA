> **TMRCA.py**
>
> https://github.com/Shuhua-Group/TMRCA
>
> For calculation of within-population TMRCA & cross-population divergence time
>
> developed in python2.7, at your own risk if you use python3
>
> require phased VCF as input, accept missing geno (i.e., ".|.")
>
> 
>
> If you use the **TMRCA.py** in your research work, please cite at least one of the following paper(s):
>
> - 



**Simple run: **

``` bash
## to calculate TMRCA using all samples in the input file at the given region(s)
$ python2 TMRCA.py --gzvcf input.vcf.gz --region region.txt --ape Human_panTro5.chr@.diff.txt.gz

## to get help
$ python2 TMRCA.py -h
```

**Inputs:**

> **`required`**
>
> **--gzvcf**: phased VCF.gz file, in GT format, i.e., "1|0". accept missing geno (.|.). male chrX should be homozygote.
>
> **--region**: bed file, one region indicated in each line. 3 columns: `chr` `start` `end` . no header line, tab or space delimited, additional columns will be ignored. chromosome ID(s) should be coded in the same way as those in the VCF file (i.e., "1" is different from "chr1")
>
> **--ape**: difference between human and chimpanzee (or any other outgroup) genome. each line indicates a genetic position with different alleles between human genome and chimpanzee genome. 2 columns: `chr` `pos` , additional columns will be ignored, no header, tab or space delimited. one file for one chromosome. chromosome ID should be replaced by "@" in the command line (i.e., --ape chr@.txt). get all the files prepared under same folder.
>
> **`optional`**
>
> **--samples**: file for sample info. each line indicates the sample ID, population ID, and also which haplotype to be used (or both haplotypes). 2 columns: `sampleID_1/2` `popID` . no header line. For example, "sampleX_2 popX" stands for the 2nd haplotype from sampleX, and this sample belongs to popX. **default**: "all" samples in the input vcf data will be used and considered as one pop. 
>
> **--Tind**: `T` / `F` . whether to estimate pairwise TMRCA between individuals/haplotypes. **default**: F. it may take some additional time if "T". 
>
> **--pairs**: indicate pairs of populations between which the time should be estimated. one pair of populations in each line. 2 columns: `popN` `popM` . pop name(s) should be in the **--samples** file. no header line, tab or space delimited, additional columns will be ignored. **Default**: all possible pairs of populations will be considered. 
>
> **--divT**: divergence time between humna and chimpanzee (or any other outgroup) in years. **Default**: 13e6
>
> **--out**: output file prefix. **default**: out . 



---

By: Yuwen Pan, 2022  |  Contact: panyuwen.x@gmail.com