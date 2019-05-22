# fastKaptive

## Introduction

`fastKaptive` is a Blast-based R written software to search and extract capsule synthesis loci in bacterial genomes. It has been used and tested on genomes of the Enterobacteriales, as described in [1]. The code was designed to extract to locus regions, based on a model of *Escherichia coli* [2]:

* *cps* locus region, associated with the synthesis of O-antigens in *E. coli*,
* *kps* locus region, associated with the synthesis of K-antigens in *E. coli*.

The approach is based on a previously published and tested method, `Kaptive` [3], but has been reoptimised for scanning large sets of genomes by avoiding `tblastx` and using `blastn` instead. 

## Dependencies
The software uses `blastn` and `gzip` as dependencies, and it order to work they should be installed prior to using `fastKaptive` and the commands should be visible in `PATH`. Also, these R packages need to be installed beforehand via the `install.packages()` command:
* `seqinr`
* `igraph `
* `optparse `

## How to use it?
The code was developed on Mac OS X and should work on Linux-based platforms via a command line. The programme should be launched from a working directory by revoking:
```
Rscript --vanilla ${fK_PATH}/main.R (options)
```
where `fK_PATH` is the location of the `fastKaptive` directory.

Using the `-h` option, one can display the options:
```
Options:
	-l LOCUS.NAME, --locus.name=LOCUS.NAME
		Which locus to search: cps [default] or kps

	-e, --extended.db
		Use an extended locus search database  [default FALSE]

	-D, --distant.search
		Search in distant bacteria using blastn instead of megablast [default FALSE]

	-p PATH.FILE, --path.file=PATH.FILE
		File with path names to all assemblies [default paths.input.txt]

	-o OUTPUT.FOLDER, --output.folder=OUTPUT.FOLDER
		Name of the output folder [default extract.output]

	-h, --help
		Show this help message and exit
```

The first option `-l` specifies the locus which should be searched, either *cps* or *kps*. Any other option will return an error.

The second option `-e` specifies whether an extended database should be used. By default the code uses a high-quality locus references from *Klebsiella*, *Escherichia/Shigella* and *Salmonella* for *cps* and only *Escherichia* for *kps*. However, an extended reference database, obtained by searching 27,365 genomes of Enterobacteriales in [1], and extracting good quality loci from draft genomes, amounting to 994 references for *cps* and 106 references for *kps*.

The third option `-D` specifies whether `blastn` option should be used with `blast` instead of `megablast`: if used, a more time-consuming search will ensue, resulting in more accurate searches in other species of Enterobacteriales. It is recommended to use when searching outside of the Enterobacteriaceae family.

The fourth option `-p` specifies the location of a file with assembly paths. If these paths are incorrect, the code will return an error. **This is a mandatory option.**

The fifth option `-o` specifies the name of the output folder.

## References

1. Holt KE, Lassalle F, Wyres KL, Wick R & Mostowy RJ, *Diversity and evolution of surface polysaccharide synthesis loci in Enterobacteriales*, submitted
2. Whitfield C *Biosynthesis and assembly of capsular polysaccharides in Escherichia coli*, Annu Rev Biochem. 2006;75(1):39-68. [doi:10.1146/annurev.biochem.75.103004.142545](https://doi.org/10.1146/annurev.biochem.75.103004.142545).
3. Wyres KL, Wick RR, Gorrie C, Jenney A, Follador R, Thomson NR & Holt KE, *Identification of Klebsiella capsule synthesis loci from whole genome data*. Microbial Genomics. 2016;2(12):e000102. [doi:10.1099/mgen.0.000102](https://doi.org/10.1099/mgen.0.000102).


