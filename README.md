# segcpp - a standalone program to run segregation analysis

Segregation analysis studies the distribution of genetic variants within and between families of related individuals, in order to make inferences about the heritability of variation

segcpp inputs one or more uncompressed vcf files and a single pedigree file, and outputs a tab-separated file with all the information needed to analyze variant segregation in families

## command-line arguments
```
Input/output parameters:
  -v [ --vcf-file ] arg                 vcf input file(s) and optional filters. 
                                        Format: file.vcf[:pR,pX,p%,pQ:fcR,fcX,fc%,fcQ:nfcR,nfcX,nfc%,ncfQ]
  -p [ --ped-file ] arg                 pedigree input file
  -o [ --out-file ] arg                 output file
  -i [ --add-info ] arg                 supplementary INFO field(s) to output
  -f [ --add-format ] arg               supplementary FORMAT field(s) to output

Treshold values for variant filtration:
  --proband-var-reads arg (=1)          min. variant reads in probands
  --proband-coverage arg (=1)           min. coverage in probands
  --proband-mutfreq arg (=0)            min. mutation frequency in probands
  --proband-genotype-quality arg (=0)   min. genotype quality in probands
  --familial-control-var-reads arg (=1) min. variant reads in familial controls
  --familial-control-coverage arg (=1)  min. coverage in familial controls
  --familial-control-mutfreq arg (=0)   min. mutation frequency in familial controls
  --familial-control-genotype-quality arg (=0)
                                        min. genotype quality in familial controls
  --external-control-var-reads arg (=3) min. variant reads in external controls
  --external-control-coverage arg (=6)  min. coverage in external controls
  --external-control-mutfreq arg (=0.2) min. mutation frequency in external controls
  --external-control-genotype-quality arg (=20)
                                        min. genotype quality in external controls
  --show-probands-with-variants arg (=10)
                                        max. proband names to show
  --show-familial-controls-with-variants arg (=10)
                                        max. familial control names to show
  --show-external-controls-with-variants arg (=10)
                                        max. external control names to show
  --use-genotyper-reads                 use AD reads (instead of ADS) for coverage filtering
  --use-indel-filtering                 use filters for indel variants
  --use-format-low-stringency           use a low format stringency for parsing (implies use-genotyper-reads)
  --show-all-variants                   show wildtype calls in probands and variants found in controls

General options:
  -h [ --help ]                         produce an help message and exit
  --verbose [=arg(=1)] (=0)             ouput more details [0,1,2]
  --version                             print program version and details and exit
```

## general concepts

#### family-variant:  
The basic unit of information in the program output is a family-variant, i.e. a single variant in a single family. A variant is defined by a unique combination of a single chromosome, position, reference allele and alternate allele. A family is defined by the first column of the input pedigree file. Each row of the program outputs corresponds to a single family-variant. This means that if a variant is present in >1 family (as defined by the pedigree file), then that variant will appear on a separate row for each of the families that contain this variant.

#### filtering:  
segcpp allows for a measure of QC by giving users the option of filtering individual sample genotypes based on a few simple sequencing metrics:  
1. variant reads (i.e. reads with ALT allele, integer)
2. coverage (i.e. total reads, integer)
3. variant read frequency (i.e. fraction of reads with ALT allele, value from 0-1)
4. genotype quality  (i.e. per-sample vcf FORMAT GQ field, integer)

Filtering a sample genotype does *not* mean that it will be removed from the final outputs. Instead, if a given sample genotype fails to meet the appropriate filters for a given variant, its genotype will appear with the prefix "Filtered". So for example a heterzygous genotype which fails to meet the filtering criteria will appear as "Filtered Heterozygote". Note that "Filtered" genotypes are counted separately in the variant counting columns of the program's output (see [below](#sample-counts-and-ids))

## file inputs
* vcf file: [Standard VCF format](https://samtools.github.io/hts-specs/VCFv4.2.pdf) containing list of variants, sample genotypes and associated annotations. Can be specified more than once. File-specific filters can also be specified here (see [below](#genotype-filtering-arguments))
* ped file: [Standard pedigree file format](https://gatk.broadinstitute.org/hc/en-us/articles/360035531972-PED-Pedigree-format), with no header

## other arguments

#### genotype filtering arguments
Filters must be set separately for affected samples ("probands"), in-family non-affecteds ("familial controls") and non-family non-affecteds ("external controls"). For more fine-grained control each of these values can be set separately for each input vcf file, using the --vcf-file:[additional filters] notation as in the [above section](#command-line-arguments). If doing so, the additional filters *must* be supplied in the order [above](#filtering), in the form
`--vcf-file /path/to/file:proband variant reads,proband coverage,proband variant read frequency,proband GQ:familial control variant reads,FC coverage, FC variant read frequency,FC GQ:external control variant reads,EC coverage,EC variant read frequency,EC GQ`  

By default, filtering is not applied to indel (insertion/deletion) variants. This behaviour can be changed with the `--use-indel-filtering` command switch

#### display arguments
In addition to providing a detailed breakdown of variant counts for each family-variant, segcpp displays a list of samples with each variant. To avoid file bloating, this list will be replaced by an ellipsis ("...") in the event that there are >10 samples (by default) with the variant. These lists are displayed separately for probands, familial controls and external controls. Default values can be changed with `--show-probands-with-variants`, `--show-familial-controls-with-variants` and `--show-external-controls-with-variants` respectively.

#### parsing arguments
By default, segcpp expects per-allele per-stand read depths to be specified in the vcf FORMAT ADS field, i.e. forward REF,reverse REF,forward ALT,reverse ALT. The command-line switch `--use-genotyper-reads` will instead force it to look for the more standard AD FORMAT field, which does not distinguish by strand. Later versions of segcpp will include this latter behaviour as default.

More generally, to prevent fatal file format-based warnings, use the `--use-format-low-stringency` command switch

#### other arguments
By default, segcpp will only output a family-variant line if that variant appears in at least one proband in that family. Use the `--show-all-variants` if you wish instead to display *all* variants which appear in *any* family sample regardless of clinical status. NOTE: if performing a [case-control study](#case-control-studies), use this option.


## understanding outputs

The output file columns fall into several unofficial categories, which appear in the following order from left to right:  

#### family-variant definition
Basic information about the family-variant: familyID, chromosome, position, REF and ALT alleles

#### variant-level annotations
All INFO field annotations from the input vcf file(s), specified for inclusion by the -i|--add-info program argument. Column headers are parsed directly from the INFO field definitions in the input vcf file headers.

#### sample counts and ids
These columns contain the information needed for segregation analysis itself. They break down the counts of variant presence along several axes:  
* clinical status (affected / non-affected / all)  
* family membership (family / external / all)  
* variant status (variant / wildtype / filtered)  

Each unique combination of axes has its own column in this section: e.g. "Family affected w/ variant", "Total subjects wildtype", "External controls filtered", etc.

#### sample-specific information
The last set of columns provide sample-specific information extracted from the vcf FORMAT fields. The first column in this section, "Family members", lists all samples in the current family. The next six columns (by default) provide the following information for the *first* sample listed in the "Family members" column:
* "Family member zygosity": the specific genotype for the sample. Can be "Heterozygous", "Homozygous", "Wildtype" or any of the previous with the prefix "Filtered"
* "Family member genotype qual": the sample's vcf FORMAT GQ value
* "Family member genotyper ref cov": the number of reads with the REF allele
* "Family member genotyper var cov": the number of reads with the ALT allele
* "Family member genotyper cov": the total number of reads
* "Family member genotyper mut freq": the fraction of reads with the ALT allele

Subsequent sets of six columns present the same six pieces of information for each family sample, in the order listed in the "Family members" column. If a sample's six columns are blank this means that there was no call for that variant in that sample. Note that there are no column headers for all samples after the first one.

## software requirements

As segcpp is written in C++, it must be compiled locally. This requires a compiler such as GNU make, as well as the boost portable C++ libraries.  

Additionally, the helper script scripts/format_vep.sh requires a local version of bcftools, with the split-vep plugin installed, and with the bcftools executable in the user's $PATH environment variable

## system requirements

segcpp does not support multithreading, so it only requires a single CPU to run.

However, memory requirements scale with the number of variants in probands (or in all samples if running in `--show-all-variants` mode), as well as the total number of annotations to be included. segcpp has not been extensively benchmarked; however we have managed to run anlayses of hundreds of whole-genome samples on past versions of this software using a single node with 96GB RAM in under 24 hours.

Should the analysis require more resources than are available, we recommend splitting vcf inputs and running segcpp on each piece.

## miscellaneous
#### case-control studies
While segcpp was designed with familial analyses in mind, it can easily accomodate case-control studies as well. To do so, simply set up your pedigree file such that all samples have the same familyId, with cases having the phenotype "affected" (2) and controls having the phenotype "unaffected" (1). Run segcpp with the `--show-all-variants` switch. The counting columns should provide the necessary raw inputs for subsequent statistical analyses.
