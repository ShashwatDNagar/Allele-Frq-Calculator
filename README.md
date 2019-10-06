# Allele Frequency Calculator

Script for calculating allele frequencies for the alternate allele in VCF files *given different grouping schemes*.

## Requirements
Tested on *Python 3.7.3*

## Usage

```
usage: get_vcf_frqs.py [-h] --vcf VCF --groups GROUPS --output OUTPUT
                       [--use-fid]

A population-specific SNP frequency calculator that let's you specify
different population groups.

required arguments:
  --vcf VCF        Input VCF file - may be gzipped.
  --groups GROUPS  Grouping file - tab-separated file that defines group
                   membership. Each row is one individual followed by
                   group(s). See the Github () for a sample file.
  --output OUTPUT

optional arguments:
  -h, --help       show this help message and exit
  --use-fid        Use Family ID (the second to last word after splitting
                   individual identifiers on '_'. Defaults to using the
                   individual ID (last word after split). This is the ID that
                   is used for group mapping - choose wisely.
```

## Input files

1.  VCF files (compressed or uncompressed)
2.  Grouping files
	```
	sampleID	groupingScheme1	groupingScheme2
	1	Group1	Category1
	2	Group3	Category2
	3	Group2	Category3
	4	Group2	Category3
	5	Group2	Category3
	```
	Header *mandatory*.

# Output format

```
Chr	Position	rsID	Ref	Alt	<grouping scheme 1, category 1>_freq	<grouping scheme 1, category 1>_count	<grouping scheme 2, category 1>_freq	<grouping scheme 2, category 1>_count ...
```

