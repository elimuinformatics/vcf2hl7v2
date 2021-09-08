# v0.0.1 (TBD)

Converts VCF variants into HL7 V2 format that conform to the [HL7 Lab Results Interface Implementation Guide](https://www.hl7.org/documentcenter/public/standards/dstu/V251_IG_LRI_R1_STU3_2018JUN.pdf). 

## Added
* In scope are somatic and germline simple variants (SNVs, MNVs, Indels), along with zygosity and phase relationships, for autosomes, sex chromosomes, and mitochondrial DNA.
* Supports VCF file (text-based or bgzipped) and optionally tabix files for query.
* Supports genome build ('GRCh37' or 'GRCh38');
* Optionally supports a query region in the form of .bed or dictionary that indicates the region(s) to convert.
* Optionally supports inclusion of  'region-studied' observations that detail which portions of the conversion region were studied.
* Optionally supports clinical annotations (e.g. clinical significance, phenotype), supplied in a separate annotation file.
