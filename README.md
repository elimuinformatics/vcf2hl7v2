## **VCF to HL7 V2 Converter**

### Introduction

VCF-formatted files are the lingua franca of next-generation sequencing, whereas HL7 Version 2 messaging format is the predominant means by which labs send structured results to Electronic Health Records (EHRs). EHRs can directly import HL7 V2 formatted results that conform to the [HL7 Lab Results Interface Implementation Guide](https://www.hl7.org/documentcenter/public/standards/dstu/V251_IG_LRI_R1_STU3_2018JUN.pdf). Here, we provide an open source utility for converting variants from VCF format into HL7 V2 format. Details of the translation logic are on the [manual page](https://github.com/elimuinformatics/vcf2hl7v2/blob/master/docs/Manual.md).

### Install
Before installing vcf2hl7v2 you need to install cython and wheel.
```
pip install cython wheel
```
Now, install vcf2hl7v2 binary from pip.
```
pip install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple/ vcf2hl7v2
```

### Quick Examples
(some sample VCF files are [here](https://github.com/elimuinformatics/vcf2hl7v2/tree/master/vcf2hl7v2/test))

```python
>>> import vcf2hl7v2
>>> vcf_hl7v2_converter = vcf2hl7v2.Converter('sample.vcf', 'GRCh37')
>>> vcf_hl7v2_converter.convert()
```

### Logging

You can use python standard [logging](https://docs.python.org/3/library/logging.html) to enable logs. Two loggers ('vcf2hl7v2.general') and ('vcf2hl7v2.invalidrecord') are available to configure.
* **vcf2hl7v2.general**: standard library logs. 
* **vcf2hl7v2.invalidrecord**: logs all the records from vcf file which are in conversion region but are not converted to HL7 V2 format.

```python
>>> import logging
# create logger
>>> logger = logging.getLogger('vcf2hl7v2.invalidrecord')
>>> logger.setLevel(logging.DEBUG)
# create console handler and set level to debug
>>> ch = logging.FileHandler('invalidrecord.log')
>>> ch.setLevel(logging.DEBUG)
# create formatter
>>> formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
# add formatter to ch
>>> ch.setFormatter(formatter)
# add ch to logger
>>> logger.addHandler(ch)
```

### Scope

Software converts somatic and germline simple variants (SNVs, MNVs, Indels), along with zygosity and phase relationships, for autosomes, sex chromosomes, and mitochondrial DNA. Select clinical annotations (e.g. clinical significance, phenotype), supplied in a separate annotation file are incorporated into the conversion.

* Not supported
    * **Structural variants**: Software does not support conversion of structural variants (where INFO.SVTYPE is present). 
    * **Alt contigs**: Software does not support conversion of variants aligned to Alt contigs. We recommend caution in using this software against VCFs generated with an alternate-locus aware variant caller, as variants mapped to Alt contigs will not be converted.
    * **Query liftover**: Software assumes that regions (conversion region, studied region) and VCF are based on the same genomic build. 
    * **Chromosome synonyms (e.g. '1' vs. 'chr1')**: Software assumes that chromosome representation is consistent between regions (e.g. in BED files) and VCF. For instance, if VCF uses 'chr1', then BED file must also use 'chr1' 

### License and Limitations

Software is available for use under an [Apache 2.0 license](https://opensource.org/licenses/Apache-2.0), and is intended solely for experimental use, to help further Genomics-EHR integration exploration. Software is expressly not ready to be used with identifiable patient data or in delivering care to patients. Code issues should be tracked here. Comments and questions can also be directed to info@elimu.io or srikarchamala@gmail.com.
