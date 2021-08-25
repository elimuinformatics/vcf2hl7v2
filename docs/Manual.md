
vcf2hl7v2 Manual
===================

Introduction
-------------------------


Conceptually, the utility takes a VCF as input and outputs a set of HL7 V2 OBX observations. These OBX observations can be then incorporated into an overarching HL7 V2 ORU message.

Conversion logic
================

Conversion region
-----------------

Variants in the VCF file are intersected against an optional file listing variants to convert and a file listing clinical annotations. The following table summarizes the scope of VCF records converted based on these regions.

| Conversion Region      | Studied Region | Annotations     |  Output     |
| :---:        |    :----:   |          :---: |  :---       |
| Not Supplied      | Not Supplied       | Not Supplied   | <ul><li>Convert all variants in VCF.</li><li>HL7 V2 message contains no region-studied OBX observation group.</li></ul>            |
| Not Supplied      | Not Supplied       | Supplied   | <ul><li>Convert all variants in VCF for which annotations are provided.</li><li>HL7 V2 message contains no region-studied OBX observation group.</li></ul>            |
| Not Supplied      | Supplied       | Not Supplied   | <ul><li>Convert all variants in VCF.</li><li>HL7 V2 message contains one region-studied OBX observation group per studied chromosome.</li><ul><li>Studied region(s) reflected in ranges-examined component(s).</li></ul></ul>            |
| Not Supplied      | Supplied       |  Supplied   | <ul><li>Convert all variants in VCF for which annotations are provided.</li><li>HL7 V2 message contains one region-studied OBX observation group per studied chromosome.</li><ul><li>Studied region(s) reflected in ranges-examined component(s).</li></ul></ul>            |
|  Supplied      | Not Supplied       |  Not Supplied   | <ul><li>Convert all variants in conversion region.</li><li>HL7 V2 message contains no region-studied OBX observation group.</li></ul>            |
|  Supplied      | Not Supplied       |   Supplied   | <ul><li>Convert all variants in conversion region for which annotations are provided.</li><li>HL7 V2 message contains no region-studied OBX observation group.</li></ul>           |
|  Supplied      | Supplied       |  Not Supplied   | <ul><li>Convert all variants in conversion region.</li><li>HL7 V2 message contains one region-studied OBX observation group per studied chromosome intersected with conversion region.</li><ul><li>Studied region(s), intersected with conversion region, reflected in ranges-examined component(s).</li></ul></ul>
|  Supplied      | Supplied       |  Supplied   | <ul><li>Convert all variants in conversion region for which annotations are provided.</li><li>HL7 V2 message contains one region-studied OBX observation group per studied chromosome intersected with conversion region.</li><ul><li>Studied region(s), intersected with conversion region, reflected in ranges-examined component(s).</li></ul></ul>  


General conversion
------------------

### Exclude VCF rows

The following VCF rows are excluded from conversion:

-   VCF REF is not a simple character string

-   VCF ALT is not a simple character string, comma-separated character string, or \'.\'.

-   VCF FILTER does not equal \'PASS\' or \'.\'.

-   VCF INFO.SVTYPE is present. (Structural variants are excluded).

-   VCF FORMAT.GT is null (\'./.\', \'.\|.\', \'.\', etc).

### Region-studied observations

<ul>
    <li>For each region studied, create a group of OBXs, grouped by SubID (1st group is 1a, second group is 1b, etc). (See HL7 Lab Results Interface Implementation Guide Section 5.4 for more details on nesting representation).

<li>Within a group, include these OBXs
    <ul>
    <li>[1..1] 48013-7^Genomic reference sequence^LN (SubID=1a)
    <li>[1..*] 51959-5^Range(s) of DNA sequence examined^LN (SubID=1a.a, 1a.b, 1a.c, …. 1a.aa..1a.az, 1a.ba..1a.bz, etc) (assumes 1-based)
    </ul>

<li>After groups, include this OBX:
    <ul>
    <li>[1..1] 51968-6^Discrete variation analysis overall interpretation^LN (SubID=1). If any variants, then OBX-5 is 'Positive'. Otherwise, OBX-5 is 'Negative'.
    </ul>
</ul>

**Example**:

**first group (subID=1a)**

```
OBX|1|ST|48013-7^Genomic reference sequence^LN|1a|NC_000005.9|
OBX|2|NR|51959-5^Range(s) of DNA sequence examined^LN|1a.a|112043201^112181936|
```

**second group (subID=1b)**

```
OBX|3|ST|48013-7^Genomic reference sequence^LN|1b|NC_000017.10|
OBX|4|NR|51959-5^Range(s) of DNA sequence examined^LN|1b.a|41196311^41277500|
OBX|5|NR|51959-5^Range(s) of DNA sequence examined^LN|1b.b|41550000^41600000|
```

**third group (subID=1c)**

```
OBX|6|ST|48013-7^Genomic reference sequence^LN|1c|NC_000019.9|
OBX|7|NR|51959-5^Range(s) of DNA sequence examined^LN|1c.a|38924339^39078204|
```

```
OBX|8|TX|51968-6^Discrete variation analysis overall interpretation^LN|1|Positive|
```

### Variant observations

<ul>
<li>For each variant, create a group of OBXs, grouped by SubID (1st group is 2a, second group is 2b, etc)</li>
<li>Within a group, include these OBXs:</li>
    <ul>
    <li>[1..1] 83005-9^Variant Category^LN (SubID=2a): 'Simple'</li>
    <li>[1..1] 47998-0^Variant Display Name^LN (SubID=2a): populate with contextual SPDI, built from refSeq (48013-7), start (81254-5), ref allele (69547-8), alt allele (69551-0) as refSeq:start-1:ref:alt</li>
    <li>[1..1] 48018-6^Gene Studied^LN (SubID=2a): Populate with gene from annotation file. If no gene is present, set equal to ''HGNC:0000^NoGene^HGNC".</li>
    <li>[1..1] 48004-6^DNA Change c.HGVS^LN (SubID=2a):</li>
        <ul>
        <li>If both cHGVS and transcriptRefSeq: transcriptRefSeq+':'+cHGVS</li>
        <li>elseIf cHGVS but not transcriptRefSeq: cHGVS</li>
        <li>else: populate with 47998-0 value</li>
        </ul>
    <li>[0..1] 48005-3^Amino Acid Change p.HGVS^LN (SubID=2a): if pHGVS is supplied, populate with transcriptRefSeq+':'+pHGVS</li>
        <ul>
        <li>If both pHGVS and proteinRefSeq: proteinRefSeq+':'+pHGVS</li>
        <li>elseIf pHGVS but not proteinRefSeq: pHGVS</li>
        <li>else: omit OBX</li>
        </ul>
    <li>[1..1] 48013-7^Genomic reference sequence^LN (SubID=2a)</li>
    <li>[1..1] 69547-8^Genomic ref allele^LN (SubID=2a)</li>
    <li>[1..1] 81254-5^Genomic allele start-end^LN (SubID=2a)</li>
    <li>[1..1] 69551-0^Genomic alt allele^LN (SubID=2a)</li>
    <li>[1..1] 48002-0^Genomic Source Class [Type]^LN (SubID=2a):</li>
        <ul>
        <li>if sourceClass = 'germline': LA6683-2^Germline^LN</li>
        <li>elseIf sourceClass = 'somatic': LA6684-0^Somatic^LN</li>
        </ul>
    <li>[1..1] 53037-8^Genetic Sequence Variation Clinical Significance^LN (SubID=2a):</li>
        <ul>
        <li>If clinSig field is populated in annotation file: clinSig</li>
        <li>else: 'not specified'</li>
        </ul>
    <li>[1..1] 69548-6^Genetic Variant Assessment^LN (SubID=2a): 'Present'</li>
    <li>[0..1] 81259-4^Probable Associated Phenotype^LN (SubID=2a): include if phenotype supplied</li>
        <ul>
        <li>If phenotype field is populated in annotation file: phenotype</li>
        <li>else: omit OBX</li>
        </ul>
    <li>[0..1] 81258-6^Allelic frequency^LN (SubID=2a): Equals FORMAT.AD/FORMAT.DP</li>
    <li>If sourceClass = germline:</li>
        <ul>
        <li>[0..1] 53034-5^Allelic state^LN (SubID=2a): based on FORMAT:GT</li>
            <ul>
            <li>LA6703-8^Heteroplasmic^LN</li>
            <li>LA6704 6^Homoplasmic^LN</li>
            <li>LA6705-3^Homozygous^LN</li>
            <li>LA6706-1^Heterozygous^LN</li>
            <li>LA6707-9^Hemizygous^LN</li>
            </ul>
        </ul>
    </ul>
</ul>

**Example**:

```
chr5    112841059    .    T    A    .    .    .    GT    1/1
```

**first group (subID=2a)**

```
OBX|9|ST|83005-9^Variant Category^LN|2a|Simple|
OBX|10|ST|47998-0^Variant Display Name^LN|2a|NC_000005.10:112841058:T:A|
OBX|11|CWE|48018-6^Gene Studied^LN|2a|HGNC:583^APC^HGNC|
OBX|12|ST|48004-6^DNA Change c.HGVS^LN|2a|NM_001127510.3:c.5465T>A|
OBX|13|ST|48005-3^Amino Acid Change p.HGVS^LN|2a|p.Val1822Asp|
OBX|14|ST|48013-7^Genomic reference sequence^LN|2a|NC_000005.10|
OBX|15|ST|69547-8^Genomic ref allele^LN|2a|T|
OBX|16|NR|81254-5^Genomic allele start-end^LN|2a|112122384|
OBX|17|ST|69551-0^Genomic alt allele^LN|2a|G|
OBX|18|ST|48002-0^Genomic Source Class [Type]^LN|2a|LA6683-2^Germline^LN|
OBX|19|ST|53037-8^Genetic Sequence Variation Clinical Significance^LN|2a|Benign|
OBX|20|ST|69548-6^Genetic Variant Assessment^LN|2a|Present|
OBX|21|CWE|81259-4^Probable Associated Phenotype^LN|2a|72900001^Familial multiple polyposis syndrome^SCT|
OBX|22|CNE|53034-5^Allelic state^LN|2a|LA6705-3^Homozygous^LN|
```

<hr>

```
chr19    38499670    .    C    T    .    .    .    GT    0/1
```

**second group (subID=2b)**

```
OBX|23|ST|83005-9^Variant Category^LN|2b|Simple|
OBX|24|ST|47998-0^Variant Display Name^LN|2b|NC_000019.10:38499669:C:T|
OBX|25|ST|48018-6^Gene Studied^LN|2b|HGNC:10483^RYR1^HGNC|
OBX|26|ST|48004-6^DNA Change c.HGVS^LN|2b|NM_001042723.2:c.7063C>T|
OBX|27|ST|48005-3^Amino Acid Change p.HGVS^LN|2b|p.Arg2355Trp|
OBX|28|ST|48013-7^Genomic reference sequence^LN|2b|NC_000019.10|
OBX|29|ST|69547-8^Genomic ref allele^LN|2b|C|
OBX|30|NR|81254-5^Genomic allele start-end^LN|2b|41254965|
OBX|31|ST|69551-0^Genomic alt allele^LN|2b|CT|
OBX|32|ST|48002-0^Genomic Source Class [Type]^LN|2b|LA6683-2^Germline^LN|
OBX|33|ST|53037-8^Genetic Sequence Variation Clinical Significance^LN|2b|Pathogenic|
OBX|34|ST|69548-6^Genetic Variant Assessment^LN|2b|Present|
OBX|35|CWE|81259-4^Probable Associated Phenotype^LN|2b|405501007^Malignant hyperthermia^SCT|
OBX|36|CNE|53034-5^Allelic state^LN|2b|LA6706-1^Heterozygous^LN|
```
