import hl7apy.core as hl7
from .common import *
import pandas as pd


GENE_KB = pd.read_csv(__file__[:-15]+'geneKB_Gene_table.csv')


def _create_variant_obx_segment(self, obx_1, obx_2, obx_3, obx_4, obx_5):
    obx = hl7.Segment("OBX")
    obx.obx_1 = obx_1
    obx.obx_2 = obx_2
    obx.obx_3 = obx_3
    obx.obx_4 = obx_4
    obx.obx_5 = obx_5
    self.message.add(obx)
    self.index += 1


def _get_gene_studied(ref_seq, pos, build):
    if build == 'GRCh37':
        df = GENE_KB[
                (GENE_KB['Hg19/build37_refSeq'] == ref_seq) &
                (GENE_KB['Hg19/build37_start'] <= pos) &
                (GENE_KB['Hg19/build37_end'] >= pos)
            ]
    elif build == 'GRCh38':
        df = GENE_KB[
                (GENE_KB['Hg38/build38_refSeq'] == ref_seq) &
                (GENE_KB['Hg38/build38_start'] <= pos) &
                (GENE_KB['Hg38/build38_end'] >= pos)
            ]

    if (len(df['gene']) == 0):
        return 'HGNC:0000^NoGene^HGNC'

    for row in df['gene']:
        return str(row)


def _get_annotations(record, annotations, spdi_representation):
    if isinstance(annotations, str):
        if annotations == 'Not Supplied':
            return [spdi_representation, None, "not specified", None]
    annotation = annotations[
                    (annotations['CHROM'] == f'chr{record.CHROM}') &
                    (annotations['POS'] == record.POS)
                ]
    if len(annotation) == 1:
        if(not annotation['cHGVS'].empty and
           not annotation['transcriptRefSeq'].empty and
           not pd.isna(annotation.iloc[0]['cHGVS']) and
           not pd.isna(annotation.iloc[0]['transcriptRefSeq'])):
            dna_change = (f'{annotation.iloc[0]["transcriptRefSeq"]}:' +
                          f'{annotation.iloc[0]["cHGVS"]}')
        elif(not annotation['cHGVS'].empty and
             (annotation['transcriptRefSeq'].empty or
             pd.isna(annotation.iloc[0]['transcriptRefSeq'])) and
             not pd.isna(annotation.iloc[0]['cHGVS'])):
            dna_change = f'{annotation.iloc[0]["cHGVS"]}'
        else:
            dna_change = spdi_representation
        if(not annotation['pHGVS'].empty and
           not annotation['proteinRefSeq'].empty and
           not pd.isna(annotation.iloc[0]['pHGVS']) and
           not pd.isna(annotation.iloc[0]['proteinRefSeq'])):
            amino_acid_change = (f'{annotation.iloc[0]["proteinRefSeq"]}:' +
                                 f'{annotation.iloc[0]["pHGVS"]}')
        elif(not annotation['pHGVS'].empty and
             (annotation['proteinRefSeq'].empty or
             pd.isna(annotation.iloc[0]['proteinRefSeq'])) and
             not pd.isna(annotation.iloc[0]['pHGVS'])):
            amino_acid_change = f'{annotation.iloc[0]["pHGVS"]}'
        else:
            amino_acid_change = None
        if(not annotation['clinSig'].empty and
           not pd.isna(annotation.iloc[0]['clinSig'])):
            clinSig = annotation.iloc[0]['clinSig']
        else:
            clinSig = "not specified"
        if(not annotation['phenotype'].empty and
           not pd.isna(annotation.iloc[0]['phenotype'])):
            phenotype = annotation.iloc[0]['phenotype']
        else:
            phenotype = None
        return [dna_change, amino_acid_change, clinSig, phenotype]
    else:
        return None


class _HL7V2_Helper:
    def __init__(self, patientID, seed):
        self.message = hl7.Message("ORU_R01")
        self.patientID = patientID
        self.seed = seed
        self.final_rs_index = None
        self.index = seed
        self.region_studied_index = 1
        self.variant_index = 1

    def add_regionstudied_obv(self, ref_seq, reportable_query_regions):
        if reportable_query_regions.empty:
            return
        obx = hl7.Segment("OBX")
        obx.obx_1 = str(self.index)
        obx.obx_2 = "ST"
        obx.obx_3 = "48013-7^Genomic reference sequence^LN"
        obx.obx_4 = f'1{get_alphabet_index(self.region_studied_index)}'
        obx.obx_5 = f'{ref_seq}'
        self.message.add(obx)
        self.index += 1
        a_index = 1
        for _, row in reportable_query_regions.df.iterrows():
            obx = hl7.Segment("OBX")
            obx.obx_1 = str(self.index)
            obx.obx_2 = "NR"
            obx.obx_3 = "51959-5^Range(s) of DNA sequence examined^LN"
            obx.obx_4 = (f'1{get_alphabet_index(self.region_studied_index)}.' +
                         f'{get_alphabet_index(a_index)}')
            obx.obx_5 = f'{row["Start"]}^{row["End"]}'
            self.message.add(obx)
            a_index += 1
            self.index += 1
        self.region_studied_index += 1

    def add_final_region_studied_obx_segment(self):
        obx = hl7.Segment("OBX")
        obx.obx_1 = str(self.index)
        obx.obx_2 = "TX"
        obx.obx_3 = ("51968-6^Discrete variation " +
                     "analysis overall interpretation^LN")
        obx.obx_4 = '1'
        obx.obx_5 = 'Negative'
        self.message.add(obx)
        self.final_rs_index = self.index - self.seed
        self.index += 1

    def add_variant_obv(
            self, record,
            ref_seq, ratio_ad_dp, source_class, annotations, ref_build):
        alleles = get_allelic_state(record, ratio_ad_dp)
        if alleles["CODE"] == '' and alleles["ALLELE"] == '':
            allelic_state = None
        else:
            allelic_state = f'{alleles["CODE"]}^{alleles["ALLELE"].title()}^LN'
        a_index = f'2{get_alphabet_index(self.variant_index)}'
        spdi_representation = (f'{ref_seq}:{record.POS - 1}:{record.REF}:' +
                               f'{"".join(list(map(str, list(record.ALT))))}')
        gene_studied = _get_gene_studied(ref_seq, int(record.POS), ref_build)
        result = _get_annotations(record, annotations, spdi_representation)
        if result is None:
            return
        dna_change, amino_acid_change, clinSig, phenotype = result
        _create_variant_obx_segment(
            self, obx_1=str(self.index), obx_2="ST",
            obx_3="83005-9^Variant Category^LN", obx_4=a_index, obx_5="Simple",
        )
        _create_variant_obx_segment(
            self, obx_1=str(self.index),
            obx_2="ST", obx_3="47998-0^Variant Display Name^LN",
            obx_4=a_index, obx_5=spdi_representation,
        )
        _create_variant_obx_segment(
            self, obx_1=str(self.index), obx_2="CWE",
            obx_3="48018-6^Gene Studied^LN", obx_4=a_index, obx_5=gene_studied,
        )
        _create_variant_obx_segment(
            self, obx_1=str(self.index), obx_2="ST",
            obx_3="48004-6^DNA Change c.HGVS^LN",
            obx_4=a_index, obx_5=dna_change,
        )
        if amino_acid_change is not None:
            _create_variant_obx_segment(
                self, obx_1=str(self.index),
                obx_2="ST", obx_3="48005-3^Amino Acid Change p.HGVS^LN",
                obx_4=a_index, obx_5=amino_acid_change,
            )
        _create_variant_obx_segment(
            self, obx_1=str(self.index),
            obx_2="ST", obx_3="48013-7^Genomic reference sequence^LN",
            obx_4=a_index, obx_5=f'{ref_seq}',
        )
        _create_variant_obx_segment(
            self, obx_1=str(self.index),
            obx_2="ST", obx_3="69547-8^Genomic ref allele^LN",
            obx_4=a_index, obx_5=str(record.REF),
        )
        _create_variant_obx_segment(
            self, obx_1=str(self.index),
            obx_2="NR", obx_3="81254-5^Genomic allele start-end^LN",
            obx_4=a_index, obx_5=str(record.POS),
        )
        _create_variant_obx_segment(
            self, obx_1=str(self.index),
            obx_2="ST", obx_3="69551-0^Genomic alt allele^LN",
            obx_4=a_index, obx_5=("".join(list(map(str, list(record.ALT))))),
        )
        _create_variant_obx_segment(
            self, obx_1=str(self.index),
            obx_2="ST", obx_3="48002-0^Genomic Source Class [Type]^LN",
            obx_4=a_index, obx_5=source_class,
        )
        _create_variant_obx_segment(
            self, obx_1=str(self.index),
            obx_2="ST", obx_3=("53037-8^Genetic Sequence " +
                               "Variation Clinical Significance^LN"),
            obx_4=a_index, obx_5=clinSig,
        )
        _create_variant_obx_segment(
            self, obx_1=str(self.index),
            obx_2="ST", obx_3="69548-6^Genetic Variant Assessment^LN",
            obx_4=a_index, obx_5="Present",
        )
        if phenotype is not None:
            _create_variant_obx_segment(
                self, obx_1=str(self.index),
                obx_2="ST", obx_3="81259-4^Probable Associated Phenotype^LN",
                obx_4=a_index, obx_5=phenotype,
            )
        if allelic_state is not None:
            _create_variant_obx_segment(
                self, obx_1=str(self.index),
                obx_2="CNE", obx_3="53034-5^Allelic state^LN",
                obx_4=a_index, obx_5=allelic_state,
            )
        self.variant_index += 1

    def export_hl7v2_message(self, output_filename):
        with open(output_filename, 'w') as output_file:
            for segment in self.message.obx:
                output_file.write(segment.to_er7())
                output_file.write('|\n')
