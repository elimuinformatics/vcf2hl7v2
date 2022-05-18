import hl7apy.core as hl7
from .common import *
import pandas as pd


def _create_variant_obx_segment(self, obx_1, obx_2, obx_3, obx_4, obx_5):
    obx = hl7.Segment("OBX")
    obx.obx_1 = obx_1
    obx.obx_2 = obx_2
    obx.obx_3 = obx_3
    obx.obx_4 = obx_4
    obx.obx_5 = obx_5
    if len(str(obx)) > 50:
        prev_len =\
            len(str(obx_1)) + len(str(obx_2)) + len(str(obx_3)) +\
            len(str(obx_4))
        required_len = 50 - prev_len
        obx_5_str = str(obx_5)
        obx_5_str = obx_5_str[:required_len]
        obx.obx_5 = obx_5_str
    self.message.add(obx)
    self.index += 1


class _HL7V2_Helper:
    def __init__(self, patientID, seed):
        self.message = hl7.Message("ORU_R01")
        self.patientID = patientID
        self.seed = seed
        self.final_rs_index = None
        self.index = seed
        self.region_studied_index = 1
        self.variant_index = 1
        self.section1_index = 1

    def add_section1_components(self, xml_reader):
        obx = hl7.Segment("OBX")
        obx.obx_1 = str(self.index)
        obx.obx_2 = "ST"
        obx.obx_3 = "51967-8^Genetic Disease Assessed^UFMOLLRR"
        obx.obx_4 = '1'
        obx.obx_5 = f'{xml_reader.get("diagnosis")}'
        self.message.add(obx)
        self.index += 1
        self.section1_index += 1

    #     obx = hl7.Segment("OBX")
    #     obx.obx_1 = str(self.index)
    #     obx.obx_2 = "TX"
    #     obx.obx_3 = "51968-6^Genetic Analysis Overall Interpretation^LN"
    #     obx.obx_4 = '1'
    #     obx.obx_5 = f'{xml_reader.get("interpretation")}'
    #     self.message.add(obx)
    #     self.index += 1
    #     self.section1_index += 1

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

    def add_final_region_studied_obx_segment(self, report):
        obx = hl7.Segment("OBX")
        obx.obx_1 = str(self.index)
        obx.obx_2 = "TX"
        obx.obx_3 = ("51968-6^Genetic Analysis " +
                     "Overall Interpretation^UFMOLLRR")
        obx.obx_4 = '1'
        obx.obx_5 = 'Negative'
        self.message.add(obx)
        self.final_rs_index = self.index - self.seed
        self.index += 1

        if report is not None:
            obx = hl7.Segment("OBX")
            obx.obx_1 = str(self.index)
            obx.obx_2 = "TX"
            obx.obx_3 = ("51969-4^Genetic Analysis Report^UFMOLLRR")
            obx.obx_4 = '1'
            report_segment = ""
            for line in report:
                report_segment += line.replace('\n', '').strip()
                report_segment += " ~"
            obx.obx_5 = report_segment
            self.message.add(obx)
            self.index += 1

    def add_variant_obv(self, record, ref_seq, ratio_ad_dp, source_class,
                        annotation_record, spdi_representation, ref_build,
                        variant_analysis_method, xml):
        allelic_state, allelic_frequency = get_as_af(record, ratio_ad_dp, xml)
        variant_display_name =\
            get_variant_display_name(record, spdi_representation, xml)
        CHROM = get_chromosome(record, xml)
        POS = get_position(record, xml)
        REF = get_reference(record, xml)
        ALT = get_alternate(record, xml)

        if source_class is not None:
            if source_class.title() == Genomic_Source_Class.GERMLINE.value:
                source_class_value = 'LA6683-2^Germline^LN'
            elif source_class.title() == Genomic_Source_Class.SOMATIC.value:
                source_class_value = 'LA6684-0^Somatic^LN'

        vam_obx_value =\
            (f"{SEQUENCING_TO_CODE.get(variant_analysis_method.lower())}^" +
             f"{variant_analysis_method}^LN")

        a_index = f'2{get_alphabet_index(self.variant_index)}'
        _create_variant_obx_segment(
            self, obx_1=str(self.index), obx_2="CWE",
            obx_3="83005-9^Variant Category^LN", obx_4=a_index,
            obx_5="LA26801-3^Simple^LN")
        _create_variant_obx_segment(
            self, obx_1=str(self.index),
            obx_2="ST", obx_3="47998-0^Variant Display Name^LN",
            obx_4=a_index, obx_5=variant_display_name)
        _create_variant_obx_segment(
            self, obx_1=str(self.index), obx_2="CWE",
            obx_3="48018-6^Gene Studied^UFMOLLRR", obx_4=a_index,
            obx_5=annotation_record['gene_studied'])
        if annotation_record['transcript_ref_seq'] is not None:
            _create_variant_obx_segment(
                self, obx_1=str(self.index), obx_2="CWE",
                obx_3="51958-7^Transcript Reference Sequence^LN",
                obx_4=a_index,
                obx_5=annotation_record['transcript_ref_seq'])
        _create_variant_obx_segment(
            self, obx_1=str(self.index), obx_2="ST",
            obx_3="48004-6^DNA Change c.HGVS^LN",
            obx_4=a_index, obx_5=annotation_record['dna_change'])
        if annotation_record['amino_acid_change'] is not None:
            _create_variant_obx_segment(
                self, obx_1=str(self.index),
                obx_2="ST", obx_3="48005-3^Amino Acid Change p.HGVS^LN",
                obx_4=a_index, obx_5=annotation_record['amino_acid_change'])
        if annotation_record['molecular_consequence'] is not None:
            _create_variant_obx_segment(
                self, obx_1=str(self.index),
                obx_2="CWE", obx_3="48006-1^Molecular Consequence^UFMOLLRR",
                obx_4=a_index,
                obx_5=annotation_record['molecular_consequence'])
        _create_variant_obx_segment(
            self, obx_1=str(self.index),
            obx_2="CWE", obx_3="48013-7^Genomic reference sequence^UFMOLLRR",
            obx_4=a_index, obx_5=f'{ref_seq}')
        if annotation_record['genomic_dna_change'] is not None:
            _create_variant_obx_segment(
                self, obx_1=str(self.index),
                obx_2="ST", obx_3="81290-9^Genomic DNA Change g.HGVS^LN",
                obx_4=a_index, obx_5=annotation_record['genomic_dna_change'])
        _create_variant_obx_segment(
            self, obx_1=str(self.index),
            obx_2="ST", obx_3="69547-8^Genomic ref allele^LN",
            obx_4=a_index, obx_5=str(REF))
        _create_variant_obx_segment(
            self, obx_1=str(self.index),
            obx_2="NR", obx_3="81254-5^Genomic allele start-end^LN",
            obx_4=a_index, obx_5=str(POS))
        _create_variant_obx_segment(
            self, obx_1=str(self.index),
            obx_2="ST", obx_3="69551-0^Genomic alt allele^LN",
            obx_4=a_index, obx_5=ALT)
        _create_variant_obx_segment(
            self, obx_1=str(self.index),
            obx_2="ST", obx_3="48000-4^Chromosome^LN",
            obx_4=a_index, obx_5=CHROM)
        if annotation_record["dna_region"] is not None:
            _create_variant_obx_segment(
                self, obx_1=str(self.index),
                obx_2="ST", obx_3="47999-8^DNA Region^LN",
                obx_4=a_index, obx_5=annotation_record["dna_region"])
        # if annotation_record["protein_ref_seq"] is not None:
        #     _create_variant_obx_segment(
        #         self, obx_1=str(self.index),
        #         obx_2="ST", obx_3="^Protein Reference Sequence^LN",
        #         obx_4=a_index, obx_5=annotation_record["protein_ref_seq"])
        # if annotation_record["db_snp_id"] is not None:
        #     _create_variant_obx_segment(
        #         self, obx_1=str(self.index),
        #         obx_2="CWE", obx_3="81255-2^dbSNP [ID]^LN",
        #         obx_4=a_index, obx_5=annotation_record["db_snp_id"])
        if source_class is not None:
            _create_variant_obx_segment(
                self, obx_1=str(self.index),
                obx_2="CWE", obx_3="48002-0^Genomic Source Class [Type]^LN",
                obx_4=a_index, obx_5=source_class_value)
        _create_variant_obx_segment(
            self, obx_1=str(self.index),
            obx_2="CWE", obx_3="81304-8^Variant Analysis Method [Type]^LN",
            obx_4=a_index,
            obx_5=vam_obx_value)
        _create_variant_obx_segment(
            self, obx_1=str(self.index),
            obx_2="CWE", obx_3=("53037-8^Genetic Sequence " +
                                "Variation Clinical Significance^LN"),
            obx_4=a_index, obx_5=annotation_record['clin_sig'])
        _create_variant_obx_segment(
            self, obx_1=str(self.index),
            obx_2="CWE", obx_3="69548-6^Genetic Variant Assessment^LN",
            obx_4=a_index, obx_5="LA9633-4^Present^LN")
        if annotation_record['phenotype'] is not None:
            _create_variant_obx_segment(
                self, obx_1=str(self.index),
                obx_2="CWE", obx_3="81259-4^Probable Associated Phenotype^LN",
                obx_4=a_index, obx_5=annotation_record['phenotype'])
        # if annotation_record['phenotype_description'] is not None:
        #     _create_variant_obx_segment(
        #         self, obx_1=str(self.index),
        #         obx_2="ST", obx_3="^Phenotype Description^LN",
        #         obx_4=a_index,
        #         obx_5=annotation_record['phenotype_description'])
        if(allelic_state is not None and
           source_class is not None and
           source_class.title() == Genomic_Source_Class.GERMLINE.value):
            _create_variant_obx_segment(
                self, obx_1=str(self.index),
                obx_2="CWE", obx_3="53034-5^Allelic state^LN",
                obx_4=a_index, obx_5=allelic_state)
        if allelic_frequency is not None:
            _create_variant_obx_segment(
                self, obx_1=str(self.index),
                obx_2="NM", obx_3="81258-6^Allelic frequency^LN",
                obx_4=a_index, obx_5=f'{round(allelic_frequency, 2)}')
        if annotation_record["read_depth"] is not None:
            _create_variant_obx_segment(
                self, obx_1=str(self.index),
                obx_2="NM", obx_3="82121-5^Allelic Read Depth^LN",
                obx_4=a_index, obx_5=annotation_record["read_depth"])
        self.variant_index += 1

    def export_hl7v2_message(self, output_filename):
        with open(output_filename, 'w', encoding='utf-8') as output_file:
            for segment in self.message.obx:
                output_file.write(segment.to_er7().replace(r"\R\\", "~"))
                output_file.write('|\n')
