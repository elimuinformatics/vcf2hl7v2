import logging
from .gene_ref_seq import _get_ref_seq_by_chrom
from .hl7v2_helper import _HL7V2_Helper
from .common import *

invalid_record_logger = logging.getLogger("vcf2hl7v2.invalidrecord")
general_logger = logging.getLogger("vcf2hl7v2.general")


def _valid_record(record):
    if not (validate_chrom_identifier(record.CHROM)):
        invalid_record_logger.debug(
            ("Reason: VCF CHROM is not recognized, " +
             "Record: %s, considered sample: %s"),
            record,
            record.samples[0].data)
        return False
    if(record.is_sv):
        invalid_record_logger.debug(
            ("Reason: VCF INFO.SVTYPE is present. " +
             "(Structural variants are excluded), " +
             "Record: %s, considered sample: %s"),
            record,
            record.samples[0].data)
        return False
    if(record.FILTER is not None and len(record.FILTER) != 0):
        invalid_record_logger.debug(
            ("Reason: VCF FILTER does not equal " +
             "'PASS' or '.', Record: %s, considered sample: %s"),
            record,
            record.samples[0].data)
        return False
    if(len(record.samples) == 0 or record.samples[0].gt_type is None):
        invalid_record_logger.debug(
            ("Reason: VCF FORMAT.GT is null ('./.', '.|.', '.', etc), " +
             "Record: %s, considered sample: %s"),
            record,
            record.samples[0].data)
        return False
    if(len(record.ALT) != 1 or
       (record.ALT[0].type != 'SNV' and record.ALT[0].type != 'MNV')):
        invalid_record_logger.debug(
            ("Reason: ALT is not a simple character string, comma-separated " +
             "character string or '.' (Rows where ALT is a token " +
             "(e.g. 'DEL') are excluded), Record: %s, considered sample: %s"),
            record,
            record.samples[0].data)
        return False
    if(not record.REF.isalpha()):
        invalid_record_logger.debug(
            ("Reason: REF is not a simple character string, " +
             "Record: %s, considered sample: %s"),
            record,
            record.samples[0].data)
    if record.CHROM == "M" and (
        (len(
            record.samples[0].gt_alleles) == 1 and
            record.samples[0].gt_alleles[0] == "0") or len(
            record.samples[0].gt_alleles) == 2):
        invalid_record_logger.debug(
            ("Reason: Mitochondrial DNA with GT = 0 or its diploid, " +
             "Record: %s, considered sample: %s"),
            record,
            record.samples[0].data)
        return False
    return True


def _get_chrom(chrom_index):
    switcher = {
        23: 'X',
        24: 'Y',
        25: 'M'
    }
    return switcher.get(chrom_index, str(chrom_index))


def _fix_regions_chrom(region):
    if region:
        region.Chromosome = region.Chromosome.apply(
            extract_chrom_identifier)


def _add_record_variants(
        record, ref_seq, patientID,
        hl7v2_helper, ratio_ad_dp, source_class,
        annotations, ref_build, vcf_type, variant_analysis_method, xml=False):
    spdi_representation = get_spdi_representation(record, ref_seq, xml)
    annotation_record =\
        get_annotations(
            record, annotations, spdi_representation, vcf_type, xml)

    if xml or (annotation_record is not None and _valid_record(record)):
        hl7v2_helper.add_variant_obv(record, ref_seq, ratio_ad_dp,
                                     source_class, annotation_record,
                                     spdi_representation, ref_build,
                                     variant_analysis_method, xml)


def _add_region_studied(
        region_studied, hl7v2_helper, chrom, ref_seq, patientID):
    if(region_studied and
       not region_studied[chrom].empty):
        general_logger.info("Adding region Studied OBX for %s", chrom)
        general_logger.debug("Region Examined %s", region_studied[chrom])
        hl7v2_helper.add_regionstudied_obv(
            ref_seq, region_studied[chrom])


def create_region_studied_obxs(
        vcf_reader_list, xml_reader, ref_build, patientID, has_tabix,
        conversion_region, region_studied, ratio_ad_dp, source_class,
        annotations, report, output_filename, hl7v2_helper):
    _fix_regions_chrom(conversion_region)
    _fix_regions_chrom(region_studied)
    if conversion_region:
        if region_studied:
            region_studied = region_studied.intersect(conversion_region)
            general_logger.debug(
                "Final Conmputed Reportable Query Regions: %s", region_studied)
    general_logger.info("Start adding the Region studied OBXs")
    if xml_reader is None:
        if has_tabix:
            for chrom_index in range(1, 26):
                chrom = _get_chrom(chrom_index)
                ref_seq = _get_ref_seq_by_chrom(
                    ref_build, extract_chrom_identifier(chrom))
                _add_region_studied(
                    region_studied, hl7v2_helper, chrom, ref_seq, patientID)
            hl7v2_helper.add_final_region_studied_obx_segment(report)
            return
        else:
            variants = vcf_reader_list
    else:
        if xml_reader.get("variant") is not None:
            if not isinstance(xml_reader.get("variant"), list):
                xml_reader["variant"] = [xml_reader.get("variant")]
            variants = xml_reader.get("variant")
        else:
            return
    chrom_index = 1
    prev_add_chrom = ""
    for variant in variants:
        CHROM = extract_chrom_identifier(get_chromosome(variant, xml_reader))
        if(not (conversion_region and
           conversion_region[CHROM].empty)):
            if prev_add_chrom != CHROM and (
                    region_studied):
                chrom = _get_chrom(chrom_index)
                while prev_add_chrom != CHROM:
                    current_ref_seq = _get_ref_seq_by_chrom(
                        ref_build, chrom)
                    _add_region_studied(
                        region_studied, hl7v2_helper, chrom,
                        current_ref_seq, patientID)
                    prev_add_chrom = chrom
                    chrom_index += 1
                    chrom = _get_chrom(chrom_index)
    hl7v2_helper.add_final_region_studied_obx_segment(report)


def create_variant_obxs(
        vcf_reader_list, vcf_reader, xml_reader, ref_build, patientID,
        has_tabix, conversion_region, ratio_ad_dp, source_class, annotations,
        vcf_type, variant_analysis_method, output_filename, hl7v2_helper):
    _fix_regions_chrom(conversion_region)
    general_logger.info("Start adding the Variant OBXs")
    if xml_reader is None:
        if has_tabix:
            for chrom_index in range(1, 26):
                chrom = _get_chrom(chrom_index)
                ref_seq = _get_ref_seq_by_chrom(
                    ref_build, extract_chrom_identifier(chrom))
                if conversion_region and not conversion_region[chrom].empty:
                    for _, row in conversion_region[chrom].df.iterrows():
                        vcf_iterator = None
                        try:
                            vcf_iterator = vcf_reader.fetch(
                                chrom, int(row['Start']), int(row['End']))
                        except ValueError:
                            pass
                        if vcf_iterator:
                            for record in vcf_iterator:
                                record.CHROM = extract_chrom_identifier(
                                    record.CHROM)
                                _add_record_variants(
                                    record, ref_seq, patientID, hl7v2_helper,
                                    ratio_ad_dp, source_class, annotations,
                                    ref_build, vcf_type,
                                    variant_analysis_method)
                elif not conversion_region:
                    vcf_iterator = None
                    try:
                        vcf_iterator = vcf_reader.fetch(chrom)
                    except ValueError:
                        pass
                    if vcf_iterator:
                        for record in vcf_iterator:
                            record.CHROM = extract_chrom_identifier(
                                record.CHROM)
                            _add_record_variants(
                                record, ref_seq, patientID, hl7v2_helper,
                                ratio_ad_dp, source_class,
                                annotations, ref_build, vcf_type,
                                variant_analysis_method)
            return
        else:
            variants = vcf_reader_list
    else:
        if xml_reader.get("variant") is not None:
            if not isinstance(xml_reader.get("variant"), list):
                xml_reader["variant"] = [xml_reader.get("variant")]
                print(xml_reader["variant"])
            variants = xml_reader.get("variant")
        else:
            return

    for variant in variants:
        CHROM = extract_chrom_identifier(get_chromosome(variant, xml_reader))
        if(not (conversion_region and
           conversion_region[CHROM].empty)):
            ref_seq = _get_ref_seq_by_chrom(ref_build, CHROM)
            POS = get_position(variant, xml_reader)
            REF = get_reference(variant, xml_reader)
            end = POS + len(REF) - 1
            if(not conversion_region or
               conversion_region[CHROM, POS - 1: end].empty is False):
                _add_record_variants(
                    variant, ref_seq, patientID, hl7v2_helper,
                    ratio_ad_dp, source_class,
                    annotations, ref_build, vcf_type,
                    variant_analysis_method, xml_reader is not None)


def _get_hl7v2_message(
        vcf_reader, xml_reader, ref_build, patientID, has_tabix,
        conversion_region, source_class, region_studied,
        ratio_ad_dp, annotations, seed, vcf_type,
        variant_analysis_method, report, output_filename):
    hl7v2_helper = _HL7V2_Helper(patientID, seed)
    general_logger.debug("Finished Initializing empty HL7V2 message")
    vcf_reader_list = list(vcf_reader) if vcf_reader else []
    create_region_studied_obxs(
        vcf_reader_list, xml_reader, ref_build, patientID, has_tabix,
        conversion_region, region_studied, ratio_ad_dp, source_class,
        annotations, report, output_filename, hl7v2_helper)
    if xml_reader is not None:
        hl7v2_helper.add_section1_components(xml_reader)
    create_variant_obxs(
        vcf_reader_list, vcf_reader, xml_reader, ref_build, patientID,
        has_tabix, conversion_region, ratio_ad_dp, source_class,
        annotations, vcf_type, variant_analysis_method, output_filename,
        hl7v2_helper)
    if hl7v2_helper.message.obx[-1].obx_4[0].value[0] == "2":
        hl7v2_helper.message.obx[hl7v2_helper.final_rs_index].obx_5 =\
            'Positive'
    general_logger.info(
        f'Export the HL7V2 message object to the file {output_filename}')
    hl7v2_helper.export_hl7v2_message(output_filename)
    general_logger.info("Completed conversion")
    return hl7v2_helper.message
