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
    if(
            len(record.ALT) != 1 or
            (record.ALT[0].type != 'SNV' and record.ALT[0].type != 'MNV')):
        invalid_record_logger.debug(
            ("Reason: ALT is not a simple character string, comma-separated " +
             "character string or '.' (Rows where ALT is a token " +
             "(e.g. 'DEL') are excluded), Record: %s, considered sample: %s"),
            record,
            record.samples[0].data)
        return False
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
        record,
        ref_seq,
        patientID,
        hl7v2_helper,
        ratio_ad_dp,
        source_class,
        annotations,
        ref_build):
    if(_valid_record(record)):
        hl7v2_helper.add_variant_obv(
                        record,
                        ref_seq,
                        ratio_ad_dp,
                        source_class,
                        annotations,
                        ref_build
                    )


def _add_region_studied(
        region_studied,
        hl7v2_helper,
        chrom,
        ref_seq,
        patientID):
    if(
            region_studied and
            not region_studied[chrom].empty):
        general_logger.info("Adding region Studied OBX for %s", chrom)
        general_logger.debug("Region Examined %s", region_studied[chrom])
        hl7v2_helper.add_regionstudied_obv(
            ref_seq, region_studied[chrom])


def create_region_studied_obxs(
        vcf_reader,
        ref_build,
        patientID,
        has_tabix,
        conversion_region,
        region_studied,
        ratio_ad_dp,
        source_class,
        annotations,
        output_filename,
        hl7v2_helper):
    _fix_regions_chrom(conversion_region)
    _fix_regions_chrom(region_studied)
    if conversion_region:
        if region_studied:
            region_studied = region_studied.intersect(conversion_region)
            general_logger.debug(
                "Final Conmputed Reportable Query Regions: %s", region_studied)
    general_logger.info("Start adding the Region studied OBXs")
    if has_tabix:
        for chrom_index in range(1, 26):
            chrom = _get_chrom(chrom_index)
            ref_seq = _get_ref_seq_by_chrom(
                ref_build, extract_chrom_identifier(chrom))
            _add_region_studied(region_studied,
                                hl7v2_helper, chrom, ref_seq, patientID)
    else:
        chrom_index = 1
        prev_add_chrom = ""
        for record in vcf_reader:
            record.CHROM = extract_chrom_identifier(record.CHROM)
            if not (
                    conversion_region and
                    conversion_region[record.CHROM].empty):
                if prev_add_chrom != record.CHROM and (
                        region_studied or nocall_region):
                    chrom = _get_chrom(chrom_index)
                    while prev_add_chrom != record.CHROM:
                        current_ref_seq = _get_ref_seq_by_chrom(
                            ref_build, chrom)
                        _add_region_studied(
                            region_studied,
                            hl7v2_helper,
                            chrom,
                            current_ref_seq,
                            patientID)
                        prev_add_chrom = chrom
                        chrom_index += 1
                        chrom = _get_chrom(chrom_index)
    hl7v2_helper.add_final_region_studied_obx_segment()


def create_variant_obxs(
        vcf_reader,
        ref_build,
        patientID,
        has_tabix,
        conversion_region,
        ratio_ad_dp,
        source_class,
        annotations,
        output_filename,
        hl7v2_helper):
    _fix_regions_chrom(conversion_region)
    general_logger.info("Start adding the Variant OBXs")
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
                                record,
                                ref_seq,
                                patientID,
                                hl7v2_helper,
                                ratio_ad_dp,
                                source_class,
                                annotations,
                                ref_build
                            )
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
                            record,
                            ref_seq,
                            patientID,
                            hl7v2_helper,
                            ratio_ad_dp,
                            source_class,
                            annotations,
                            ref_build
                        )
    else:
        for record in vcf_reader:
            record.CHROM = extract_chrom_identifier(record.CHROM)
            if not (
                    conversion_region and
                    conversion_region[record.CHROM].empty):
                ref_seq = _get_ref_seq_by_chrom(ref_build, record.CHROM)
                if (
                        not conversion_region or
                        conversion_region[
                            record.CHROM,
                            record.POS - 1: record.POS
                        ].empty is False):
                    _add_record_variants(
                        record,
                        ref_seq,
                        patientID,
                        hl7v2_helper,
                        ratio_ad_dp,
                        source_class,
                        annotations,
                        ref_build
                    )


def _get_hl7v2_message(
        vcf_reader,
        ref_build,
        patientID,
        has_tabix,
        conversion_region,
        region_studied,
        ratio_ad_dp,
        source_class,
        annotations,
        seed,
        output_filename):
    hl7v2_helper = _HL7V2_Helper(patientID, seed)
    general_logger.debug("Finished Initializing empty HL7V2 message")
    create_region_studied_obxs(
        vcf_reader,
        ref_build,
        patientID,
        has_tabix,
        conversion_region,
        region_studied,
        ratio_ad_dp,
        source_class,
        annotations,
        output_filename,
        hl7v2_helper
    )
    if annotations is not None:
        create_variant_obxs(
            vcf_reader,
            ref_build,
            patientID,
            has_tabix,
            conversion_region,
            ratio_ad_dp,
            source_class,
            annotations,
            output_filename,
            hl7v2_helper
        )
    if (hl7v2_helper.index + seed - hl7v2_helper.final_rs_index) > 1:
        hl7v2_helper.message.obx[hl7v2_helper.final_rs_index].obx_5 =\
            'Positive'
    general_logger.info(
        f'Export the HL7V2 message object to the file {output_filename}')
    hl7v2_helper.export_hl7v2_message(output_filename)
    general_logger.info("Completed conversion")
