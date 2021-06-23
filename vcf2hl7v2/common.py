import pandas as pd
import logging
import re

general_logger = logging.getLogger("vcf2hl7v2.general")


def get_allelic_state(record, ratio_ad_dp):
    allelic_state = ''
    allelic_code = ''
    # Using  the first sample
    sample = record.samples[0]
    alleles = sample.gt_alleles
    if record.CHROM != 'M':
        if len(alleles) >= 2 and sample.gt_type == 1:
            allelic_state = 'heterozygous'
            allelic_code = 'LA6706-1'
        elif len(alleles) >= 2:
            allelic_state = 'homozygous'
            allelic_code = 'LA6705-3'
        elif sample.gt_type is not None and len(alleles) == 1:
            allelic_state = 'hemizygous'
            allelic_code = 'LA6707-9'
        else:
            _error_log_allelicstate(record)
    elif (
            sample.gt_type is not None and
            len(alleles) == 1 and
            alleles[0] == '1'):
        if hasattr(sample.data, 'AD') and hasattr(sample.data, 'DP'):
            try:
                if(
                        isinstance(sample.data.AD, list) and
                        len(sample.data.AD) > 0):
                    ratio = float(
                        sample.data.AD[0]) / float(sample.data.DP)
                else:
                    ratio = float(sample.data.AD) / float(sample.data.DP)
                if ratio > ratio_ad_dp:
                    allelic_state = "homoplasmic"
                    allelic_code = "LA6704-6"
                else:
                    allelic_state = "heteroplasmic"
                    allelic_code = "LA6703-8"
            except Exception as e:
                general_logger.debug(e)
                _error_log_allelicstate(record)
                pass
        else:
            _error_log_allelicstate(record)
    else:
        _error_log_allelicstate(record)
    return {'ALLELE': allelic_state, 'CODE': allelic_code}


def extract_chrom_identifier(chrom):
    chrom = chrom.upper().replace("CHR", "")
    if chrom == "MT":
        chrom = "M"
    return chrom


def validate_chrom_identifier(chrom):
    chrom = extract_chrom_identifier(chrom)
    pattern = '^[1-9]$|^1[0-9]$|^2[0-2]$|^[XYM]$'
    result = re.match(pattern, chrom)
    return bool(result)


def validate_ratio_ad_dp(ratio_ad_dp):
    if not (ratio_ad_dp):
        return False
    if not isinstance(ratio_ad_dp, float):
        return False
    if ratio_ad_dp < 0 or ratio_ad_dp >= 1:
        return False
    return True


def validate_has_tabix(has_tabix):
    if not isinstance(has_tabix, bool):
        return False
    return True


def validate_seed(seed):
    if isinstance(seed, int) and not isinstance(seed, bool):
        return True
    return False


def get_alphabet_index(n):
    a_index = ""
    while(n != 0):
        remainder = n % 26
        n = n // 26
        if remainder == 0:
            remainder = 26
            n = n - 1
        a_index = chr(96 + remainder) + a_index
    return a_index


def _error_log_allelicstate(record):
    general_logger.error(
        "Cannot Determine AllelicState for: %s , considered sample: %s",
        record,
        record.samples[0].data)
