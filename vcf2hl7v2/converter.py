import vcf
import pyranges
from .hl7v2_message_generator import _get_hl7v2_message
import logging
import sys
import pandas as pd
from .common import *
general_logger = logging.getLogger('vcf2hl7v2.general')
"""Converter for a VCF version >4.1 file."""


class Converter(object):
    """Creates a new Converter Object to convert a VCF file.

    Parameters
    ----------

    **vcf_filename** (required): Path to a text-based or bgzipped VCF
    file. Valid path and filename without whitespace must be provided.
    VCF file must conform to VCF Version 4.1 or later. FORMAT.GT must be
    present. Multi-sample VCFs are allowed, but only the first sample
    will be converted.
    has_tabix (required if VCF file is bgzipped): Set to 'True' if
    there is a tabix index. Tabix file must have the same name as the
    bgzipped VCF file, with a '.tbi' extension, and must be in the same
    folder.

    **ref_build** (required): Genome Reference Consortium genome
    assembly to which variants in the VCF were called. Must be one of
    'GRCh37' or 'GRCh38'.

    **patient_id** (optional): Supplied patient ID is inserted into
    generated HL7V2 output. Alphanumeric string without whitespace. if \
    not provided,
    header of first sample column is used.

    **Conversion region** (optional): Subset of the VCF file to be
    converted into HL7V2. If absent, the entire VCF file is converted. Can
    be supplied as either a parameter (conv_region_dict) or as a BED
    file (conv_region_filename):

       **conv_region_dict** : Array of regions (e.g. '{"Chromosome":
       ["X", "X", "M"],"Start": [50000, 55000, 50000],"End": [52000,
       60600, 60025]}'). Values for Chromosome must align with values in
       VCF #CHROM field. Ranges must be `0-based \
       <https://www.biostars.org/p/84686/>`_
       (or 0-start, half-open) and based on GRCh37 or GRCh38 \
       reference sequences.

       **conv_region_filename**: Valid path and filename without
       whitespace must be provided. Must be a valid BED file with first 3
       columns: <chr> <start> <stop>. Values in <chr> field must align
       with values in VCF #CHROM field. Ranges must be based on GRCh37 or
       GRCh38 reference sequences. Note that BED files are
       `0-based <https://www.biostars.org/p/84686/>`_ (or 0-start,
       half-open) whereas VCF files and HL7V2 output are 1-based (or
       1-start, fully-closed).

    **region_studied_filename** (optional): Subset of patient's genome
    that was studied in the generation of the VCF file. Valid path and
    filename without whitespace must be provided. Must be a valid BED
    file, with first 3 columns: <chr> <start> <stop>. Values in <chr>
    field must align with values in VCF #CHROM field. Ranges must be
    based on GRCh37 or GRCh38 reference sequences. Note that BED files
    are `0-based <https://www.biostars.org/p/84686/>`_ (or 0-start,
    half-open) whereas VCF files and HL7V2 output are 1-based (or 1-start,
    fully-closed).

    **ratio_ad_dp** (optional)(default value = 0.99): This ratio \
    determine whether to assign Homoplasmic or Heteroplasmic.
    If allelic depth (FORMAT.AD) / read depth (FORMAT.DP) is \
    greater than ratio_ad_dp then allelic state is
    homoplasmic; else heteroplasmic.

    **source_class** (optional)(default value = germline):

    **seed** (optional)(default value = 1000):

    **annotation_filename** (optional):
    Returns
    -------

    Object
    An Instance of Converter that helps to convert vcf file.

    """

    def __init__(
            self, vcf_filename=None, ref_build=None, patient_id=None,
            has_tabix=False, conv_region_filename=None, conv_region_dict=None,
            region_studied_filename=None, ratio_ad_dp=0.99,
            source_class='germline', seed=1000, annotation_filename=None):

        super(Converter, self).__init__()
        if not (vcf_filename):
            raise Exception('You must provide vcf_filename')
        if not ref_build or ref_build not in ["GRCh37", "GRCh38"]:
            raise Exception(
                'You must provide build number ("GRCh37" or "GRCh38")')
        self.vcf_filename = vcf_filename
        try:
            self._vcf_reader = list(vcf.Reader(filename=vcf_filename))
        except FileNotFoundError:
            raise
        except BaseException:
            self._generate_exception("Please provide valid 'vcf_filename'")
        self.annotation_filename = annotation_filename
        if self.annotation_filename is None:
            self.annotations = 'Not Supplied'
        else:
            try:
                self.annotations = pd.read_csv(
                            self.annotation_filename,
                            names=[
                                'CHROM', 'POS', 'REF', 'ALT',
                                'transcriptRefSeq', 'cHGVS', 'proteinRefSeq',
                                'pHGVS', 'clinSig', 'phenotype'
                            ],
                            sep='\t'
                        )
                if len(self.annotations) == 0:
                    self.annotations = None
            except FileNotFoundError:
                raise
            except BaseException:
                self._generate_exception(
                        "Please provide valid 'annotation_filename'"
                    )
        if not patient_id:
            patient_id = self._vcf_reader.samples[0]
        if conv_region_filename:
            try:
                self.conversion_region = pyranges.read_bed(
                    conv_region_filename)
            except FileNotFoundError:
                raise
            except BaseException:
                self._generate_exception(
                    "Please provide valid 'conv_region_filename'")
        elif conv_region_dict:
            try:
                self.conversion_region = pyranges.from_dict(conv_region_dict)
            except FileNotFoundError:
                raise
            except BaseException:
                self._generate_exception(
                    "Please provide valid 'conv_region_dict'")
        else:
            self.conversion_region = None
        if region_studied_filename:
            try:
                self.region_studied = pyranges.read_bed(
                    region_studied_filename)
            except FileNotFoundError:
                raise
            except BaseException:
                self._generate_exception(
                    "Please provide valid 'region_studied_filename'")
        else:
            self.region_studied = None

        if not validate_has_tabix(has_tabix):
            raise Exception("Please provide a valid 'has_tabix'")

        if not validate_ratio_ad_dp(ratio_ad_dp):
            raise Exception("Please provide a valid 'ratio_ad_dp'")

        if not validate_seed(seed):
            raise Exception("Please provide a valid seed")

        if source_class not in ["germline", "somatic"]:
            raise Exception(
                'Please provide a valid Source Class ("germline" or "somatic")'
            )

        self.ratio_ad_dp = ratio_ad_dp
        self.has_tabix = has_tabix
        self.patient_id = patient_id
        self.ref_build = ref_build
        self.conv_region_filename = conv_region_filename
        if source_class == 'germline':
            self.source_class = 'LA6683-2^Germline^LN'
        elif source_class == 'somatic':
            self.source_class = 'LA6684-0^Somatic^LN'
        self.seed = seed
        general_logger.info("Converter class instantiated successfully")

    def convert(self, output_filename='hl7v2.txt'):
        """ Generates HL7 V2 format data as output_filename \
        or hl7v2.txt if it is not provided

        Parameters
        ----------
        output_filename:
            Path to output hl7v2 text file.

        """
        general_logger.info("Starting VCF to HL7V2 Conversion")
        _get_hl7v2_message(
            self._vcf_reader, self.ref_build, self.patient_id, self.has_tabix,
            self.conversion_region, self.region_studied, self.ratio_ad_dp,
            self.source_class, self.annotations, self.seed, output_filename)
        general_logger.info("Completed VCF to HL7V2 Conversion")

    def _generate_exception(self, msg):
        general_logger.error(msg, exc_info=True)
        raise Exception(msg, sys.exc_info)
