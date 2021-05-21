import doctest
import unittest
import os
import vcf2hl7v2
from os.path import join, dirname
import shutil
import logging
import filecmp

suite = doctest.DocTestSuite(vcf2hl7v2)


class TestTemporary(unittest.TestCase):
    @classmethod
    def setUpClass(self):
        self.TEST_RESULT_DIR = os.path.join(
            os.path.dirname(__file__), 'output')
        if os.path.exists(self.TEST_RESULT_DIR):
            shutil.rmtree(self.TEST_RESULT_DIR)
        os.mkdir(self.TEST_RESULT_DIR)

    @classmethod
    def tearDownClass(self):
        shutil.rmtree(self.TEST_RESULT_DIR)

    def test_NB6TK328(self):
        o_vcf_2_hl7v2 = vcf2hl7v2.Converter(
                vcf_filename=os.path.join(
                                        os.path.dirname(__file__),
                                        'NB6TK328_filtered.vcf'
                                    ),
                ref_build='GRCh38',
                patient_id='NB6TK328',
                has_tabix=False,
                conv_region_filename=os.path.join(
                                        os.path.dirname(__file__),
                                        'NB6TK328_conversion_region.bed'
                                    ),
                conv_region_dict=None,
                region_studied_filename=os.path.join(
                                        os.path.dirname(__file__),
                                        'NB6TK328_region_studied.bed'
                                    ),
                ratio_ad_dp=0.95,
                source_class='germline',
                seed=5000,
                annotation_filename=os.path.join(
                                        os.path.dirname(__file__),
                                        'NB6TK328_annotations.txt'
                                    ),
            )
        output_filename = os.path.join(os.path.dirname(
            __file__), self.TEST_RESULT_DIR, 'NB6TK328.txt')
        expected_output_filename = os.path.join(
            os.path.dirname(__file__), 'expected_NB6TK328.txt')
        o_vcf_2_hl7v2.convert(output_filename)
        self.assertEqual(
                filecmp.cmp(output_filename, expected_output_filename),
                True
            )
