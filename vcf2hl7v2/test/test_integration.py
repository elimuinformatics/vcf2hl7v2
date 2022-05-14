import doctest
import unittest
import os
import vcf2hl7v2
from os.path import join, dirname
import shutil
import logging
import filecmp

suite = doctest.DocTestSuite(vcf2hl7v2)


class Testvcf2hl7v2Inputs(unittest.TestCase):

    def test_required_filename(self):
        with self.assertRaises(Exception) as context:
            vcf2hl7v2.Converter()
        self.assertTrue(
            'You must provide a vcf or a xml file' in str(context.exception))

    def test_required_ref_build(self):
        with self.assertRaises(Exception) as context:
            vcf2hl7v2.Converter(os.path.join(
                os.path.dirname(__file__), 'vcf_example1.vcf'))
        self.assertEqual(
            'You must provide build number ("GRCh37" or "GRCh38")', str(
                context.exception))

    def test_invalid_ref_build(self):
        with self.assertRaises(Exception) as context:
            vcf2hl7v2.Converter(os.path.join(
                os.path.dirname(__file__), 'vcf_example1.vcf'), 'b38')
        self.assertEqual(
            'You must provide build number ("GRCh37" or "GRCh38")', str(
                context.exception))

    def test_valid_ref_build_37(self):
        o_vcf_2_hl7v2 = vcf2hl7v2.Converter(os.path.join(
            os.path.dirname(__file__), 'vcf_example1.vcf'), 'GRCh37')
        self.assertEqual(type(o_vcf_2_hl7v2), vcf2hl7v2.Converter)

    def test_valid_ref_build_38(self):
        o_vcf_2_hl7v2 = vcf2hl7v2.Converter(os.path.join(
            os.path.dirname(__file__), 'vcf_example1.vcf'), 'GRCh38')
        self.assertEqual(type(o_vcf_2_hl7v2), vcf2hl7v2.Converter)

    def test_conv_region_only(self):
        conv_region_filename = os.path.join(os.path.dirname(
            __file__), 'RegionsToConvert_example3.bed')
        o_vcf_2_hl7v2 = vcf2hl7v2.Converter(
            os.path.join(
                os.path.dirname(__file__),
                'vcf_example3.vcf'),
            'GRCh37',
            'abc',
            conv_region_filename=conv_region_filename)
        self.assertEqual(type(o_vcf_2_hl7v2), vcf2hl7v2.Converter)

    def test_conv_region_dict(self):
        conv_region_dict = {
            "Chromosome": ["X", "X", "M"],
            "Start": [50000, 55000, 50000],
            "End": [52000, 60600, 60025]
        }
        o_vcf_2_hl7v2 = vcf2hl7v2.Converter(
            os.path.join(
                os.path.dirname(__file__),
                'vcf_example3.vcf'),
            'GRCh37',
            'abc',
            conv_region_dict=conv_region_dict)
        self.assertEqual(type(o_vcf_2_hl7v2), vcf2hl7v2.Converter)

    def test_conv_region_region_studied(self):
        region_studied_filename = os.path.join(
            os.path.dirname(__file__), 'RegionsStudied_example3.bed')
        conv_region_filename = os.path.join(os.path.dirname(
            __file__), 'RegionsToConvert_example3.bed')
        o_vcf_2_hl7v2 = vcf2hl7v2.Converter(
            os.path.join(
                os.path.dirname(__file__),
                'vcf_example3.vcf'),
            'GRCh37',
            'abc',
            conv_region_filename=conv_region_filename,
            region_studied_filename=region_studied_filename)
        self.assertEqual(type(o_vcf_2_hl7v2), vcf2hl7v2.Converter)

    def test_no_conv_region_region_studied(self):
        region_studied_filename = os.path.join(
            os.path.dirname(__file__), 'RegionsStudied_example3.bed')
        o_vcf_2_hl7v2 = vcf2hl7v2.Converter(
            os.path.join(
                os.path.dirname(__file__),
                'vcf_example3.vcf'),
            'GRCh37',
            'abc',
            region_studied_filename=region_studied_filename)
        self.assertEqual(type(o_vcf_2_hl7v2), vcf2hl7v2.Converter)


class TestTranslation(unittest.TestCase):
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

    def test_wo_patient_id(self):
        self.maxDiff = None
        o_vcf_2_hl7v2 = vcf2hl7v2.Converter(os.path.join(
            os.path.dirname(__file__), 'vcf_example1.vcf'), 'GRCh37',
            source_class='germline')
        output_filename = os.path.join(os.path.dirname(
            __file__), self.TEST_RESULT_DIR, 'hl7v2_wo_patient_example1.txt')
        expected_output_filename = os.path.join(
            os.path.dirname(__file__), 'expected_example1_wo_patient.txt')
        o_vcf_2_hl7v2.convert(output_filename)
        self.assertEqual(
                filecmp.cmp(output_filename, expected_output_filename),
                True
            )

    def test_with_patient_id(self):
        self.maxDiff = None
        o_vcf_2_hl7v2 = vcf2hl7v2.Converter(os.path.join(os.path.dirname(
            __file__), 'vcf_example1.vcf'), 'GRCh37', 'HG00628',
            source_class='germline')
        output_filename = os.path.join(os.path.dirname(
            __file__), self.TEST_RESULT_DIR, 'hl7v2_with_patient_example1.txt')
        expected_output_filename = os.path.join(os.path.dirname(
            __file__), 'expected_example1_with_patient.txt')
        o_vcf_2_hl7v2.convert(output_filename)
        self.assertEqual(
                filecmp.cmp(output_filename, expected_output_filename),
                True
            )

    def test_region_studied(self):
        self.maxDiff = None
        region_studied_filename = os.path.join(
            os.path.dirname(__file__), 'RegionsStudied_example3.bed')
        conv_region_filename = os.path.join(os.path.dirname(
            __file__), 'RegionsToConvert_example3.bed')
        output_filename = os.path.join(os.path.dirname(
            __file__), self.TEST_RESULT_DIR, 'hl7v2_example3.txt')
        expected_output_filename = os.path.join(
            os.path.dirname(__file__), 'expected_example3.txt')
        o_vcf_2_hl7v2 = vcf2hl7v2.Converter(
            os.path.join(
                os.path.dirname(__file__),
                'vcf_example3.vcf'),
            'GRCh38',
            'HG00628',
            conv_region_filename=conv_region_filename,
            region_studied_filename=region_studied_filename,
            source_class='germline')
        o_vcf_2_hl7v2.convert(output_filename)
        self.assertEqual(
                filecmp.cmp(output_filename, expected_output_filename),
                True
            )

    def test_region_studied_dict(self):
        conv_region_dict = {
            "Chromosome": ["X", "X", "M"],
            "Start": [50000, 55000, 50000],
            "End": [52000, 60600, 60025]
        }
        self.maxDiff = None
        region_studied_filename = os.path.join(
            os.path.dirname(__file__), 'RegionsStudied_example3.bed')
        output_filename = os.path.join(os.path.dirname(
            __file__), self.TEST_RESULT_DIR, 'hl7v2_example3_dict.txt')
        expected_output_filename = os.path.join(
            os.path.dirname(__file__), 'expected_example3.txt')
        o_vcf_2_hl7v2 = vcf2hl7v2.Converter(
            os.path.join(
                os.path.dirname(__file__),
                'vcf_example3.vcf'),
            'GRCh38',
            'HG00628',
            conv_region_dict=conv_region_dict,
            region_studied_filename=region_studied_filename,
            source_class='germline')
        o_vcf_2_hl7v2.convert(output_filename)
        self.assertEqual(
                filecmp.cmp(output_filename, expected_output_filename),
                True
            )

    # Check if region studied observation outside the vcf files are also
    # included in hl7v2 report.
    def test_multiple_region_studied(self):
        self.maxDiff = None
        region_studied_filename = os.path.join(
            os.path.dirname(__file__), 'RegionsStudied_example4.bed')
        conv_region_filename = os.path.join(os.path.dirname(
            __file__), 'RegionsToConvert_example4.bed')
        output_filename = os.path.join(os.path.dirname(
            __file__), self.TEST_RESULT_DIR, 'hl7v2_example4.txt')
        expected_output_filename = os.path.join(
            os.path.dirname(__file__), 'expected_example4.txt')
        o_vcf_2_hl7v2 = vcf2hl7v2.Converter(
            os.path.join(
                os.path.dirname(__file__),
                'vcf_example4.vcf'),
            'GRCh38',
            'HG00628',
            conv_region_filename=conv_region_filename,
            region_studied_filename=region_studied_filename,
            source_class='germline')
        o_vcf_2_hl7v2.convert(output_filename)
        self.assertEqual(
                filecmp.cmp(output_filename, expected_output_filename),
                True
            )

    def test_region_studied_only(self):
        region_studied_filename = os.path.join(
            os.path.dirname(__file__), 'RegionsStudied_example4.bed')
        output_filename = os.path.join(os.path.dirname(
            __file__), self.TEST_RESULT_DIR, 'hl7v2_example4_test.txt')
        o_vcf_2_hl7v2 = vcf2hl7v2.Converter(
            os.path.join(
                os.path.dirname(__file__),
                'vcf_example4.vcf'),
            'GRCh38',
            'HG00628',
            region_studied_filename=region_studied_filename,
            source_class='germline')
        o_vcf_2_hl7v2.convert(output_filename)

    def test_empty_hl7v2_txt(self):
        conv_region_filename = os.path.join(os.path.dirname(
            __file__), 'RegionsToConvert_empty_example4.bed')
        output_filename = os.path.join(os.path.dirname(
            __file__), self.TEST_RESULT_DIR, 'hl7v2_example4_test.txt')
        o_vcf_2_hl7v2 = vcf2hl7v2.Converter(
            os.path.join(
                os.path.dirname(__file__),
                'vcf_example4.vcf'),
            'GRCh38',
            'HG00628',
            conv_region_filename=conv_region_filename,
            source_class='germline')
        o_vcf_2_hl7v2.convert(output_filename)

    def test_tabix(self):
        self.maxDiff = None
        region_studied_filename = os.path.join(
            os.path.dirname(__file__), 'RegionsStudied_example4.bed')
        conv_region_filename = os.path.join(os.path.dirname(
            __file__), 'RegionsToConvert_example4.bed')
        output_filename = os.path.join(os.path.dirname(
            __file__), self.TEST_RESULT_DIR, 'hl7v2_example4_tabix.txt')
        expected_output_filename = os.path.join(
            os.path.dirname(__file__), 'expected_example4.txt')
        o_vcf_2_hl7v2 = vcf2hl7v2.Converter(
            os.path.join(
                os.path.dirname(__file__),
                'vcf_example4.vcf.gz'),
            'GRCh38',
            'HG00628',
            has_tabix=True,
            conv_region_filename=conv_region_filename,
            region_studied_filename=region_studied_filename,
            source_class='germline')
        o_vcf_2_hl7v2.convert(output_filename)
        self.assertEqual(
                filecmp.cmp(output_filename, expected_output_filename),
                True
            )

    def test_annotation(self):
        o_vcf_2_hl7v2 = vcf2hl7v2.Converter(
                filename=os.path.join(
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

    # def test_qci_xml(self):
    #     output_filename = os.path.join(os.path.dirname(
    #         __file__), self.TEST_RESULT_DIR, 'qci_hl7v2.txt')
    #     expected_output_filename = os.path.join(
    #         os.path.dirname(__file__), 'expected_qci_xml.txt')
    #     o_vcf_2_hl7v2 = vcf2hl7v2.Converter(
    #         filename=os.path.join(os.path.dirname(__file__), 'QCI.xml'),
    #         ref_build='GRCh37',
    #         patient_id='HG00628',
    #         source_class='somatic')
    #     o_vcf_2_hl7v2.convert(output_filename)
    #     self.assertEqual(
    #             filecmp.cmp(output_filename, expected_output_filename),
    #             True
    #         )

    # def test_qiagen_vcf(self):
    #     output_filename = os.path.join(os.path.dirname(
    #         __file__), self.TEST_RESULT_DIR, 'qiagen_vcf_hl7v2.txt')
    #     expected_output_filename = os.path.join(
    #         os.path.dirname(__file__), 'expected_qiagen_vcf.txt')
    #     o_vcf_2_hl7v2 = vcf2hl7v2.Converter(
    #         filename=os.path.join(os.path.dirname(__file__), 'qiagen.vcf'),
    #         ref_build='GRCh37',
    #         vcf_type='QiaGen',
    #         source_class='germline')
    #     o_vcf_2_hl7v2.convert(output_filename)
    #     self.assertEqual(
    #             filecmp.cmp(output_filename, expected_output_filename),
    #             True
    #         )


class TestLogger(unittest.TestCase):
    def setUp(self):
        # create file handler and set level to debug
        self.log_general_filename = os.path.join(
            os.path.dirname(__file__), self.LOG_DIR, 'general.log')
        self.log_invalid_record_filename = os.path.join(
            os.path.dirname(__file__), self.LOG_DIR, 'invalid_record.log')
        genearl_fh = logging.FileHandler(self.log_general_filename)
        invalid_record_fh = logging.FileHandler(
            self.log_invalid_record_filename)
        genearl_fh.setLevel(logging.DEBUG)
        invalid_record_fh.setLevel(logging.DEBUG)
        # create formatter
        formatter = logging.Formatter(
            '%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        # add formatter to fh
        genearl_fh.setFormatter(formatter)
        invalid_record_fh.setFormatter(formatter)
        # store as a class variable
        self.genearl_fh = genearl_fh
        self.invalid_record_fh = invalid_record_fh

    @classmethod
    def setUpClass(self):
        self.LOG_DIR = os.path.join(os.path.dirname(__file__), 'log')
        if os.path.exists(self.LOG_DIR):
            shutil.rmtree(self.LOG_DIR)
        os.mkdir(self.LOG_DIR)

    # TODO: Delete the log folder after running
    # all the tests, below method throws error becasue
    # log files are still in use.
    # @classmethod
    # def tearDownClass(self):
    #     shutil.rmtree(self.LOG_DIR)

    def test_logger_forks(self):
        region_studied_filename = os.path.join(
            os.path.dirname(__file__), 'RegionsStudied_example3.bed')
        conv_region_filename = os.path.join(os.path.dirname(
            __file__), 'RegionsToConvert_example3.bed')
        output_filename = os.path.join(os.path.dirname(
            __file__), self.LOG_DIR, 'logging_hl7v2.txt')
        # create logger
        general_logger = logging.getLogger('vcf2hl7v2.general')
        invalid_record_logger = logging.getLogger('vcf2hl7v2.invalidrecord')
        general_logger.setLevel(logging.DEBUG)
        invalid_record_logger.setLevel(logging.DEBUG)
        # add ch to logger
        general_logger.addHandler(self.genearl_fh)
        invalid_record_logger.addHandler(self.invalid_record_fh)
        o_vcf_2_hl7v2 = vcf2hl7v2.Converter(
            os.path.join(
                os.path.dirname(__file__),
                'vcf_example3.vcf'),
            'GRCh38',
            'HG00628',
            conv_region_filename=conv_region_filename,
            region_studied_filename=region_studied_filename,
            source_class='germline')
        o_vcf_2_hl7v2.convert(output_filename)
        self.assertEqual(os.path.exists(self.log_general_filename), True)
        self.assertEqual(os.path.exists(
            self.log_invalid_record_filename), True)


suite.addTests(
    unittest.TestLoader().loadTestsFromTestCase(Testvcf2hl7v2Inputs))
suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestTranslation))
suite.addTests(unittest.TestLoader().loadTestsFromTestCase(TestLogger))


if __name__ == '__main__':
    unittest.main()
