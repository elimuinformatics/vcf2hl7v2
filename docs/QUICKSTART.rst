.. _quickstart:

Quickstart
==========

Eager to get started? This page gives a good introduction in how to get started
with vcf2hl7v2.

Installation
---------------------

Installing vcf2hl7v2 is pretty simple. Here is a step by step plan on how to do it.

.. note::
    vcf2hl7v2 is available on TestPyPI as ``vcf2hl7v2``.

Before installing vcf2hl7v2 you need to install cython and wheel:

.. code-block:: bash
    
    pip install cython wheel  
    pip install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple/ vcf2hl7v2

This will install vcf2hl7v2 library.
 
Examples
---------------------

(some sample VCF files are `here <https://github.com/elimuinformatics/vcf2hl7v2/tree/master/vcf2hl7v2/test>`_)

Quick example

.. code-block:: python

    import vcf2hl7v2
    vcf_hl7v2_converter = vcf2hl7v2.Converter('sample.vcf', 'GRCh37')
    vcf_hl7v2_converter.convert()

More Example to instantiate Converter

-  Converts all variants in VCF. HL7 V2 message contains no region-studied
   OBX observation group.

.. code-block:: python

       vcf2hl7v2.Converter('vcftests.vcf','GRCh37', 'aabc')

-  Converts all variants in VCF. HL7 V2 message assign homoplasmic vs.
   heteroplasmic based on:

   If allelic depth (FORMAT.AD)/ read depth (FORMAT.DP) is greater than 0.89
   then allelic state is homoplasmic; else heteroplasmic.

   **Note** : default value of ratio_ad_dp = 0.99 and ratio_ad_dp is considered valid only when value lies between 0 and 1

.. code-block:: python

        vcf2hl7v2.Converter('vcftests.vcf','GRCh37', 'aabc', ratio_ad_dp = 0.89)

-  Converts all variants in VCF. HL7 V2 message contains one region-studied
   OBX observation group per studied chromosome.

.. code-block:: python

       vcf2hl7v2.Converter('vcftests.vcf','GRCh37', 'cabc', region_studied_filename='WGS_b37_region_studied.bed')

-  Converts all variants in VCF. HL7 V2 message contains one region-studied
   OBX observation group per studied chromosome.

.. code-block:: python

       vcf2hl7v2.Converter('vcftests.vcf','GRCh37', 'dabc', region_studied_filename='WGS_b37_region_studied.bed')

-  Converts all variants in conversion region. HL7 V2 message contains no
   region-studied OBX observation group.

.. code-block:: python

       vcf2hl7v2.Converter('vcftests.vcf','GRCh37', 'eabc', conv_region_filename='WGS_b37_convert_everything.bed')

-  Converts all variants in conversion region. HL7 V2 message contains one
   region-studied OBX observation group per studied chromosome, intersected with
   conversion region.

.. code-block:: python

       vcf2hl7v2.Converter('vcftests.vcf','GRCh37', 'habc', conv_region_filename='WGS_b37_convert_everything.bed', region_studied_filename='WGS_b37_region_studied.bed')

-  The conversion region will be ignored and all the variants in VCF for 
   which annotations are provided will be converted. HL7 V2 message contains one
   region-studied OBX observation group per studied chromosome.

.. code-block:: python

       vcf2hl7v2.Converter('NB6TK328_filtered.vcf','GRCh38', 'NB6TK328', conv_region_filename='NB6TK328_conversion_region.bed', region_studied_filename='NB6TK328_region_studied.bed',
       annotation_filename='NB6TK328_annotations.txt')

-  Conversion of a bgzipped VCF

.. code-block:: python

       vcf2hl7v2.Converter('vcf_example4.vcf.gz','GRCh37', 'kabc', has_tabix=True)

Logging
---------------------
You can use python standard `logging <https://docs.python.org/3/library/logging.html>`_ to enable logs. Two logger ('vcf2hl7v2.general') and ('vcf2hl7v2.invalidrecord') are avialble to configure.

-  **vcf2hl7v2.general**: standard library logs. 

-  **vcf2hl7v2.invalidrecord**: logs all the records from vcf file which are in conversion region but are not converted to HL7 V2 format.

.. code-block:: python

    >>> import logging
    # create logger
    >>> logger = logging.getLogger('vcf2hl7v2.invalidrecord')
    >>> logger.setLevel(logging.DEBUG)
    # create console handler and set level to debug
    >>> ch = logging.FileHandler('invalidrecord.log')
    >> ch.setLevel(logging.DEBUG)
    # create formatter
    >>> formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    # add formatter to ch
    >>> ch.setFormatter(formatter)
    # add ch to logger
    >>> logger.addHandler(ch)
