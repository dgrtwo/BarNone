.. Basic usage of BarNone script

Usage
===================================

BarNone is usually used directly as a script, with the following options:

.. program-output:: BarNone -h

..  Options
    -------

    Read Composition
    ~~~~~~~~~~~~~~~~
    
    All within-read positions are indexed starting at 1.
    
    -s START, --start START
        Start position of barcode
    
    -l LENGTH, --length LENGTH
        Length of barcode. If shorter barcodes exist in the catalog, those substrings (starting from the beginning of the read) will also be checked.
    
    --uptag UPTAG
        Sequence designating an up-tag.
    
    --downtag 

Input
-----

Two files are used for input: one containing the sequencing results, and another containing the barcode catalog file.

The sequencing file can be in any of the following formats:

* `fastq <http://maq.sourceforge.net/fastq.shtml>`_
* `fasta <http://blast.ncbi.nlm.nih.gov/blastcgihelp.shtml>`_
* `qseq <http://jumpgate.caltech.edu/wiki/QSeq>`_
* Raw text: one read per line.

The barcode catalog file should be tab-delimited, with content analogous to::

    YBL103C	AGTCTACCCACATGCTTTAG	ATTCATAGGACACTTGCCGG
    YBL104C	ACTCATCAGGTATGGGAACG	AATTCTAGCGCATGGTCCGG
    YBL105C	AGTGCTAATACGGCATAACG	ACACTACCTAGATGGAGCGG

Each line has three components: a unique identifier for the gene, an uptag sequence, and a downtag sequence. The uptag and downtag sequences are each optional; if one is missing, it will be counted as having 0 reads in the output.

Output
------

Counts
~~~~~~

BarNone produces a tab-delimited file describing the number of reads that matched to each barcode.

If the file was not multiplexed, the output counts file has three columns, such as::

    YAL008W	352	0
    YBR255W	505	1411
    YGR164W	429	1255

This describes, for each strain, how many up and down barcodes (respectively) were found.

If the file was multiplexed using the :command:`--multiplexfile` option, the counts file contains two columns for each multiplexing barcode, and a header row describing each. For example::

    Strain	Sample1|UP	Sample1|DOWN	Sample2|UP	Sample2|DOWN	Sample3|UP	Sample3|DOWN
    YAL008W	352	0	471	0	478	0
    YBR255W	505	1411	724	1319	754	1476
    YGR164W	429	1255	608	1594	596	1495

Revised Barcode Catalog 
~~~~~~~~~~~~~~~~~~~~~~~

BarNone is capable of constructing a new barcode catalog based on common mismatches found in a sequencing file. If a mismatched version of a barcode is found more frequently than the original, BarNone replaces that barcode with the more common version. This is specified using the :command:`--revisedcatalog` option.

Mismatches
~~~~~~~~~~

One can also obtain a much more detailed report about the specific mismatches using the :command:`--mismatchfile` option. The specified file will be written in a format analogous to::

    YOR041C	ACATCCCGATCAGGTGACTG	23215	ACATCCCGATCAGGCGACTG (62)/CCATCCCGATCAGGTGACTG (34)/ACATCCCGATCAGNTGACTG (33)
    YEL008W	GCATCACAGCAATGGGCATA	2989	GCATCACAGCACTGGGCATA (20)/GCATCACAGCAACGGGCATA (8)/GCATCACAGCAATGGGCATC (8)
    YMR011W	GCCTCACTTAAAGCATACGA	6140	GCCTCACTTAACGCATACGA (33)/ACCTCACTTAAAGCATACGA (14)/GCATCACTTAAAGCATACGA (11)
    YDR440W	GCGGCCCTACAATTTATGAA	17	CCGGCCCTACAATTTATGAA (27840)/CCGGCCCTACACTTTATGAA (52)/ACGGCCCTACAATTTATGAA (46)

Each line contains the name of the strain, the original barcode, the number of times the original barcode appeared in the reads, and then a slash-delimited list of the close barcodes that were matched with the number of times each of them occured.

Of these four barcodes, YDR440W has a close match that is much more common than its original barcode, meaning that it would be replaced in a revised catalog.

Examples
--------

Nutrient Starvation Data
~~~~~~~~~~~~~~~~~~~~~~~~

Using the Gresham et al 2011 nutrient starvation data available `here <http://genomics-pubs.princeton.edu/StarvationGenetics/download.shtml>`_ (in this example, the file `Expt10.txt <ftp://gen-ftp.princeton.edu/StarvationGenetics/Expt10.txt.bz2>`_), would require a command like::

    BarNone -f fastq -m 2 Expt10.txt counts.txt barcode_catalog.txt

This would match all barcodes in the file to a matching read in the catalog within an edit distance of 2 (set by the ``-m`` option).

An example of such a barcode catalog file is available in the :file:`examples/smith_et_al_2009_barcodes.txt`.
    