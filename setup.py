from setuptools import setup
import os, sys, subprocess

def build_config():
    dir_path = os.path.realpath('') 
    mitosort = dir_path + '/MitoSort/'
    config = dir_path + '/config/'
    if not os.path.exists(config):
        os.makedirs(config)

    # MAKE SURE DEPENDENCIES ARE INSTALLED
    bwa = subprocess.check_output(['which', 'bwa'])
    samtools = subprocess.check_output(['which', 'samtools'])
    qsub = subprocess.check_output(['which', 'qsub'])
    
    f = open(config + 'config.txt', 'w')
    f.write('[PARAMS]\n')
    f.write('# PARAMETERS\nHETEROPLASMY_FRACTION = 0.01\nMAPPING_QUALITY = 20\nBASE_QUALITY = 20\n')
    f.write('STRAND_BIAS_FILTER = 5\nLOG_SCALE_COVERAGE_PLOT = True\nMIN_HET_FOR_HAPLOGROUPING = 0.5\n')
    f.write('MIN_HET_TO_BACKFILL = 0.10\nRCRS_LEN = 16569\nCOVERAGE_RESOLUTION = 21\nMAX_PILEUP_DEPTH = 1000000\n\n')
    f.write('[REF]\n')
    f.write('# REFERENCES\nMITO_REF = %(MITOSORT)s/resources/rCRS.fa\n')
    f.write('GENOME_REF = %(MITOSORT)s/resources/whole_genome/links/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa\n')
    f.write('MITO_GENES = %(MITOSORT)s/resources/mtdna_genes.txt\n')
    f.write('MITO_GENES_PLOT = %(MITOSORT)s/resources/mtdna_genes.plot.txt\n')
    # MITO_CODING_GENES IS NOT USED / DOES NOT EXIST IN RESOURCES
    f.write('MITO_CODING_GENES = %(MITOSORT)s/resources/mtdna_coding_genes.txt\n')
    # MITO_CODON_TABLE IS NOT USED / DOES NOT EXIST IN RESOURCES 
    f.write('MITO_CODON_TABLE = %(MITOSORT)s/resources/mtdna_codon_table.txt\n')
    f.write('MITOMAP_COUNTS = %(MITOSORT)s/resources/Mitomap_Allele_Counts_2016_June25.csv\n')
    f.write('MITOMAP_POLY_CODING = %(MITOSORT)s/resources/PolymorphismsCoding_lt_MITOMAP_lt_Foswiki.csv\n')
    f.write('MITOMAP_POLY_CONTROL = %(MITOSORT)s/resources/PolymorphismsControl_lt_MITOMAP_lt_Foswiki.csv\n')
    f.write('MITOMAP_MUT_CODING_CONTROL = %(MITOSORT)s/recoures/MutationsCodingControl_lt_MITOMAP_lt_Foswiki.csv\n')
    f.write('MITOMAP_MUT_RNA = %(MITOSORT)s/resources/MutationsRNA_lt_MITOMAP_lt_Foswiki.csv\n')
    f.write('MITOMAP_MUT_SOMATIC = %(MITOSORT)s/resources/MutationsSomatic_lt_MITOMAP_lt_Foswiki.csv\n')
    f.write('MITIMPACT = %(MITOSORT)s/resources/MitImpact_db_2.5.select_columns.txt\n')
    f.write('PHYLOTREE = %(MITOSORT)s/resources/phylotree/mtDNA_tree_Build_17-rCRS_oriented_version.tab.reformatted.txt\n')
    f.write('MITOSORT = ' + mitosort + '\n\n')
    f.write('[TOOLS]\n# TOOL PATHS\n')
    # mitosort in config 2x?
    f.write('MITOSORT = ' + mitosort + '\n')
    f.write('PYTHON = ' + sys.executable + '\n')
    f.write('MATPLOTLIB = ' + '/'.join(matplotlib.__file__.split('/')[:-2]) + '\n')
    f.write('PYSAM = ' + '/'.join(pysam.__file__.split('/')[:-2]) + '\n')
    f.write('VIENNARNA = ' + '/'.join(RNA.__file__.split('/')[:-2]) + '\n')
    f.write('BWA = ' + bwa)
    f.write('SAMTOOLS = ' + samtools)
    #java is not written in here
    f.write('HTML_FILES = %(MITOSORT)s/html\n\n')
    f.write('[CLUSTER]\n# CLUSTER PARAMETERS\n')
    f.write('QSUB = ' + qsub)
    f.write('QUEUE = 1-day\nMEMORY= h_vmem=8G\nEMAIL = ENTER EMAIL HERE\n')
    
if __name__ == '__main__':
    
    
    DESCRIPTION = 'Mitosort is a python tool for mitochondrial DNA analysis.'  #with open('README.md') as f:
    #     LONG_DESCRIPTION = f.read()

 
    setup(
        name='mitosort',
        version='1.0',
        description=DESCRIPTION,
        long_description='',
        license='', # LOOK INTO THIS
        author='Jared Evans',
        author_email=' ',
        packages=['MitoSort'],
        install_requires=['pandas', 'numpy', 
                          'matplotlib', 'pysam']
        )

    build_config()
