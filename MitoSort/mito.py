'''

Main module for the MitoSort package

This module checks for appropriate packages, creates necessary directories,
and initializes samples for analysis.

Created on Aug 11, 2016

@author: Jared Evans evans.jared@mayo.edu 
'''
#from timeit import default_timer as timer
import sys, os, csv, subprocess
import logging
import getopt
import textwrap
import ConfigParser
import pysam
from shutil import copyfile, copytree, copy
from sample import Sample
import analysis
import parsing


def create_directories(params):
    ''' Create output directory structure. '''
    #directories = ("","alignment","coverage","coverage/wig","coverage/plot", "coverage/stats","docs","docs/config","docs/scripts","docs/logs","docs/html","variants","variants/vcf","variants/haplogroups","variants/annotation","variants/pileup","variants/plot") -- what's this for?
    directories = ("","alignment","coverage","coverage/wig","coverage/plot", "coverage/stats","docs","docs/config","docs/scripts","docs/logs","variants","variants/vcf","variants/haplogroups","variants/annotation","variants/pileup","variants/plot")
    output = params['output']
    for d in directories:
        full_path = output + "/" + d
        if not os.path.exists(full_path):
            try:
                os.mkdir(full_path)
            except:
                logging.error("Creating output directory: " + full_path)
                sys.exit(1)
        if full_path == output + "/":
            try:
                if os.listdir(full_path):
                    logging.error("Output directory is not empty!")
                    sys.exit(1)
            except:
                logging.error("checking contents of output dir.")
                sys.exit(1)

    copyfile(params['config_file'],
             params['output'] + "/docs/config/config.txt")
    copyfile(params['sample_info_file'],
             params['output'] + "/docs/config/sample_info.txt")
    
                
def env_check():
    ''' Check environment and tool versions. '''
    # may need to add to this
    if pysam.__version__ < 0.8:
        logging.error("You are using pysam " + pysam.__version__ + " and need to update to pysam 0.8+")
        sys.exit(1)

def usage():
    logging.error("Missing Required Parameters!")
    print textwrap.dedent('''
    EXAMPLE USAGE:
    python MitoSort.py -s sample_info.txt -o /path/to/output_dir/

    OPTIONS:
    -s,--sample_info        A tab-separated text file defining sample information (required)
    -o,--output             Full path to output directory (required)
    -c,--config             Config file (required)
    -p,--parallel           Parallel mode. Submits jobs to an SGE cluster using qsub. Default: False
    -i,--ipsc               Add iPSC annotations. Default: False
    -w,--wgs                Samples are from Whole Genome Sequencing. Default: False
    -f,--pedigree           Pedigree file for family analysis.
    -t,--somatic            Perform a Somatic analysis with Normal and Tumor samples. Default: False
    -h,--help               help message
    ''')
    sys.exit(2)


def create_sample_objects(sample_info_list, ipsc, somatic, pedigree):
    '''
    Builds and returns a list of sample objects from the sample info file.
    '''
    sample_objects = []
    
    for s in sample_info_list:
        sobj = Sample()
        sobj.load_sample_info(s, ipsc, somatic)
        if(pedigree):
            sobj.group = pedigree_info[sobj.sname]
            sobj.is_pedigree = True
            sobj.sample_type = "family"
        sample_objects.append(sobj)

    return sample_objects

def main():
    # set logging format
    #start = timer()
    logging.basicConfig(format='%(levelname)s\t%(asctime)s - %(message)s',
                        level=logging.INFO)
    logging.info("Starting Mitochondrial Analysis")
    
    if len(sys.argv) < 2:
        usage()
    try:
        opts, args = getopt.getopt(sys.argv[1:], "ho:s:c:piwbf:tv", 
                                   ["help", "output=", "sample_info=", 
                                    "config=", "parallel", "ipsc", 
                                    "wgs", "pedigree=", "somatic", 
                                    "verbose"])
    except getopt.GetoptError as err:
        print str(err)
        usage()
        sys.exit(2)
    
    params = parsing.parse_args(opts, args)

    create_directories(params)

    config = parsing.parse_config(params['config_file'])

    pedigree_info = parsing.parse_pedigree_file(params['pedigree_file'])

    sample_info = parsing.parse_sample_info(params['sample_info_file'], 
                                    params['ipsc'],
                                    params['somatic'])
    sample_objects = create_sample_objects(sample_info,
                                           params['ipsc'],
                                           params['somatic'],
                                           pedigree_info)
    output = params['output']
    wgs = params['wgs']
    group_analysis = params['group_analysis']
    ipsc = params['ipsc']

    if params['parallel']:
        # submit to cluster
        job_ids = []
        [job_ids.append(analysis.cluster_submit_sample_analysis(sobj, config, output, wgs)) for sobj in sample_objects]
                                            
        analysis.cluster_submit_allsample_reports(sample_objects, config, 
                                                  output, group_analysis, job_ids)
        logging.info("Finished Submitting Jobs for Mitochondrial Analysis")    
 
    else:
        # run locally
        [analysis.run_sample_analysis(sobj, config, output, wgs) for sobj in sample_objects]
    
        analysis.create_allsample_reports(sample_objects, config, output, 
                                          group_analysis, ipsc)

    logging.info("Finished Mitochondrial Analysis")
    #end = timer()
    #print end-start
if __name__ == '__main__':
    sys.exit(main())

