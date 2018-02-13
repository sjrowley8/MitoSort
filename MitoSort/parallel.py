'''
Functions for running the analysis in parallel on an SGE cluster

Created on Dec 8, 2016

@author: Jared Evans evans.jared@mayo.edu
'''

import os, sys, subprocess
import ConfigParser
import logging
import pickle
from sample import Sample
from analysis import run_sample_analysis
from analysis import create_allsample_reports
from mito import parse_sample_info




def cluster_submit_sample_analysis(sobj, config, output_dir, wgs):
    '''
    Submit the sample analysis job and return job id
    '''
    # write a shell script
    mitosort_path = config.get("TOOLS","MITOSORT")
    python = config.get("TOOLS","PYTHON")
    
    shell_script = output_dir + "/docs/scripts/mito_analysis." + sobj.sname + ".sh"
    f = open(shell_script,"w")
    f.write("#!/bin/sh\n\n")
    f.write("echo $(date)\n")
    f.write("PYTHONPATH=" + mitosort_path + ":$PYTHONPATH \n")
    f.write(python + " -c 'import parallel; parallel.cluster_run_sample_analysis(\"" + sobj.sname + "\",\"" + output_dir + "\",\"" + str(sobj.is_ipsc) + "\",\"" + str(sobj.is_somatic) + "\",\"" + str(wgs) + "\") '\n")
    f.write("echo $(date)\n")
    f.close()

    qsub = config.get("CLUSTER","QSUB")
    queue = config.get("CLUSTER","QUEUE")
    memory = config.get("CLUSTER","MEMORY")
    email = config.get("CLUSTER","EMAIL")

    qsub_cmd = [qsub,"-V","-wd",output_dir+"/docs/logs/","-q",queue,"-M",email,"-m","a","-l",memory,"-l","h_stack=10M",shell_script]
    qsub_output = subprocess.check_output(qsub_cmd)
    print qsub_output.strip()
    job_id = qsub_output.split(" ")[2]

    return job_id


def cluster_run_sample_analysis(sample_name, output_dir, ipsc, somatic, wgs, sobjs):
    '''
    Run the analysis for a sample, to be called from the SGE job
    '''
    # load necessary objects
    logging.basicConfig(format='%(levelname)s\t%(asctime)s - %(message)s', level=logging.INFO)
    config_file = output_dir + "/docs/config/config.txt"
    sample_info_file = output_dir + "/docs/config/sample_info.txt"
    config = ConfigParser.SafeConfigParser()
    config.read(config_file)
    
    if ipsc == "True":
        ipsc = True
    else:
        ipsc = False
        
    if somatic == "True":
        somatic = True
    else:
        somatic = False
        
    if wgs == "True":
        wgs = True
    else:
        wgs = False
    
    sample_info = parse_sample_info(sample_info_file, ipsc, somatic)
    this_sobj = None
    for s in sample_info:
        sobj = Sample()
        sobj.load_sample_info(s, ipsc, somatic)
        if sobj.sname == sample_name:
            this_sobj = sobj
            break
    
    # run usual analysis from mito module
    run_sample_analysis(this_sobj, config, output_dir, wgs)
    # dump necessary objects to disk for use in the next job
    #this_sobj.dump_variants_to_file(output_dir)
    output_stream = open(output_dir + "/docs/scripts/" + this_sobj.sname + ".sample.pkl", 'wb')
    pickle.dump(this_sobj, output_stream)
    output_stream.close()
    
def cluster_submit_allsample_reports(sample_objects, config, output_dir, group_analysis, job_ids):
    '''
    submit the allsample reports SGE job
    ''' 
    # write a shell script
    mitosort_path = config.get("TOOLS","MITOSORT")
    python = config.get("TOOLS","PYTHON")
    
    ipsc = False
    somatic = False
    pedigree = False
    if group_analysis:
        for sobj in sample_objects:
            if sobj.is_ipsc:
                ipsc = True
            if sobj.is_somatic:
                somatic = True
            if sobj.is_pedigree:
                pedigree = True
    
    shell_script = output_dir + "/docs/scripts/mito_analysis.allsample_reports.sh"
    f = open(shell_script,"w")
    f.write("#!/bin/sh\n\n")
    f.write("echo $(date)\n")
    f.write("PYTHONPATH=" + mitosort_path + ":$PYTHONPATH \n")
    f.write(python + " -c 'import parallel; parallel.cluster_run_allsample_reports(\"" + output_dir + "\",\"" + str(ipsc) + "\",\"" + str(somatic) + "\",\"" + str(pedigree) + "\") '\n")
    f.write("echo $(date)\n")
    f.close()

    qsub = config.get("CLUSTER","QSUB")
    queue = config.get("CLUSTER","QUEUE")
    memory = config.get("CLUSTER","MEMORY")
    email = config.get("CLUSTER","EMAIL")

    s = ","
    hold_jobs = s.join(job_ids)

    qsub_cmd = [qsub,"-V","-wd",output_dir+"/docs/logs/","-q",queue,"-hold_jid",hold_jobs,"-M",email,"-m","a","-l",memory,"-l","h_stack=10M",shell_script]
    qsub_output = subprocess.check_output(qsub_cmd)
    print qsub_output.strip()
    
    
def cluster_run_allsample_reports(output_dir, ipsc, somatic, pedigree):
    '''
    run the allsample reports code, called from an SGE job
    '''
    # load necessary objects
    logging.basicConfig(format='%(levelname)s\t%(asctime)s - %(message)s', level=logging.INFO)
    config_file = output_dir + "/docs/config/config.txt"
    sample_info_file = output_dir + "/docs/config/sample_info.txt"
    config = ConfigParser.SafeConfigParser()
    config.read(config_file)
    
    group_analysis = False
    
    if ipsc == "True":
        ipsc = True
        group_analysis = True
    else:
        ipsc = False
        
    if somatic == "True":
        somatic = True
        group_analysis = True
    else:
        somatic = False
        
    if pedigree == "True":
        pedigree = True
        group_analysis = True
    else:
        pedigree = False
    
    sample_info = parse_sample_info(sample_info_file, ipsc, somatic)
    sample_objects = []
    for s in sample_info:
        sobj = Sample()
        sobj.load_sample_info(s, ipsc, somatic)
        # load object data from pickle dump
        #sobj.load_variants_from_file(output_dir)
        output_stream = open(output_dir + "/docs/scripts/" + sobj.sname + ".sample.pkl", 'rb')
        sobj = pickle.load(output_stream)
        output_stream.close()
        os.remove(output_dir + "/docs/scripts/" + sobj.sname + ".sample.pkl")
        sample_objects.append(sobj)
    
    # run usual analysis from mito module
    create_allsample_reports(sample_objects, config, output_dir, group_analysis, ipsc)
    
