
import csv, sys, os
import logging
import ConfigParser

def parse_config(config_file):
    '''
    Reads in the configuration file, and copies it to a newly
    created directory in the output folder.
    '''
    config = ConfigParser.SafeConfigParser()
    config.read(config_file)
    return config

def parse_sample_info(sample_info_file, ipsc, somatic):
    ''' Parse sample_info file and return list of lists. '''
    sample_info_list = []
    try:
        with open(sample_info_file,'r') as f:
            reader = csv.reader(f,delimiter='\t',quoting=csv.QUOTE_NONE)
            for row in reader:
                if ipsc:
                    if len(row) == 4:
                        if not (len(sample_info_list) == 0 and row[0] == "line"):
                            sample_info_list.append(row)
                    else:
                        logging.error("For iPSC analysis sample_info file should have 4 columns:\tline\tclone\tpassage_number\tfile")
                        sys.exit(1)
                elif somatic:
                    if len(row) == 4 or len(row) == 5:
                        if not (len(sample_info_list) == 0 and row[0] == "sample"):
                            sample_info_list.append(row)
                    else:
                        logging.error("For Somatic analysis sample_info file should have 4 columns:\tsample\tsubject\tsample_type\tfile")
                        sys.exit(1)
                else:
                    if len(row) == 2 or len(row) == 3:
                        if not (len(sample_info_list) == 0 and row[0] == "sample"):
                            sample_info_list.append(row)
                    else:
                        logging.error("Sample_info file should have 2 columns:\tsample\tfile")
                        sys.exit(1)

                
    except:
        logging.error("Can't read file: " + sample_info_file)
        sys.exit(1)
    return sample_info_list


def parse_pedigree_file(pedigree_file):
    ''' Parse pedigree file and return dictionary of sampleID->familyID '''
    pedigree_family = {}
    if pedigree_file:
        try:
            with open(pedigree_file,'r') as f:
                reader = csv.reader(f,delimiter=' ',quoting=csv.QUOTE_NONE)
                for row in reader:
                    if len(row) == 6:
                        pedigree_family[row[1]] = row[0]
                    else:
                        logging.error("Pedigree file should have 6 columns with no header:\tFamilyID\tSampleID\tFatherID\tMotherID\tSex\tAffected")
                        sys.exit(1)   
        except:
            logging.error("Can't read file: " + pedigree_file)
            logging.error(sys.exc_info()[0])
            sys.exit(1)
    return pedigree_family
    
def parse_args(options, arguments):
    '''
    Parses out the arguments and options given from a terminal command and
    checks their validity. Returns a dictionary with these parameters.
    '''
    # Default Params
    params = {'output':None,'sample_info_file':None,
            'config_file':None, 'verbose':False,
            'parallel':False, 'ipsc':False, 'wgs':False,
            'pedigree_file':None, 'somatic':False,
            'analysis_types':0, 'group_analysis':False}
    for o, a in options:
        if o in ("-v", "--verbose"):
            params['verbose'] = True
        elif o in ("-h", "--help"):
            usage()
        elif o in ("-o", "--output"):
            params['output'] = a
        elif o in ("-s", "--sample_info"):
            params['sample_info_file'] = a
        elif o in ("-c", "--config"):
            params['config_file'] = a
        elif o in ("-p", "--parallel"):
            params['parallel'] = True
        elif o in ("-i", "--ipsc"):
            params['ipsc'] = True
            params['group_analysis'] = True
            params['analysis_types'] += 1
        elif o in ("-w", "--wgs"):
            params['wgs'] = True
        elif o in ("-f", "--pedigree"):
            params['pedigree_file'] = a
            params['group_analysis'] = True
            params['analysis_types'] += 1
        elif o in ("-t", "--somatic"):
            params['somatic'] = True
            params['group_analysis'] = True
            params['analysis_types'] += 1
        else:
            # no arguments given
            usage()

    # required parameters
    if not params['output'] or not params['sample_info_file'] or not params['config_file']:
        usage()
    if params['analysis_types'] > 1:
        logging.error("Cannot have multiple analysis types in the same run (somatic, pedigree, iPSC)!")
        usage()
    
    check_files_exist([params['sample_info_file'],
                       params['config_file']])
    
    return params


def check_files_exist(file_list):
    ''' Check if list of files exist '''
    for input_file in file_list:
        if not os.path.isfile(input_file):
            logging.error("Input file doesn't exist: " + input_file)
            sys.exit(1)

