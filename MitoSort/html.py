'''
Functions to generate html report

Created on Dec 8, 2016

@author: Jared Evans evans.jared@mayo.edu
'''

import logging
import ConfigParser
import time
import sys
import csv

def write_html_report(sample_objects, config, output_dir, output_file, html_dir):
    '''
        write main html report
    '''
    html_dir = "docs/html/"
    
    sample_count = len(sample_objects)
    read_len = 0
    haplogroups = {}
    group_analysis = False
    ipsc = False
    for sobj in sample_objects:
        read_len = sobj.read_length
        haplogroups[sobj.haplogroup] = 1
        if sobj.is_ipsc or sobj.is_somatic or sobj.is_pedigree:
            group_analysis = True
        if sobj.is_ipsc:
            ipsc = True
    f = open(output_file,"w")
    f.write(html_head(html_dir, group_analysis, ipsc))
    haplogroup_count = len(haplogroups)
    min_het = "%.2f" % (float(config.get("PARAMS","HETEROPLASMY_FRACTION"))*100)
    f.write(html_summary(sample_count, haplogroup_count, read_len, min_het))
    if ipsc:
       f.write(html_ipsc_risk(sample_objects, output_dir))
    f.write(html_coverage_plots(sample_objects))
    f.write(html_coverage_stats(sample_objects, output_dir))
    f.write(html_haplogroups(sample_objects, output_dir))
    if group_analysis:
        f.write(html_private_mutation_plots(sample_objects))
    f.write(html_private_mutations(sample_objects, output_dir, ipsc))
    f.write(html_foot(html_dir))


def html_head(html_dir, group_analysis, ipsc):
    '''
        return head of html file
    '''
    html = '''<!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="utf-8">
        <meta http-equiv="X-UA-Compatible" content="IE=edge">
        <meta name="viewport" content="width=device-width, initial-scale=1">
        <!-- The above 3 meta tags *must* come first in the head; any other head content must come *after* these tags -->
        <meta name="description" content="">
        <meta name="author" content="">
        <!--<link rel="icon" href="favicon.ico">-->
    
        <title>MitoSort Mitochondrial Analysis Results</title>
    
        <!-- Bootstrap core CSS -->
        <link href="{0}css/bootstrap.min.css" rel="stylesheet">
    
        <!-- IE10 viewport hack for Surface/desktop Windows 8 bug -->
        <link href="{0}css/ie10-viewport-bug-workaround.css" rel="stylesheet">
    
        <!-- Custom styles for this template -->
        <link rel="stylesheet" type="text/css" href="docs/html/css/datatables.min.css"/>
        <link href="{0}mito.css" rel="stylesheet">
    
        <!-- HTML5 shim and Respond.js for IE8 support of HTML5 elements and media queries -->
        <!--[if lt IE 9]>
        <script src="https://oss.maxcdn.com/html5shiv/3.7.3/html5shiv.min.js"></script>
        <script src="https://oss.maxcdn.com/respond/1.4.2/respond.min.js"></script>
        <![endif]-->
    </head>
    
    <body data-spy="scroll" data-target="#myScrollspy" data-offset="20">
    
        <nav class="navbar navbar-inverse navbar-fixed-top">
            <div class="container-fluid">
                <div class="navbar-header">
                    <button type="button" class="navbar-toggle collapsed" data-toggle="collapse" data-target="#navbar" aria-expanded="false" aria-controls="navbar">
                        <span class="sr-only">Toggle navigation</span>
                        <span class="icon-bar"></span>
                        <span class="icon-bar"></span>
                        <span class="icon-bar"></span>
                    </button>
                    <a class="navbar-brand" href="#">MitoSort v0.3</a>
                </div>
                <div id="navbar" class="navbar-collapse collapse">
                    <ul class="nav navbar-nav navbar-right">
                        <li><a href="docs/MitoSort_help_manual.pdf">Help</a></li>
                    </ul>
                </div>
            </div>
        </nav>
        
        <div class="container-fluid">
            <div class="row">
                <div id="myScrollspy" class="col-sm-3 col-md-2 sidebar">
                    <ul class="nav nav-sidebar">
                        <li class="active"><a href="#summary">Summary <span class="sr-only">(current)</span></a></li>
    '''
    if ipsc:
        html += '''
                        <li><a href="#ipsc_risk">iPSC Risk</a></li>
        '''
    html += '''
                        <li><a href="#coverage_plots">Coverage Plots</a></li>
                        <li><a href="#coverage_stats">Coverage Stats</a></li>
                        <li><a href="#haplogroups">Haplogroups</a></li>
    '''
    if group_analysis:
        html += '''
                        <li><a href="#private_mutation_plots">Private Mutation Plots</a></li>
        '''
    html += '''
                        <li><a href="#private_mutations">Private Mutations</a></li>
                    </ul>
                </div>
            </div>
            <div class="col-sm-9 col-sm-offset-3 col-md-10 col-md-offset-2 main">
    '''    
    return html.format(html_dir)

def html_summary(sample_count, haplogroup_count, read_len, min_het):
    '''
        return summary section of html file
    '''
    date = time.strftime("%b %d, %Y")
    html = '''
                <h1 id="summary" class="page-header">Summary</h1>
                <div class="container-fluid">
                    <div class="col-sm-4">
                        <div class="panel panel-primary">
                            <div class="panel-heading">
                                <h3 class="panel-title">Sample Info</h3>
                            </div>
                            <div class="panel-body">
                                <table class="table">
                                    <tr>
                                        <td>Samples</td>
                                        <td>{0}</td>
                                    </tr>
                                    <tr>
                                        <td>Haplogroups</td>
                                        <td>{1}</td>
                                    </tr>
                                    <tr>
                                        <td>Read Length</td>
                                        <td>{2}</td>
                                    </tr>
                                </table>
                            </div>
                        </div>
                    </div>
                    <div class="col-sm-4">
                        <div class="panel panel-primary">
                            <div class="panel-heading">
                                <h3 class="panel-title">Analysis Info</h3>
                            </div>
                            <div class="panel-body">
                                <table class="table">
                                    <tr>
                                        <td>Date</td>
                                        <td>{3}</td>
                                    </tr>
                                    <tr>
                                        <td>Min Heteroplasmy Detection</td>
                                        <td>{4}%</td>
                                    </tr>
                                </table>
                            </div>
                        </div>
                    </div>
                    <div class="col-sm-4">
                        <div class="panel panel-primary">
                            <div class="panel-heading">
                                <h3 class="panel-title">Reference Info</h3>
                            </div>
                            <div class="panel-body">
                                <table class="table">
                                    <tr>
                                        <td>Sequence</td>
                                        <td>rCRS</td>
                                    </tr>
                                    <tr>
                                        <td>PhyloTree</td>
                                        <td>17</td>
                                    </tr>
                                </table>
                            </div>
                        </div>
                    </div>
                </div>
    '''
    return html.format(sample_count, haplogroup_count, read_len, date, min_het)

def html_ipsc_risk(sample_objects, output_dir):
    '''
        return html table for iPSC Risk scores
    '''
    html = '''
                <h2 id="ipsc_risk" class="sub-header">iPSC Risk</h2>
                <div class="container-fluid">
                    <div class="table-responsive">
                        <table id="ipsc_risk_table" class="table table-striped table-bordered" cellspacing="0" width="100%">
                            <thead>
                                <tr>
                                    <th>Line</th>
                                    <th>Clone</th>
                                    <th>Passage_Number</th>
                                    <th>Risk_Category</th>
                                    <th>Risk_Score</th>
                                    <th>Highest_Risk_Variant</th>
                                </tr>
                            </thead>
                            <tbody>
    '''

    html+= '''
                        <div class="col-sm-12">
                        <a href="iPSC_risk.png"><img src="iPSC_risk.png" class="img-responsive" alt="ipsc risk plot"></a>
                    </div>  
    '''

    for sobj in sample_objects:
        risk_row = []
        try:
            with open(output_dir+"/iPSC_Risk.csv",'r') as f:
                reader = csv.reader(f,delimiter=",",quoting=csv.QUOTE_MINIMAL)
                for row in reader:
                    if row[0] == sobj.group and row[1] == sobj.clone and row[2] == sobj.passage_number:
                        risk_row = row
        except:
            logging.error("Reading file " + output_dir+"/iPSC_Quality.csv")
            logging.error(sys.exc_info()[0])
            sys.exit(1)
            
        html_table = '''
                                <tr>
                                    <td>{0}</td>
                                    <td>{1}</td>
                                    <td>{2}</td>
                                    <td>{3}</td>
                                    <td>{4}</td>
                                    <td>{5}</td>
                                </tr>
        '''
        html += html_table.format(risk_row[0],risk_row[1],risk_row[2],risk_row[3],risk_row[4],risk_row[5])
    
    html += '''
                            </tbody>
                        </table>
                    </div>
                </div>
    '''
    return html

def html_coverage_plots(sample_objects):
    '''
        return coverage plots html
    '''
    html = '''
                <h2 id="coverage_plots" class="sub-header">Coverage Plots</h2>
                <div class="container-fluid">    
    '''
    for sobj in sample_objects:
            html_img = '''
                    <div class="col-sm-6">
                        <a href="coverage/plot/{0}.coverage.png"><img src="coverage/plot/{0}.coverage.png" width="500" class="img-responsive" alt="coverage plot"></a>
                    </div>
            '''
            html += html_img.format(sobj.sname)
                     
    html += '''
                </div>
    '''
    return html

def html_coverage_stats(sample_objects, output_dir):
    '''
        return coverage stats html
    '''
    html = '''
                <h2 id="coverage_stats" class="sub-header">Coverage Stats</h2>
                <div class="container-fluid">
                    <div class="table-responsive">
                        <table id="cov_stats_table" class="table table-striped table-bordered" cellspacing="0" width="100%">
                            <thead>
                                <tr>
                                    <th>Sample</th>
                                    <th>Mito_Reads</th>
                                    <th>Nucl_Reads</th>
                                    <th>Mito_Copies</th>
                                    <th>Average</th>
                                    <th>Minimum</th>
                                    <th>First_Quartile</th>
                                    <th>Median</th>
                                    <th>Third_Quartile</th>
                                    <th>Maximum</th>
                                    <th>No_Coverage</th>
                                </tr>
                            </thead>
                            <tbody>
    '''
    for sobj in sample_objects:
        stats_row = []
        try:
            with open(output_dir+"/coverage/stats/"+sobj.sname+".coverage_stats.csv",'r') as f:
                reader = csv.reader(f,delimiter=",",quoting=csv.QUOTE_MINIMAL)
                for row in reader:
                    if row[0] == sobj.sname:
                        stats_row = row
        except:
            logging.error("Reading file " + output_dir+"/coverage/stats/"+sobj.sname+".coverage_stats.csv")
            logging.error(sys.exc_info()[0])
            sys.exit(1)
            
        html_table = '''
                                <tr>
                                    <td>{0}</td>
                                    <td>{1}</td>
                                    <td>{2}</td>
                                    <td>{3}</td>
                                    <td>{4}</td>
                                    <td>{5}</td>
                                    <td>{6}</td>
                                    <td>{7}</td>
                                    <td>{8}</td>
                                    <td>{9}</td>
                                    <td>{10}</td>
                                </tr>
        '''
        #html += html_table.format(stats_row)
        html += html_table.format(stats_row[0],stats_row[1],stats_row[2],stats_row[3],stats_row[4],stats_row[5],stats_row[6],stats_row[7],stats_row[8],stats_row[9],stats_row[10])
    
    html += '''
                            </tbody>
                        </table>
                    </div>
                </div>
    '''
    return html

def html_haplogroups(sample_objects, output_dir):
    '''
        return haplogroup html
    '''
    html = '''
                <h2 id="haplogroups" class="sub-header">Haplogroups</h2>
                <div class="container-fluid">
                    <div class="table-responsive">
                        <table id="haplogroups_table" class="table table-striped table-bordered" cellspacing="0" width="100%">
                            <thead>
                                <tr>
                                    <th>Sample</th>
                                    <th>Haplogroup</th>
                                    <th>Quality_Score</th>
                                    <th>Polymorphisms_Found</th>
                                    <th>Polymorphisms_Missing</th>
                                    <th>Global_Private_Mutations</th>
                                    <th>Local_Private_Mutations</th>
                                </tr>
                            </thead>
                            <tbody>
    '''
    for sobj in sample_objects:
        stats_row = []
        try:
            with open(output_dir+"/variants/haplogroups/"+sobj.sname+".haplogroup_summary.csv",'r') as f:
                reader = csv.reader(f,delimiter=",",quoting=csv.QUOTE_MINIMAL)
                for row in reader:
                    if row[0] == sobj.sname:
                        stats_row = row
        except:
            logging.error("Reading file " + output_dir+"/variants/haplogroups/"+sobj.sname+".haplogroup_summary.csv")
            logging.error(sys.exc_info()[0])
            sys.exit(1)
            
        html_table = '''
                                <tr>
                                    <td>{0}</td>
                                    <td>{1}</td>
                                    <td>{2}</td>
                                    <td>{3}</td>
                                    <td>{4}</td>
                                    <td>{5}</td>
                                    <td>{6}</td>
                                </tr>
        '''
        html += html_table.format(stats_row[0],stats_row[1],stats_row[2],stats_row[3],stats_row[4],stats_row[5],stats_row[6])
        
    html += '''
                            </tbody>
                        </table>
                    </div>
                </div>
    '''
    return html

def html_private_mutation_plots(sample_objects):
    '''
        return private mut plot html
    '''
    html = '''
                <h2 id="private_mutation_plots" class="sub-header">Private Mutation Plots</h2>
                <div class="container-fluid">
    '''
    groups = {}
    for sobj in sample_objects:
        groups[sobj.group] = 1
    for group in groups.keys():
        html_img = '''
                    <div class="col-sm-12">
                        <a href="variants/plot/{0}.private_mutations.png"><img src="variants/plot/{0}.private_mutations.png" class="img-responsive" alt="private mutation plot"></a>
                    </div>  

        '''
        html += html_img.format(group)
    html += '''
                </div>
    '''
    return html

def html_private_mutations(sample_objects, output_dir, ipsc):
    '''
        return private mut html table
    '''
    html = '''
                <h2 id="private_mutations" class="sub-header">Private Mutations</h2>
                <div class="container-fluid">
            
                            <div class="table-responsive">
                                <table id="priv_mut_table" class="table table-striped table-bordered" cellspacing="0" width="100%">
                                    <thead>
                                        <tr>
                                            <th>Sample</th>
                                            <th>Pos</th>
                                            <th>Ref</th>
                                            <th>Alt</th>
                                            <th>Ref_Reads</th>
                                            <th>Alt_Reads</th>
                                            <th>Total_Depth</th>
                                            <th>Alt_Percent</th>
                                            <th>Genotype</th>
                                            <th>Variant_Type</th>
                                            <th>Private_Category</th>
                                            <th>Locus</th>
                                            <th>MITOMAP_Poly</th>
                                            <th>Codon_Change</th>
                                            <th>AA_Change</th>
                                            <th>AA_Pos</th>
                                            <th>Disease</th>
                                            <th>GenBank_AF</th>
                                            <th>PolyPhen2_Pred</th>
                                            <th>PolyPhen2_Score</th>
                                            <th>SIFT_Pred</th>
                                            <th>SIFT_Score</th>
                                            <th>PROVEAN_Pred</th>
                                            <th>PROVEAN_Score</th>
                                            <th>MutationAssessor_Pred</th>
                                            <th>MutationAssessor_Score</th>
                                            <th>CADD_Score</th>
                                            <th>RNA_Structure</th>
                                        </tr>
                                    </thead>
                                    <tbody>
    '''
    html_table = '''
                                        <tr>
                                            <td>{0}</td>
                                            <td>{1}</td>
                                            <td>{2}</td>
                                            <td>{3}</td>
                                            <td>{4}</td>
                                            <td>{5}</td>
                                            <td>{6}</td>
                                            <td>{7}</td>
                                            <td>{8}</td>
                                            <td>{9}</td>
                                            <td>{10}</td>
                                            <td>{11}</td>
                                            <td>{12}</td>
                                            <td>{13}</td>
                                            <td>{14}</td>
                                            <td>{15}</td>
                                            <td>{16}</td>
                                            <td>{17}</td>
                                            <td>{18}</td>
                                            <td>{19}</td>
                                            <td>{20}</td>
                                            <td>{21}</td>
                                            <td>{22}</td>
                                            <td>{23}</td>
                                            <td>{24}</td>
                                            <td>{25}</td>
                                            <td>{26}</td>
                                            <td>{27}</td>
                                        </tr>
    '''
    for sobj in sample_objects:
        try:
            with open(output_dir+"/variants/annotation/"+sobj.sname+".private_mutations.csv",'r') as f:
                reader = csv.reader(f,delimiter=",",quoting=csv.QUOTE_MINIMAL)
                for row in reader:
                    if ipsc:
                        if row[0]+"."+row[1]+".p"+row[2] == sobj.sname:
                            html += html_table.format(sobj.sname, row[6], row[7], row[8], row[9], row[10], row[11], row[12], row[13], row[14], row[5], row[15], row[16], row[17], row[18], row[19], row[20], row[21], row[22], row[23], row[24], row[25], row[26], row[27], row[28], row[29], row[30], row[31])
                    else:
                        if row[0] == sobj.sname:
                            html += html_table.format(sobj.sname, row[4], row[5], row[6], row[7], row[8], row[9], row[10], row[11], row[12], row[3], row[13], row[14], row[15], row[16], row[17], row[18], row[19], row[20], row[21], row[22], row[23], row[24], row[25], row[26], row[27], row[28], row[29])
                    
        except:
            logging.error("Reading file " + output_dir+"/variants/annotation/"+sobj.sname+".private_mutations.csv")
            logging.error(sys.exc_info()[0])
            sys.exit(1)
    html += '''
                                    </tbody>
                                </table>
                            </div>
                </div>
    '''
    return html

def html_all_variants(sample_objects, output_dir, ipsc):
    '''
        return all variant html table
    '''
    html = '''
                <h2 id="all_variants" class="sub-header">All Variants</h2>
                <div class="container-fluid">
                        
                            <div class="table-responsive">
                                <table id="all_var_table" class="table table-striped table-bordered" cellspacing="0" width="100%">
                                    <thead>
                                        <tr>
                                            <th>Sample</th>
                                            <th>Haplogrep_Cat</th>
                                            <th>Pos</th>
                                            <th>Ref</th>
                                            <th>Alt</th>
                                            <th>Ref_Reads</th>
                                            <th>Alt_Reads</th>
                                            <th>Total_Depth</th>
                                            <th>Alt_Percent</th>
                                            <th>Genotype</th>
                                            <th>Variant_Type</th>
                                            <th>Locus</th>
                                            <th>MITOMAP_Poly</th>
                                            <th>GenBank_AF</th>
                                        </tr>
                                    </thead>
                                    <tbody>
    '''
    html_table = '''
                                        <tr>
                                            <td>{0}</td>
                                            <td>{1}</td>
                                            <td>{2}</td>
                                            <td>{3}</td>
                                            <td>{4}</td>
                                            <td>{5}</td>
                                            <td>{6}</td>
                                            <td>{7}</td>
                                            <td>{8}</td>
                                            <td>{9}</td>
                                            <td>{10}</td>
                                            <td>{11}</td>
                                            <td>{12}</td>
                                            <td>{13}</td>
                                        </tr>
    '''
    for sobj in sample_objects:
        try:
            with open(output_dir+"/variants/annotation/"+sobj.sname+".all_variants.csv",'r') as f:
                reader = csv.reader(f,delimiter=",",quoting=csv.QUOTE_MINIMAL)
                for row in reader:
                    if ipsc:
                        if row[0]+"."+row[1]+".p"+row[2] == sobj.sname:
                            html += html_table.format(sobj.sname, row[5], row[6], row[7], row[8], row[9], row[10], row[11], row[12], row[13], row[14], row[15], row[16], row[17])
                    else:
                        if row[0] == sobj.sname:
                            html += html_table.format(sobj.sname, row[3], row[4], row[5], row[6], row[7], row[8], row[9], row[10], row[11], row[12], row[13], row[14], row[15])
                    
        except:
            logging.error("Reading file " + output_dir+"/variants/annotation/"+sobj.sname+".all_variants.csv")
            logging.error(sys.exc_info()[0])
            sys.exit(1)
    html += '''
                                    </tbody>
                                </table>
                            </div>
                </div>
    '''
    return html

def html_foot(html_dir):
    '''
        return the foot of the html 
    '''
    html = '''
            </div>
        </div>
    
        <!-- Bootstrap core JavaScript
        ================================================== -->
        <!-- Placed at the end of the document so the pages load faster -->
        <script src="https://ajax.googleapis.com/ajax/libs/jquery/1.12.4/jquery.min.js"></script>
        <script>window.jQuery || document.write('<script src="{0}js/jquery.min.js"><\/script>')</script>
        <script src="{0}js/bootstrap.min.js"></script>
        <script type="text/javascript" src="docs/html/js/datatables.min.js"></script>
        <script type="text/javascript" src="docs/html/mito.js"></script>
        <!-- IE10 viewport hack for Surface/desktop Windows 8 bug -->
        <script src="{0}js/ie10-viewport-bug-workaround.js"></script>
    </body>
    </html>
    '''
    return html.format(html_dir)
