from mito import *
import logging 
import parsing
import analysis
import subprocess
import ConfigParser

output_dir = '/data2/bsi/tertiary/Nelson_Timothy_m000917/s203481.mtDNA_tool_development/full_output'
config = parsing.parse_config('/home/m169420/config.txt')
job_ids = []
shell_script = output_dir + "/docs/scripts/mito_analysis.allsample_reports.sh"

qsub = config.get("CLUSTER","QSUB")
queue = config.get("CLUSTER","QUEUE")
memory = config.get("CLUSTER","MEMORY")
email = config.get("CLUSTER","EMAIL")

s = ","
#hold_jobs = s.join(job_ids)

#qsub_cmd = [qsub, "-wd",output_dir+"/docs/logs/","-q",queue,"-hold_jid",hold_jobs,"-M",email,"-m","a","-l",memory,"-l","h_stack=10M",shell_script]
qsub_cmd = [qsub, "-wd",output_dir+"/docs/logs/","-q",queue,"-M",email,"-m","a","-l",memory,"-l","h_stack=10M",shell_script]
qsub_output = subprocess.check_output(qsub_cmd)
print qsub_output.strip()
