####### Trigger a luigi pipeline #######

import ConfigParser
import sys
import subprocess

def main(samples_config):
    ''' The main function
    '''

    parser = ConfigParser.ConfigParser()
    parser.read(samples_config)

    for section in parser.sections():
        sample_name = section
        R1_fastq = parser.get(section,'R1_fastq')
        R2_fastq = parser.get(section,'R2_fastq')

        cmd = """PYTHONPATH=$PYTHONPATH:"" luigi --module single_cell_rnaseq JoinCountFiles --R1-fastq {R1} --R2-fastq {R2} --sample-name {sample} --workers 14 --worker-wait-interval 20""".format(R1=R1_fastq,R2=R2_fastq,sample=sample_name)
        print "Running command : {cmd}\n".format(cmd=cmd)
        p = subprocess.Popen(cmd,shell=True)
        p.wait()
        if p.returncode: ## Non zero returncode
            print "Command failed : {cmd}\n".format(cmd=cmd)

if __name__ == '__main__':
    luigi_config = sys.argv[1]
    main(luigi_config)
