import os
import glob
import sys
import luigi
import sqlite3
import ConfigParser
## Modules from this project
sys.path.append(os.path.join(os.path.dirname(
    os.path.realpath(__file__)),'core'))
from extract_multiplex_region import extract_region
from demultiplex_cells import create_cell_fastqs
from align_transcriptome import star_alignment,star_load_index,star_remove_index,annotate_bam_umi
from count_mt import count_umis,count_umis_wts
from merge_mt_files import merge_count_files,merge_metric_files
from combine_sample_results import combine_count_files,combine_cell_metrics
from create_annotation_tables import create_gene_tree

GENE_TREE = None 
class config(luigi.Config):
    ''' Initialize values from configuration file
    '''
    star = luigi.Parameter(description="Path to the STAR executable")
    star_params = luigi.Parameter(description="Params for STAR")
    star_load_params = luigi.Parameter(description="Params for STAR to load a genome file")
    genome_dir = luigi.Parameter(description="The path to the star index dir")
    seqtype = luigi.Parameter(description="Whether this is a targetted or wts experiment")
    primer_file = luigi.Parameter(description="The primer file,if wts this is not applicable")
    annotation_gtf = luigi.Parameter(description="Gencode annotation file")

class MyExtTask(luigi.ExternalTask):
    ''' Checks whether the file specified exists on disk
    '''
    file_loc = luigi.Parameter()
    def output(self):
        return luigi.LocalTarget(self.file_loc)

class ExtractMultiplexRegion(luigi.Task):
    ''' Task for extracted the <cell_index><mt> region
    from R2 reads
    '''

    ## The parameters for this task
    R1_fastq = luigi.Parameter()
    R2_fastq = luigi.Parameter()
    output_dir = luigi.Parameter()
    sample_name = luigi.Parameter()
    cell_index_file = luigi.Parameter()
    vector_sequence = luigi.Parameter()
    isolator = luigi.Parameter()
    cell_index_len = luigi.IntParameter()
    mt_len = luigi.IntParameter()
    num_cores = luigi.IntParameter()
    num_errors = luigi.IntParameter()
    instrument = luigi.Parameter()

    def __init__(self,*args,**kwargs):
        ''' The constructor
        '''
        super(ExtractMultiplexRegion,self).__init__(*args,**kwargs)
        ## Set up folder structure
        self.sample_dir = os.path.join(self.output_dir,self.sample_name)
        if not os.path.exists(self.sample_dir):
            os.makedirs(self.sample_dir)
        ## Create a directory for storing verification files for task completion
        self.target_dir = os.path.join(self.sample_dir,'targets')
        if not os.path.exists(self.target_dir):
            os.makedirs(self.target_dir)
        ## Create a directory for storing log files
        self.logdir = os.path.join(self.sample_dir,'logs')
        if not os.path.exists(self.logdir):
            os.makedirs(self.logdir)
        ## The verification file for this task
        self.verification_file = os.path.join(self.target_dir,
                                              self.__class__.__name__+
                                              '.verification.txt')
        self.multiplex_file = os.path.join(self.sample_dir,
                                           '%s_multiplex_region.txt'%self.sample_name)
    def requires(self):
        ''' Dependencies for this task
        The R2 fastq reads must be present
        '''
        return MyExtTask(self.R2_fastq)

    def run(self):
        ''' Run the function to extract the multiplex region,
        i.e. the <cell_index><MT> sequence
        The resultant file is a tsv <read_id> <cell_index> <MT>
        '''
        extract_region(self.vector_sequence,
                       self.num_errors,self.cell_index_len,
                       self.mt_len,self.isolator,self.R2_fastq,
                       self.multiplex_file,self.num_cores,
                       self.instrument)
        ## Create the verification file
        with open(self.verification_file,'w') as OUT:
            print >> OUT,"verification"

    def output(self):
        ''' Check for the existence of the output file
        '''
        return luigi.LocalTarget(self.verification_file)

class DeMultiplexer(luigi.Task):
    ''' Task for demultiplexing a fastq into individual cells
    '''
    ## The parameters for this task
    R1_fastq = luigi.Parameter()
    R2_fastq = luigi.Parameter()
    output_dir = luigi.Parameter()
    sample_name = luigi.Parameter()
    cell_index_file = luigi.Parameter()
    vector_sequence = luigi.Parameter()
    isolator = luigi.Parameter()
    cell_index_len = luigi.IntParameter()
    mt_len = luigi.IntParameter()
    num_cores = luigi.IntParameter()
    num_errors = luigi.IntParameter()
    instrument = luigi.Parameter()

    def __init__(self,*args,**kwargs):
        ''' The constructor
        '''
        super(DeMultiplexer,self).__init__(*args,**kwargs)
        self.sample_dir = os.path.join(self.output_dir,self.sample_name)
        self.temp_metric_file = os.path.join(self.sample_dir,
                                       '%s_read_stats.temp.txt'%self.sample_name)
        self.target_dir = os.path.join(self.sample_dir,'targets')
        ## The verification file for this task
        self.verification_file = os.path.join(self.target_dir,
                                              self.__class__.__name__+
                                              '.verification.txt')
        self.multiplex_file = os.path.join(self.sample_dir,
                                           '%s_multiplex_region.txt'%self.sample_name)

    def requires(self):
        ''' We need the the ExtractMultiplexRegion task to be finished
        '''
        return self.clone(ExtractMultiplexRegion)

    def run(self):
        '''
        '''
        if config().seqtype.upper() == 'WTS':
            create_cell_fastqs(self.sample_dir,self.temp_metric_file,
                               self.cell_index_file,self.multiplex_file,
                               self.R1_fastq,True)
        else:
            create_cell_fastqs(self.sample_dir,self.temp_metric_file,
                               self.cell_index_file,
                               self.multiplex_file,self.R1_fastq)
        ## Create the verification file
        with open(self.verification_file,'w') as OUT:
            print >> OUT,"verification"

    def output(self):
        '''
        '''
        return luigi.LocalTarget(self.verification_file)

class LoadGenomeIndex(luigi.Task):
    ''' Task for loading genome index for STAR
    '''
    ## Define some parameters
    output_dir = luigi.Parameter(significant=False)

    def __init__(self,*args,**kwargs):
        ''' Class constructor
        '''
        super(LoadGenomeIndex,self).__init__(*args,**kwargs) 
        self.target_dir = os.path.join(self.output_dir,'targets')
        ## The verification file for this task
        self.verification_file = os.path.join(self.target_dir,
                                              self.__class__.__name__+
                                              '.verification.txt')

    def requires(self):
        ''' The dependency for this task is the existence of the genome dir
        '''
        return MyExtTask(config().genome_dir)

    def run(self):
        '''
        '''
        star_load_index(config().star,config().genome_dir,config().star_load_params)
        ## Create the verification file
        with open(self.verification_file,'w') as OUT:
            print >> OUT,"verification"

    def output(self):
        ''' The output from this task is the verification task
        '''
        return luigi.LocalTarget(self.verification_file)

class Alignment(luigi.Task):
    ''' Task for running STAR for alignment
    '''
    ## Define some parameters
    R1_fastq = luigi.Parameter()
    R2_fastq = luigi.Parameter()
    output_dir = luigi.Parameter()
    sample_name = luigi.Parameter()
    cell_index_file = luigi.Parameter()
    vector_sequence = luigi.Parameter()
    isolator = luigi.Parameter()
    cell_index_len = luigi.IntParameter()
    mt_len = luigi.IntParameter()
    num_cores = luigi.IntParameter()
    num_errors = luigi.IntParameter()
    instrument = luigi.Parameter()

    cell_fastq = luigi.Parameter()
    cell_num = luigi.IntParameter()
    cell_index = luigi.Parameter()


    def __init__(self,*args,**kwargs):
        ''' Class constructor
        '''
        super(Alignment,self).__init__(*args,**kwargs)
        self.sample_dir = os.path.join(self.output_dir,self.sample_name)
        self.multiplex_file = os.path.join(self.sample_dir,
                                           '%s_multiplex_region.txt'%self.sample_name)

        self.cell_dir = os.path.join(self.sample_dir,'Cell%i_%s'%(self.cell_num,
                                                                  self.cell_index))
        self.bam = os.path.join(self.cell_dir,'Aligned.sortedByCoord.out.bam')
        self.tagged_bam = os.path.join(self.cell_dir,'Aligned.sortedByCoord.out.tagged.bam')
        self.target_dir = os.path.join(self.sample_dir,'targets')
        ## The verification file for this task
        self.verification_file = os.path.join(self.target_dir,
                                              self.__class__.__name__+
                                              '.'+str(self.cell_num)+
                                              '.verification.txt')
    def requires(self):
        '''
        '''
        yield LoadGenomeIndex(output_dir=self.output_dir)
        yield self.clone(DeMultiplexer)

    def run(self):
        ''' The commands to run
        '''
        if os.path.getsize(self.cell_fastq) > 0: ## Make sure the file is not empty
            ## Do the alignment
            star_alignment(config().star,config().genome_dir,
                           os.path.join(self.cell_dir,''),config().star_params,
                           self.cell_fastq)
            ## Check if the bam file has any records ?
            ## Add bam tags
            annotate_bam_umi(self.multiplex_file,self.bam,self.tagged_bam)

        ## Create the verification file
        with open(self.verification_file,'w') as OUT:
            print >> OUT,"verification"

    def output(self):
        ''' The result from this task is the creation of the bam file
        '''
        return luigi.LocalTarget(self.verification_file)

class CountMT(luigi.Task):
    ''' Task for counting MTs, presumably this is the final step
    which gives us a Primer/Gene x Cell count matrix file
    Will likely have some wrapper task to enapsulate this to parallelize
    by Cells
    '''
    ## Parameters
    R1_fastq = luigi.Parameter()
    R2_fastq = luigi.Parameter()
    output_dir = luigi.Parameter()
    sample_name = luigi.Parameter()
    cell_index_file = luigi.Parameter()
    vector_sequence = luigi.Parameter()
    isolator = luigi.Parameter()
    cell_index_len = luigi.IntParameter()
    mt_len = luigi.IntParameter()
    num_cores = luigi.IntParameter()
    num_errors = luigi.IntParameter()
    instrument = luigi.Parameter()

    cell_fastq = luigi.Parameter()
    cell_num = luigi.IntParameter()
    cell_index = luigi.Parameter()

    def __init__(self,*args,**kwargs):
        ''' The Class constructor
        '''
        super(CountMT,self).__init__(*args,**kwargs)
        self.sample_dir = os.path.join(self.output_dir,self.sample_name)
        self.cell_dir = os.path.join(self.sample_dir,'Cell%i_%s'%(self.cell_num,self.cell_index))
        self.bam = os.path.join(self.cell_dir,'Aligned.sortedByCoord.out.bam')
        self.tagged_bam = os.path.join(self.cell_dir,'Aligned.sortedByCoord.out.tagged.bam')
        self.outfile = os.path.join(self.cell_dir,'umi_count.txt')
        self.outfile_primer = os.path.join(self.cell_dir,'umi_count.primers.txt')
        self.metricsfile = os.path.join(self.cell_dir,'read_stats.txt')
        self.target_dir = os.path.join(self.sample_dir,'targets')
        self.logdir = os.path.join(self.sample_dir,'logs')
        ## The verification file for this task
        self.verification_file = os.path.join(self.target_dir,
                                              self.__class__.__name__+
                                              '.'+str(self.cell_num)+
                                              '.verification.txt')
        self.logfile = os.path.join(self.logdir,
                                    self.__class__.__name__ +
                                    self.sample_name +
                                    '.'+str(self.cell_num)+'.log.txt')

    def requires(self):
        ''' The requirement is the completion of the Alignment task
        '''
        return self.clone(Alignment)

    def run(self):
        '''
        '''
        if os.path.getsize(self.cell_fastq) > 0:
            if config().seqtype.upper() == 'WTS':
                count_umis_wts(GENE_TREE,self.tagged_bam,self.outfile,
                               self.metricsfile,self.logfile)
            else:
                count_umis(config().primer_file,self.tagged_bam,self.outfile_primer,
                           self.outfile,self.metricsfile,self.num_cores)

        with open(self.verification_file,'w') as OUT:
            print >> OUT,"verification"

    def output(self):
        ''' The output from this task
        '''
        return luigi.LocalTarget(self.verification_file)

class JoinCountFiles(luigi.Task):
    ''' Task for joining MT count files
    '''
    # Parameters
    R1_fastq = luigi.Parameter()
    R2_fastq = luigi.Parameter()
    output_dir = luigi.Parameter()
    sample_name = luigi.Parameter()
    cell_index_file = luigi.Parameter()
    vector_sequence = luigi.Parameter()
    isolator = luigi.Parameter()
    mt_len = luigi.IntParameter()
    num_cores = luigi.IntParameter()
    num_errors = luigi.IntParameter()
    instrument = luigi.Parameter()

    def __init__(self,*args,**kwargs):
        ''' The class constructor
        '''
        super(JoinCountFiles,self).__init__(*args,**kwargs)
        self.sample_dir = os.path.join(self.output_dir,self.sample_name)
        self.count_file = os.path.join(self.sample_dir,
                                       self.sample_name+'.umi.counts.txt')
        self.temp_metric_file = os.path.join(self.sample_dir,
                                       '%s_read_stats.temp.txt'%self.sample_name)
        self.metric_file = os.path.join(self.sample_dir,
                                       '%s_read_stats.txt'%self.sample_name)
        self.metric_file_cell = os.path.join(self.sample_dir,
                                             '%s_cell_stats.txt'%self.sample_name)
        ## The verification file for this task
        self.target_dir = os.path.join(self.sample_dir,'targets')
        self.verification_file = os.path.join(self.target_dir,
                                              self.__class__.__name__+
                                              '.verification.txt')
        self.cell_indices = []
        i=0
        with open(self.cell_index_file,'r') as IN:
            for cell_index in IN:
                if i==0:
                    self.cell_index_len = len(cell_index.strip('\n'))
                    i+=1                                          
                self.cell_indices.append(cell_index.strip('\n'))

    
    def requires(self):
        '''
        '''
        ## Schedule the dependencies first
        dependencies = []
        for i,cell_index in enumerate(self.cell_indices):
            cell_num = i+1
            cell_dir = os.path.join(self.sample_dir,'Cell%i_%s'%(
                cell_num,cell_index))
            cell_fastq = os.path.join(cell_dir,'cell_'+str(cell_num)+
                                      '_R1.fastq')
            dependencies.append(CountMT(R1_fastq=self.R1_fastq,
                                        R2_fastq=self.R2_fastq,
                                        output_dir=self.output_dir,
                                        sample_name=self.sample_name,
                                        cell_index_file=self.cell_index_file,
                                        vector_sequence=self.vector_sequence,
                                        isolator=self.isolator,
                                        cell_index_len=self.cell_index_len,
                                        mt_len=self.mt_len,
                                        num_cores=self.num_cores,
                                        num_errors=self.num_errors,
                                        instrument=self.instrument,
                                        cell_fastq=cell_fastq,
                                        cell_num=cell_num,
                                        cell_index=cell_index))
        yield dependencies        

    def run(self):
        '''
        '''
        ## Merge gene level count files first
        files_to_merge = glob.glob(os.path.join(self.sample_dir,"*/umi_count.txt"))
        merge_count_files(self.sample_dir,self.count_file,self.sample_name,True,len(self.cell_indices),files_to_merge)
        ## Join the files
        if config().seqtype.upper() == 'WTS':
            wts = True
        else:
            ## Merge primer level count files
            wts = False
            files_to_merge = glob.glob(os.path.join(self.sample_dir,"*/umi_count.primers.txt"))
            merge_count_files(self.sample_dir,self.count_file,self.sample_name,wts,len(self.cell_indices),files_to_merge)            
        ## Merge metric files
        files_to_merge = glob.glob(os.path.join(self.sample_dir,"*/read_stats.txt"))
        merge_metric_files(self.sample_dir,self.temp_metric_file,self.metric_file,self.metric_file_cell,self.sample_name,wts,len(self.cell_indices),files_to_merge)

        with open(self.verification_file,'w') as OUT:
            print >> OUT,"verification"

    def output(self):
        ''' The output from this task
        '''
        return luigi.LocalTarget(self.verification_file)

class CombineSamples(luigi.Task):
    ''' Task for combining results from multiple samples
    '''
    # Parameters
    output_dir = luigi.Parameter()
    samples_cfg = luigi.Parameter()
    cell_index_file = luigi.Parameter()
    vector_sequence = luigi.Parameter()
    isolator = luigi.Parameter()
    mt_len = luigi.IntParameter()
    num_cores = luigi.IntParameter()
    num_errors = luigi.IntParameter()

    def __init__(self,*args,**kwargs):
        ''' The class constructor
        '''
        super(CombineSamples,self).__init__(*args,**kwargs)
        self.combined_count_file = os.path.join(self.output_dir,'combined.umi.counts.txt')
        self.combined_cell_metrics_file = os.path.join(self.output_dir,'combined.cell.metrics.txt')
        ## The verification file for this task
        self.target_dir = os.path.join(self.output_dir,'targets')
        if not os.path.exists(self.target_dir):
            os.makedirs(self.target_dir)
        self.verification_file = os.path.join(self.target_dir,
                                              self.__class__.__name__+
                                              '.verification.txt')
        ## Annotation information from gencode        
        global GENE_TREE
        GENE_TREE = create_gene_tree(config().annotation_gtf)
        
    def requires(self):
        '''
        '''
        dependencies = []
        parser = ConfigParser.ConfigParser()
        parser.read(self.samples_cfg)        
        for section in parser.sections():            
            sample_name = section
            R1_fastq = parser.get(section,'R1_fastq')
            R2_fastq = parser.get(section,'R2_fastq')
            instrument = parser.get(section,'Instrument')
            dependencies.append(
                JoinCountFiles(
                    R1_fastq=R1_fastq,R2_fastq=R2_fastq,
                    output_dir=self.output_dir,sample_name=sample_name,
                    cell_index_file=self.cell_index_file,vector_sequence=self.vector_sequence,
                    isolator=self.isolator,mt_len=self.mt_len,num_cores=self.num_cores,
                    num_errors=self.num_errors,instrument=instrument
                )
            )
            
        
        yield dependencies        
        
            
    def run(self):
        '''
        '''
        combine_count_files(self.output_dir,self.combined_count_file)
        combine_cell_metrics(self.output_dir,self.combined_cell_metrics_file)
        with open(self.verification_file,'w') as IN:
            IN.write('done\n')
        
    def output(self):
        '''
        '''
        return luigi.LocalTarget(self.verification_file)
       
        
        
