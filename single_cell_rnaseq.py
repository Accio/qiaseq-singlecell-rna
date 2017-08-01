import os
import luigi

## Modules from this project
from extract_multiplex_region import extract_region
from demultiplex_cells import create_cell_fastqs
from align_transcriptome import star_alignment,annotate_bam_umi
from count_mt import count_mts


class config(luigi.config):
    '''
    '''
    star = luigi.Parameter()
    star_params = luigi.Parameter()
    genome_dir = luigi.Parameter()

class CheckDemultiplex(luigi.Target):
    '''
    '''

    def __init__(self,file_loc):
        '''
        '''
        self.file_loc = file_loc

    def exists(self):
        '''
        '''
        n=subprocess.check_output('wc -l %s'%file_loc).strip('\n')
        if n == 4:
            return True
        else:
            return False

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
    primer_file = luigi.Parameter()
    vector_sequence = luigi.Parameter()
    isolator = luigi.Parameter()
    cell_index_len = luigi.IntParameter()
    mt_len = luigi.IntParameter()
    num_cores = luigi.IntParameter()
    num_errors = luigi.IntParameter()

    def __init__(self,*args,**kwargs):
        ''' The constructor
        '''
        #super(ExtractMultiplexRegion,self).__init__(*args,**kwargs)
        ## Set up folder structure
        self.sample_dir = os.path.join(self.output_dir,self.sample_name)
        if not os.path.exists(self.sample_dir):
            os.makedirs(self.sample_dir)
        ## Create a directory for storing verification files for task completion
        self.target_dir = os.path.join(self.sample_dir,'targets')
        if not os.path.exists(self.target_dir):
            os.makedirs(self.target_dir)
        ## The verification file for this task
        self.verification_file = os.path.join(self.target_dir,
                                              self.__class__.__name__+
                                              '.verification.txt')
        self.multiplex_file = os.path.join(self.sample_dir,
                                           '%s_multiplex_region.txt'%sample_name)
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
                       self.mt_len,self.isolator,self.read2_fastq,
                       self.multiplex_file,self.num_cores)
        ## Create the verification file
        with open(self.verification_file,'w') as OUT:
            print >> OUT,"verification"

    def output(self):
        ''' Check for the existence of the output file
        '''
        return MyExtTask(self.verification_file)

class DeMultiplexer(luigi.Task):
    ''' Task for demultiplexing a fastq into individual cells
    '''

    ## The parameters for this task
    R1_fastq = luigi.Parameter()
    R2_fastq = luigi.Parameter()
    output_dir = luigi.Parameter()
    sample_name = luigi.Parameter()
    cell_index_file = luigi.Parameter()
    primer_file = luigi.Parameter()
    vector_sequence = luigi.Parameter()
    isolator = luigi.Parameter()
    cell_index_len = luigi.IntParameter()
    mt_len = luigi.IntParameter()
    num_cores = luigi.IntParameter()
    num_errors = luigi.IntParameter()

    def __init__(self,*args,**kwargs):
        ''' The constructor
        '''
        super(DeMultiplexer,self).__init__(*args,**kwargs)
        self.working_dir = os.path.join(self.output_dir,self.sample_name)
        self.metric_file = os.path.join(self.working_dir,
                                       '%s_demultiplex_stats.txt'%sample_name)
        ## The verification file for this task
        self.verification_file = os.path.join(self.target_dir,
                                              self.__class__.__name__+
                                              '.verification.txt')
        self.multiplex_file = os.path.join(self.output_dir,
                                           '%s_multiplex_region.txt'%sample_name)

    def requires(self):
        ''' We need the the ExtractMultiplexRegion task to be finished
        '''
        return self.clone(ExtractMultiplexRegion)

    def run(self):
        '''
        '''
        create_cell_fastqs(self.working_dir,self.metric_file,
                                             self.cell_index_file,
                                             self.multiplex_file,self.R1_fastq)
        ## Create the verification file
        with open(self.verification_file,'w') as OUT:
            print >> OUT,"verification"

    def output(self):
        ''' Hard to determine how many valid fastqs will be demultiplexed
        Regardless of anything 4 lines should be present in the metric file
        Checking for that.
        '''
        return MyExtTask(self.verification_file)

class Alignment(luigi.Task):
    '''
    Task for running STAR for alignment
    '''

    ## Define some parameters
    R1_fastq = luigi.Parameter()
    R2_fastq = luigi.Parameter()
    output_dir = luigi.Parameter()
    sample_name = luigi.Parameter()
    cell_index_file = luigi.Parameter()
    primer_file = luigi.Parameter()
    vector_sequence = luigi.Parameter()
    isolator = luigi.Parameter()
    cell_index_len = luigi.IntParameter()
    mt_len = luigi.IntParameter()
    num_cores = luigi.IntParameter()
    num_errors = luigi.IntParameter()

    cell_fastq = luigi.Parameter()
    cell_num = luigi.IntParameter()
    cell_index = luigi.Parameter()


    def __init__(self,*args,**kwargs):
        ''' Class constructor
        '''
        super(Alignment,self).__init__(*args,**kwargs)
        self.sample_dir = os.path.join(self.output_dir,self.sample_name)
        self.cell_dir = os.path.join(self.sample_dir,'Cell%i_%s'%(cell_num,
                                                                  cell_index))
        self.bam = os.path.join(self.cell_dir,'Aligned.sortedByCoord.out.bam')
        self.tagged_bam = os.path.join(self.cell_dir,'Aligned.sortedByCoord.out.tagged.bam')
        ## The verification file for this task
        self.verification_file = os.path.join(self.target_dir,
                                              self.__class__.__name__+
                                              '.'+str(cell_num)+
                                              '.verification.txt')
    def requires(self):
        ''' The dependency for this task is the completion of the
        Demultiplexing task
        '''
        return Demultiplexer(R1_fastq=self.R1_fastq,R2_fastq=self.R2_fastq,
                             output_dir=self.output_dir,
                             sample_name=self.sample_name,
                             cell_index_file=self.cell_index_file,
                             vector_sequence=self.vector_sequence,
                             isolator=self.isolator,
                             primer_file=self.primer_file,
                             cell_index_len=self.cell_index_len,
                             mt_len=self.mt_len,num_cores=self.num_cores,
                             num_errors=self.num_errors)

    def run(self):
        ''' The commands to run
        '''
        ## Do the alignment
        star_alignment(config().star,config().genome_dir,
                                           self.cell_dir,config().star_params,
                                           self.cell_fastq)
        ## Check if the bam file has any records ?
        ## Add bam tags
        annotate_bam_umi(self.multiplex_file,self.bam,self.tagged_bam,
                         programs().samtools)
        ## Create the verification file
        with open(self.verification_file,'w') as OUT:
            print >> OUT,"verification"

    def output(self):
        ''' The result from this task is the creation of the bam file
        '''
        return MyExtTask(self.verification_file)

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
    primer_file = luigi.Parameter()
    vector_sequence = luigi.Parameter()
    isolator = luigi.Parameter()
    cell_index_len = luigi.IntParameter()
    mt_len = luigi.IntParameter()
    num_cores = luigi.IntParameter()
    num_errors = luigi.IntParameter()

    cell_fastq = luigi.Parameter()
    cell_num = luigi.IntParameter()
    cell_index = luigi.Parameter()

    def __init__(self,*args,**kwargs):
        ''' The Class constructor
        '''
        super(CountMT,self).__init__(*args,**kwargs)
        self.sample_dir = os.path.join(self.output_dir,self.sample_name)
        self.cell_dir = os.path.join(self.sample_dir,'Cell%i_%s'%(cell_num,
        self.bam = os.path.join(self.cell_dir,'Aligned.sortedByCoord.out.bam')
        self.tagged_bam = os.path.join(self.cell_dir,'Aligned.sortedByCoord.out.tagged.bam')
                                                                  cell_index))
        self.outfile = os.path.join(self.cell_dir,'mt_count.txt')
        ## The verification file for this task
        self.verification_file = os.path.join(self.target_dir,
                                              self.__class__.__name__+
                                              '.'+str(cell_num)+
                                              '.verification.txt')

    def requires(self):
        ''' The requirement is the completion of the Alignment task
        '''
        return self.clone(Alignment)

    def run(self):
        '''
        '''

        count_mts(self.primer_file,self.tagged_bam,self.outfile)
        with open(self.verification_file,'w') as OUT:
            print >> OUT,"verification"

    def output(self):
        ''' The output from this task
        '''
        return MyExtTask(self.verification_file)

class JoinCountFiles(luigi.Task):
    ''' Task for joining MT count files
    '''

    # Parameters
    R1_fastq = luigi.Parameter()
    R2_fastq = luigi.Parameter()
    output_dir = luigi.Parameter()
    sample_name = luigi.Parameter()
    cell_index_file = luigi.Parameter()
    primer_file = luigi.Parameter()
    vector_sequence = luigi.Parameter()
    isolator = luigi.Parameter()
    cell_index_len = luigi.IntParameter()
    mt_len = luigi.IntParameter()
    num_cores = luigi.IntParameter()
    num_errors = luigi.IntParameter()

    def __init__(self,*args,**kwargs):
        ''' The class constructor
        '''
        super(PipelineWrapper,self).__init__(*args,**kwargs)
        self.sample_dir = os.path.join(self.output_dir,self.sample_name)
        ## The verification file for this task
        self.verification_file = os.path.join(self.target_dir,
                                              self.__class__.__name__+
                                              '.'+str(cell_num)+
                                              '.verification.txt')
        self.cell_indices = []
        with open(self.cell_index_file,'r') as IN:
            for cell_index in IN:
                self.cell_indices.append(cell_index.strip('\n'))

    def run(self):
        '''
        '''

        ## Schedule the dependencies first
        dependencies = []
        for i,cell_index in enumerate(self.cell_indices):
            cell_num = i+1
            cell_fastq = os.path.join(self.sample_dir,
                                      cell_number+'_'+index+'/cell_'+
                                      str(cell_number)+'_R1.fastq')
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
                                        cell_fastq=cell_fastq,
                                        cell_num=cell_num,
                                        cell_index=cell_index))
        yield dependencies

        ## Join the files
        merge_mt_files(self.sample_dir,self.count_file,len(self.cell_indices))
        with open(self.verification_file,'w') as OUT:
            print >> OUT,"verification"

    def output(self):
        ''' The output from this task
        '''

        return MyExtTask(self.verification_file)

