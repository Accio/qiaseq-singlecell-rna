import luigi


class ExtractMultiplexRegion(luigi.Task):
    '''
    Task for extracted the <cell_index><mt> region
    from R2 reads
    '''

    ## Define some parameters
    cell_index_file = luigi.Parameter()
    isolator = luigi.Paramter()

    def __init__(self,*args,**kwargs):
        super(ExtractMultiplexRegion,self).__init__(*args,**kwargs)

    def requires(self):
        pass

    def run(self):
        pass

    def output(self):
        pass

class DeMultiplexer(luigi.Task):
    '''
    Task for demultiplexing a fastq into individual cells
    '''

    ## Define some parameters
    r1_fastq = luigi.Parameter()
    r2_fastq = luigi.Parameter()

    def __init__(self,*args,**kwargs):
        super(DeMultiplexer,self).__init__(*args,**kwargs)

    def requires(self):
        return ExtractMultiplexRegion(*args)

    def run(self):
        pass

    def output(self):
        pass

class TrimFastq(luigi.Task):
    '''
    Task for trimming a fastq file
    '''

    ## Define some parameters
    fastq = luigi.Parameter()
    seq = luigi.Parameter()
    loc = luigi.Parameter()
    cell_num = luigi.Parameter()

    def __init__(self,*args,**kwargs):
        super(TrimFastq,self).__init__(*args,**kwargs)

    def requires(self):
        return DeMultiplexer(*args)

    def output(self):
        pass

class STARAlignment(luigi.Task):
    '''
    Task for running STAR for alignment
    '''

    ## Define some parameters
    cell_num = luigi.Parameter()
    r1_fastq = luigi.Parameter()

    def __init__(self,*args,**kwargs):
        super(STARAlignment,self).__init__(*args,**kwargs)

    def requires(self):
        return TrimFastq(*args)

    def output(self):
        pass

class CountMT(luigi.Task):
    '''
    Task for counting MTs, presumably this is the final step
    which gives us a Primer/Gene x Cell count matrix file
    Will likely have some wrapper task to enapsulate this to parallelize
    by Cells
    '''

    ## Define some parameters

    def __init__(self,*args,**kwargs):
        super(CountMT,self).__init__(*args,**kwargs)

    def requires(self):
        return STARAlignment(*args)

    def output(self):
        pass

