from modules.run_dna_rna_tools import nucl_acid
from modules.run_dna_rna_tools import proc_dna
from modules.run_dna_rna_tools import proc_rna
from modules.run_dna_rna_tools import gc_content
from modules.filter_fastq import quality
from modules.filter_fastq import gc_content_filter
from modules.filter_fastq import lenght_filter
from modules.filter_fastq import qual_filter
from modules.filter_fastq import read_fastq_file
from modules.filter_fastq import fastq_dict
from modules.filter_fastq import path_output
from modules.filter_fastq import write_output


def run_dna_rna_tools(*seq: str):
    '''Function can transcribe, reverse, complement,
    reverse_comp RNA/DNA and find gc_content
    '''
    result = []
    operation = seq[-1]
    for i in seq[:-1]:
        cur_nucl_acid = nucl_acid(i)
        if cur_nucl_acid == "dna":
            result.append(proc_dna(i, operation))
        elif cur_nucl_acid == "rna":
            result.append(proc_rna(i, operation))
        else:
            result.append("It is not DNA or RNA")
    if len(result) == 1:
        return result[0]
    else:
        return result


def filter_fastq(
        input_fastq:str, output_fastq:str = '', gc_bounds: tuple = (0, 100),
        length_bounds: tuple = (0, 2**32),
        quality_threshold: float = 0):
    '''Function can filter fastq files by gc_content, length and quality
    gc_bounds, length_bounds and quality_threshold can also be int or float
    As a result it makes a new fastq file with filtered sequences
    '''
    lines_fastq = read_fastq_file(input_fastq)
    seqs = fastq_dict(lines_fastq)
    filtred_fastq_seq = {}
    if not isinstance(gc_bounds, tuple):
        gc_bounds = (0, gc_bounds)
    if not isinstance(length_bounds, tuple):
        length_bounds = (0, length_bounds)
    for seq_name, (sequence, seq_name_double, seq_quality) in seqs.items():
        qual_seq = quality(seq_quality, quality_threshold)
        gc_content_seq = gc_content(sequence)
        seq_len = len(sequence)
        if gc_content_filter(gc_content_seq, gc_bounds) and lenght_filter(seq_len, length_bounds) and qual_filter(qual_seq, quality_threshold):
            filtred_fastq_seq[seq_name] = [sequence, seq_name_double, seq_quality]
    output_file = path_output(input_fastq, output_fastq)
    write_output(filtred_fastq_seq, output_file)
