from modules.run_dna_rna_tools import nucl_acid
from modules.run_dna_rna_tools import proc_dna
from modules.run_dna_rna_tools import proc_rna
from modules.run_dna_rna_tools import gc_content
from modules.filter_fastq import quality


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
        seqs: dict, gc_bounds: tuple = (0, 100),
        length_bounds: tuple = (0, 2**32),
        quality_threshold: float = 0):
    '''Function can filter fastq_seq by gc_content, length and quality
    gc_bounds, length_bounds and quality_threshold can also be int or float
    '''
    filtred_fastq_seq = {}
    if not isinstance(gc_bounds, tuple):
        gc_bounds = (0, gc_bounds)
    if not isinstance(length_bounds, tuple):
        length_bounds = (0, length_bounds)
    for key, value in seqs.items():
        qual_seq = quality(value[1], quality_threshold)
        gc_content_seq = gc_content(value[0])
        seq_len = len(value[0])
        if gc_content_seq >= gc_bounds[0] and gc_content_seq <= gc_bounds[1]:
            if seq_len >= length_bounds[0] and seq_len <= length_bounds[1]:
                if qual_seq > quality_threshold:
                    filtred_fastq_seq.update({key: value})
    return filtred_fastq_seq
