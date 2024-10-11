import os


def quality(seq: str, qual_treshhold: float = 0):
    qual_sum = 0
    len_seq = len(seq)
    for qual in seq:
        qual_sum += ord(qual) - 33
    if len_seq == 0:
        return None
    else:
        return qual_sum/len_seq


def gc_content_filter(gc_content_seq: float, gc_bounds: tuple):
    if gc_content_seq is None :
        return False
    elif gc_content_seq >= gc_bounds[0] and gc_content_seq <= gc_bounds[1]:
        return True
    else:
        return False
    

def lenght_filter(seq_len: int, length_bounds: tuple):
    if seq_len >= length_bounds[0] and seq_len <= length_bounds[1]:
        return True
    else:
        return False
    

def qual_filter(qual_seq: float, quality_threshold: float):
    if qual_seq is None:
        return False
    elif qual_seq > quality_threshold:
        return True
    else:
        return False
    

def read_fastq_file(input_fastq: str):
    with open(input_fastq) as input_file:
        lines = input_file.readlines()
    return lines


def fastq_dict(lines_fastq: list):
    seqs = {}
    for number_line in range(0, len(lines_fastq), 4):
        seqs[lines_fastq[number_line]] = [lines_fastq[number_line+1][:-1],lines_fastq[number_line+2][:-1],lines_fastq[number_line+3][:-1]]
    return seqs


def path_output(input_fastq: str, output_fastq: str):
    path = '/'.join(input_fastq.split("/")[:-1]) + '/filtered'
    if output_fastq == '':
        return os.path.join(path, 'output_' + input_fastq.split("/")[-1])
    else:
        return os.path.join(path, output_fastq)
    

def write_output(filtred_fastq_seq: dict, output_file: str):
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    with open(output_file, mode='w') as output_fastq:
        for seq_name, (sequence, seq_name_double, seq_quality) in filtred_fastq_seq.items():
            output_fastq.write(seq_name + sequence + '\n' + seq_name_double + '\n' + seq_quality + '\n')
