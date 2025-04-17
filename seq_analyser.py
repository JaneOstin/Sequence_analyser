from abc import ABC, abstractmethod
from Bio import SeqIO, SeqUtils

import os
import argparse
import logging


class BiologicalSequence(ABC):
    '''Abstract method for working with biological sequences, 
    containing 5 abstract methods'''
    @abstractmethod
    def __len__(self):
        pass

    @abstractmethod
    def __getitem__(self, index):
        pass

    @abstractmethod
    def __str__(self):
        pass

    @abstractmethod
    def __repr__(self):
        pass

    @abstractmethod
    def check_alphabet(self) -> bool:
        pass


class NucleicAcidSequence(BiologicalSequence):
    '''Class for dealing with nucleic acid sequences - parent class BiologicalSequence. 
    Contains dunder methods, as well as methods for alphabet validation, sequence flipping, 
    complement and both functions, and GC content counting'''
    def __init__(self, sequence: str):
        self.sequence = sequence

    def __len__(self):
        return len(self.sequence)

    def __getitem__(self, index):
        if isinstance(index, slice):
            return self.sequence[index.start:index.stop:index.step]
        else:
            return self.sequence[index]

    def __str__(self):
        return f'{self.sequence}'

    def __repr__(self):
        return f'Nucleic Acid Sequence: {self.sequence}'

    def check_alphabet(self) -> bool:
        if set(self.sequence.upper()).issubset(self.alphabet):
            return True
        else:
            return False

    def reverse(self):
        return self.sequence[::-1]

    def complement(self):
        try:
            return ''.join([self.vocab[i] for i in self.sequence])
        except:
            raise NotImplementedError('Make an example or check sequence for alphabet') from None

    def reverse_complement(self):
        return self.complement()[::-1]

    def gc_content(self):
        gc_seq = self.sequence.upper()
        len_seq = len(gc_seq)
        if len_seq == 0:
            return None
        else:
            count_gc = (gc_seq.count('G') + gc_seq.count('C')) * 100 / len_seq
            return count_gc


class DNASequence(NucleicAcidSequence):
    '''Class for working with DNA sequences - parent class NucleicAcidSequence. 
    Contains function for transcribing sequence'''
    alphabet = set(["A", "T", "G", "C"])
    vocab = {
        'A': 'T',
        'a': 't',
        'T': 'A',
        't': 'a',
        'G': 'C',
        'g': 'c',
        'C': 'G',
        'c': 'g'
    }

    def __init__(self, sequence: str):
        super().__init__(sequence)

    def transcribe(self):
        trans_seq = self.sequence.replace('T', 'U')
        trans_seq = trans_seq.replace('t', 'u')
        return trans_seq


class RNASequence(NucleicAcidSequence):
    '''Class for working with RNA sequences - parent class NucleicAcidSequence.'''
    alphabet = set(["A", "U", "G", "C"])
    vocab = {
        'A': 'U',
        'a': 'u',
        'U': 'A',
        'u': 'a',
        'G': 'C',
        'g': 'c',
        'C': 'G',
        'c': 'g'
    }

    def __init__(self, sequence: str):
        super().__init__(sequence)


class AminoAcidSequence(BiologicalSequence):
    '''Class for working with protein sequences - parent class BiologicalSequence.
    Contains dunder methods, as well as methods for alphabet validation and counting
    amino acid frequency in a sequence'''
    def __init__(self, sequence: str):
        self.sequence = sequence

    def __len__(self):
        return len(self.sequence)

    def __getitem__(self, index):
        if isinstance(index, slice):
            return self.sequence[index.start:index.stop:index.step]
        else:
            return self.sequence[index]

    def __str__(self):
        return f'{self.sequence}'

    def __repr__(self):
        return f'Amino Acid Sequence: {self.sequence}'

    def check_alphabet(self) -> bool:
        alphabet = set("ACDEFGHIKLMNPQRSTVWY")
        return set(self.sequence.upper()).issubset(alphabet)

    def aa_frequency(self):
        aa_freq = {}
        for amino_acid in self.sequence:
            if amino_acid in aa_freq:
                aa_freq[amino_acid] += 1
            else:
                aa_freq[amino_acid] = 1
        return aa_freq


def filter_fastq(
        input_fastq: str, output_fastq: str = '', gc_bounds: tuple = (0, 100),
        length_bounds: tuple = (0, 2**32),
        quality_threshold: float = 0):
    '''Function can filter fastq files by gc_content, length and quality
    gc_bounds, length_bounds and quality_threshold can also be int or float
    As a result it makes a new fastq file with filtered sequences
    '''
    fastq_dict = {record.id: record for record in SeqIO.parse(input_fastq, "fastq")}
    filtred_fastq_seq = []
    path = '/'.join(input_fastq.split("/")[:-1]) + '/filtered'

    logging.basicConfig(filename=os.path.join(path, 'fastq-filter.log'), filemode='w', level=logging.INFO)

    if not isinstance(gc_bounds, tuple):
        gc_bounds = (0, gc_bounds)

    if not isinstance(length_bounds, tuple):
        length_bounds = (0, length_bounds)

    if output_fastq == '':
        path = os.path.join(path, 'output_' + input_fastq.split("/")[-1])
        logging.error('No output file name, used the same as input')
    else:
        path = os.path.join(path, output_fastq)

    for seq_name, record in fastq_dict.items():
        seq_len = len(record.seq)
        qual_seq = sum(record.letter_annotations["phred_quality"]) / seq_len
        gc_content_seq = SeqUtils.gc_fraction(record.seq) * 100
        if (length_bounds[1] >= seq_len >= length_bounds[0]) and (gc_bounds[0] <= gc_content_seq <= gc_bounds[1]) and (qual_seq > quality_threshold):
            filtred_fastq_seq.append(SeqIO.SeqRecord(record.seq, id=seq_name, description=record.description, letter_annotations=record.letter_annotations))

    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as output_fastq:
        SeqIO.write(filtred_fastq_seq, output_fastq, "fastq")

    logging.info('Filtering is done!')
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter a FASTQ file by GC content, length, and quality threshold")

    parser.add_argument("input_fastq", help="Path for input file")
    parser.add_argument("-o", "--output_fastq", default='', help="Path for output file (same as input by default)")
    parser.add_argument("--gc", type=float, nargs='+', default=(0, 100), help="GC% to filter (--gc 40 60)")
    parser.add_argument("--length", type=int, nargs='+', default=(0, 2**32), help="Length of sequence (--length 50 150)")
    parser.add_argument("--quality", type=float, default=0, help="Treshold for quality  (например, --quality 30)")

    args = parser.parse_args()

    if len(args.gc) == 1:
        gc_bounds = (0, int(args.gc[0]))
    elif len(args.gc) == 2:
        gc_bounds = tuple(args.gc)
    else:
        parser.error("--gc takes 1 or 2 values")

    if len(args.length) == 1:
        length_bounds = (0, int(args.length[0]))
    elif len(args.length) == 2:
        length_bounds = tuple(args.length)
    else:
        parser.error("--length takes 1 or 2 values")

    filter_fastq(
        args.input_fastq,
        args.output_fastq,
        gc_bounds,
        length_bounds,
        args.quality
    )
