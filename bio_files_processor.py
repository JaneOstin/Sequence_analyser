from modules.filter_fastq import read_fastq_file
from modules.bio_files_processor_module import path_output
from modules.bio_files_processor_module import write_fasta_file
from modules.bio_files_processor_module import find_protein
from modules.bio_files_processor_module import write_output_protein_file
from modules.bio_files_processor_module import find_genes_proteins
from modules.bio_files_processor_module import find_genes_before_after
from modules.bio_files_processor_module import path_output_fasta
from modules.bio_files_processor_module import write_output_fasta_file


def convert_multiline_fasta_to_oneline(input_fasta:str, output_fasta:str = ''):
    '''Function can convert multiple lines of sequence into one line'''
    lines_fasta = read_fastq_file(input_fasta)
    output_file = path_output(input_fasta, output_fasta)
    write_fasta_file(lines_fasta, output_file)


def parse_blast_output(input_file: str, output_file: str = ''):
    '''Function can extract genes for each query in alphabetic oder'''
    lines_file = read_fastq_file(input_file)
    output_protein_file = path_output(input_file, output_file)
    list_protein = find_protein(lines_file)
    write_output_protein_file(list_protein, output_protein_file)


def select_genes_from_gbk_to_fasta(input_gbk: str, genes: list, n_before: int = 1, n_after: int = 1, output_fasta: str = ''):
    '''Function can find genes and its protein sequence near genes of interest in n amount from them'''
    genes_list, proteins_list = find_genes_proteins(input_gbk)
    #matching_genes = [s for s in genes_list if any(xs in s for xs in genes)]
    genes_output, proteins_output = find_genes_before_after(genes, genes_list, proteins_list, n_before, n_after)
    output_fasta_file = path_output_fasta(input_gbk, output_fasta)
    write_output_fasta_file(genes_output, proteins_output, output_fasta_file)
