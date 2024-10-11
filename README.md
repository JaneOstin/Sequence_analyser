# HW 5. Sequence analyser and bio files processor

## Content

* [What is it?](#What_is_it?)
* [Installation](#installation)
* [Running instructions](#running-instructions)
* [Examples](#examples)
* [See also](#see-also)

## What is it?

This is a toolkit for working with sequnces and bioinformation files. It can filter FASTQ files, convert multiline FASTA files to oneline, parse BLAST files, select genes from GBK file to FASTA file. 


## Installation

The source code is currently hosted on GitHub at: https://github.com/JaneOstin/Sequence_analyser

This repository contains a package with bioinformatic utilities. The main project file is `seq_analyser.py` and `bio_files_processor.py`. The auxiliary modules are located in the `modules` folder and are represented by the files `run_dna_rna_tools.py`, `filter_fastq.py` and `bio_files_processor_module.py`. A visual structure of the repository is presented below.

    ```
    -/
     |- README.md
     |- seq_analyser.py
     |- bio_files_processor.py
     |- modules/
           |- run_dna_rna_tools.py
           |- filter_fastq.py
           |- bio_files_processor_module.py
    ``` 

## Running instructions

### Sequence analyser

There are two functions in the `seq_analyser.py` file:

1. The `run_dna_rna_tools` function takes a sequence as input, checks if it is RNA/DNA, and processes it depending on the argument given. This function can do the following from the list:
    - `transcribe` - transcribe DNA
    - `reverse` - reverse the sequence
    - `complement` - build complementary sequence
    - `reverse_complement` - unfold and build complementary sequence
    - `gc_content_in_percentage` - count GC composition of the sequence 
    - The function outputs a single string or a list of strings after transformation. If the submitted string is not RNA/DNA or the argument is not in the list, then the function will report it

2. The `filter_fastq` function takes as input the path to the fastq file, as well as the optional arguments `output_fastq`, `gc_bounds`, `length_bounds`, `quality_threshold`:
    - `input_fastq` - the path to the fastq file
    - `output_fastq` - the result of function is a filtered file in the folder `/filtered/`, with output_ appended to the name, or saved exactly as specified
    - `gc_bounds` - GC composition interval (in per cent) for filtering (default is (0, 100), i.e. all reads are preserved). If a single number is passed in the argument, it is assumed to be the upper bound. Examples: `gc_bounds = (20, 80)` - save only reads with GC composition from 20 to 80%, `gc_bounds = 44.4` - save reads with GC composition less than 44.4%
    - `length_bounds` - length interval for filtering, all similar to `gc_bounds`, but defaults to (0, 2**32)
    - `quality_threshold` - threshold value of the average quality of the read for filtering, default value is 0 (phred33 scale). Reads with average quality for all nucleotides below the threshold are discarded

### Bio files processor

There are three functions in the `bio_files_processor.py` file:

1. The `convert_multiline_fasta_to_oneline` function takes as input the path to a FASTQ file in which the sequence is split into multiple lines, and returns a FASTQ file in which the sequence goes on a single line:
    - `input_fasta` - the path to the fasta file
    - `output_fasta` - the output file is named with output_ appended to the name, or saved exactly as specified

2. The `parse_blast_output` function reads the given txt file, for each query QUERY stores the name of the protein and outputs this in a new file in alphabetical order:
    - `input_file` - the path to the txt file as a result of BLAST
    - `output_file` - the output file is named with output_ appended to the name, or saved exactly as specified in format `.txt`

3. The `select_genes_from_gbk_to_fasta` function pulls genes and the sequence of their proteins next to the genes of interest from the GBK file:
    - `input_gbk` - the path to the GBK file
    - `genes` - genes of interest with respect to which function induces other genes
    - `n_before` - number of genes before (>=1). Default value is `1`
    - `n_after` - number of genes after (>=1). Default value is `1`
    - `output_fasta` - the output file is named with output_ appended to the name, or saved exactly as specified in format `.fasta`

## Examples

There are examples of using `run_dna_rna_tools` function:

```python
run_dna_rna_tools('ATG', 'transcribe') # 'AUG'
run_dna_rna_tools('ATG', 'reverse') # 'GTA'
run_dna_rna_tools('AtG', 'complement') # 'TaC'
run_dna_rna_tools('ATg', 'reverse_complement') # 'cAT'
run_dna_rna_tools('ATG', 'aT', 'reverse') # ['GTA', 'Ta']
```

After using `seq_analyser.py` and exactly `filter_fastq` function you will see a new folder. So, the structure of the repository will be:

    ```
    -/
     |- README.md
     |- seq_analyser.py
     |- bio_files_processor.py
     |- example_fastq.fastq
     |- modules/
           |- run_dna_rna_tools.py
           |- filter_fastq.py
           |- bio_files_processor_module.py
     |- filtered/
           |- output_example_fastq.fastq
    ``` 

Example of `filter_fastq` function:

    ```example_fastq.fastq with length_bounds = (20, 2**32)
    @SRX079804:1:SRR292678:1:1101:21885:21885 1:N:0:1 BH:ok
    ACAGCAACATAAACATGATGGGATGGCGTAAGCCCCCGAGATATCAGTTTACCCAGGATAAGAGATTAAATTATGAGCAACATTATTAA
    +SRX079804:1:SRR292678:1:1101:21885:21885 1:N:0:1 BH:ok
    FGGGFGGGFGGGFGDFGCEBB@CCDFDDFFFFBFFGFGEFDFFFF;D@DD>C@DDGGGDFGDGG?GFGFEGFGGEF@FDGGGFGFBGGD
    @SRX079804:1:SRR292678:1:1101:1182927:1182927 1:N:0:1 BH:changed:1
    TGAA
    +SRX079804:1:SRR292678:1:1101:1182927:1182927 1:N:0:1 BH:changed:1
    GGBH
    ``` 

    ```output_example_fastq.fastq
    @SRX079804:1:SRR292678:1:1101:21885:21885 1:N:0:1 BH:ok
    ACAGCAACATAAACATGATGGGATGGCGTAAGCCCCCGAGATATCAGTTTACCCAGGATAAGAGATTAAATTATGAGCAACATTATTAA
    +SRX079804:1:SRR292678:1:1101:21885:21885 1:N:0:1 BH:ok
    FGGGFGGGFGGGFGDFGCEBB@CCDFDDFFFFBFFGFGEFDFFFF;D@DD>C@DDGGGDFGDGG?GFGFEGFGGEF@FDGGGFGFBGGD
    ``` 

Example of `convert_multiline_fasta_to_oneline` function:

    ```input_fasta
    >5S_rRNA::NODE_272_length_223_cov_0.720238:18-129(+)
    ACGGCCATAGGACTTTGAAAGCACCGCATCCCGTCCGATCTGCGAAGTTAACCAAGATGCCGCCTGGTTAGTACCATGGTGGGGGACCACATGGGAATCCCT
    GGTGCTGTG
    ``` 

    ```output_fasta
    >5S_rRNA::NODE_272_length_223_cov_0.720238:18-129(+)
    ACGGCCATAGGACTTTGAAAGCACCGCATCCCGTCCGATCTGCGAAGTTAACCAAGATGCCGCCTGGTTAGTACCATGGTGGGGGACCACATGGGAATCCCTGGTGCTGTG
    ``` 

Example of output file `parse_blast_output` function (example of input file you can find in the repository `example_blast_results.txt`):

    ```output_file
    conjugal transfer protein TraA [Enterobacteriaceae]                            
    DinI-like family protein [Escherichia coli]                       
    DNA methylase [Enterobacteriaceae]                                               
    extended-spectrum class A beta-lactamase CTX-M-15 [Bacteria]      
    hypothetical protein [Enterobacteriaceae]                                    
    KlcA [Escherichia coli]                                               
    PilK [Escherichia coli]                                           
    ProQ/FINO family protein [Enterobacteriaceae]                    
    ``` 

Example of output file `select_genes_from_gbk_to_fasta` function (example of input file you can find in the repository `example_gbk.gbk`):

    ```output_file
    >yjiE_1
    MKNIETKWLYDFLTLEACRHFSQAAKERNLSQPAFSRRIKALEAAIGVVLFDRTTTPLQLTEEGKLFHSQTRSLLQQLECNLGELNGQSLLGVPNIKIAAAHSLSLSVLPKLVHSLTAYGGEFVYHVEAIDVVQAVNTLREGKSDFIISFRDEDLMQSPFCCLKLFESELYPVCAADAQGHPVFDITQPQVPLLNYTATSYMGRLVNRHLAEVGGITGRTIFISSMSELLKNMALNGYGIAWLPIWSIVDELQTKRLICLDAAKLTVPIQAYIYRMNTRLNRTAENLWRILQEHMPDDLIQQISMEEPARRN
    >stiP
    MQPHDTFTGSYQPGDVEFLLKPVVIEMTPVEQKEELIQSGKKHYSDMLSQEPAPTQWHLDLFHRALDRGAERLAKEVTQLAIALAKRFGDEPIVLASLVRAGVPLGVMLHQALRDMGKTSWHYGISIIRDRGIDGAALDVIEERHGTSGIVFVDGWTGKGAITGELVRALKDRPGYPEQPRLVVLADPCGCSWLAASDDDWLIPFGIMGAPVSGLISRSVWSSEGLHGCMVCEHLSEFECSRMLADTVAHFRKKLTPSSLAPLSWNTESARVLWQTSRDVIAFLADEFKVDSVNRIKPGIAEATRAVLRRVPDHVFVRSIDDPDVALLVGLAREKGIVVTEMGGTLGQYRAVTIIKKVL                
    ``` 

## See also

<img src="meme.jpg" width="800"/>

PS. It is not true picture, because while I was doing this toolkit I also didn't eat, drink and sleep as an experimental biologist 🥹

Don't hesitate to contact me on any questions about this toolkit!

Author: Antonina Kuznetsova
Email: kuznetsova.antonina@outlook.com
