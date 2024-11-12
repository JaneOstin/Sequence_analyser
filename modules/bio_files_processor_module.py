import os


def path_output(input_fasta: str, output_fasta: str):
    path = '/'.join(input_fasta.split("/")[:-1])
    if output_fasta == '':
        return os.path.join(path, 'output_' + input_fasta.split("/")[-1])
    else:
        return os.path.join(path, output_fasta)   


def write_fasta_file(lines_fasta: list, output_file: str):
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    with open(output_file, mode='w') as output_fasta:
        for line_number in range(len(lines_fasta)-1):
            if lines_fasta[line_number].startswith('>'):
                output_fasta.write(lines_fasta[line_number])
            elif lines_fasta[line_number+1].startswith('>'):
                output_fasta.write(lines_fasta[line_number])
            else:
                output_fasta.write(lines_fasta[line_number].strip('\n'))
        output_fasta.write(lines_fasta[len(lines_fasta)-1].strip('\n'))


def find_protein(lines_file: list):
    list_protein = []
    for line_number in range(len(lines_file)):
        if lines_file[line_number].startswith('Query #'):
            position = lines_file[line_number+3].find('Scientific')
            list_protein.append(lines_file[line_number+5][:position])
    return list_protein

def write_output_protein_file(list_protein: list, output_protein_file: str):
    os.makedirs(os.path.dirname(output_protein_file), exist_ok=True)
    list_protein = sorted(list_protein, key=lambda s: s.casefold())
    with open(output_protein_file, mode='w') as output_file:
        for protein in list_protein:
            output_file.write(protein + '\n')


def find_genes_proteins(input_gbk: str):
    with open(input_gbk) as input_file:
        lines_file = input_file.readlines()
    proteins_list = []
    genes_list = []
    for line_number in range(len(lines_file)):
        if lines_file[line_number].find('/gene=') == -1:
            continue
        gene_name = lines_file[line_number].split('"')[1]
        seq = ''
        line_number_protein_start = line_number + 1
        while lines_file[line_number_protein_start].find('/translation') == -1:
            line_number_protein_start += 1
        line_number_protein_end = line_number_protein_start + 1
        while lines_file[line_number_protein_end].find('"') == -1:
            line_number_protein_end += 1
        seq += lines_file[line_number_protein_start].split('"')[1].strip()
        for number in range(line_number_protein_start + 1, line_number_protein_end):
            seq += lines_file[number].strip()
        seq += lines_file[line_number_protein_end].split('"')[0].strip()
        genes_list.append(gene_name)
        proteins_list.append(seq)
    return genes_list, proteins_list


def find_genes_before_after(genes: list, genes_list: list, proteins_list: list,  n_before: int, n_after: int):
    genes_output = []
    proteins_output = []
    for gene in genes:
        if gene in genes_list:
            index = genes_list.index(gene)
            number_before = index-1
            count_before = 1
            while (number_before > -1 and count_before < n_before + 1):
                genes_output.append(genes_list[number_before])
                proteins_output.append(proteins_list[number_before])
                number_before += 1
                count_before += 1
            number_after = index+1
            count_after = 1
            while (number_after < len(genes_list) and count_after < n_after + 1) == True:
                genes_output.append(genes_list[number_after])
                proteins_output.append(proteins_list[number_after])
                number_after += 1
                count_after += 1
    temp_list = []
    for gene in genes_output:
        if gene not in temp_list:
            temp_list.append(gene)
    genes_output = temp_list
    temp_list = []
    for protein in proteins_output:
        if protein not in temp_list:
            temp_list.append(protein)
    proteins_output = temp_list
    return genes_output, proteins_output


def path_output_fasta(input_gbk: str, output_fasta: str):
    path = '/'.join(input_gbk.split("/")[:-1])
    if output_fasta == '':
        file_name = input_gbk.split("/")[-1]
        return os.path.join(path, 'output_' + file_name.split(".gbk")[0] + '.fasta')
    return os.path.join(path, output_fasta)
    

def write_output_fasta_file(genes_output: list, proteins_output: list, output_fasta_file: str):
    os.makedirs(os.path.dirname(output_fasta_file), exist_ok=True)
    with open(output_fasta_file, mode='w') as output_file:
        for index in range(len(genes_output)):
            output_file.write('>' + genes_output[index] + '\n' + proteins_output[index] + '\n')
