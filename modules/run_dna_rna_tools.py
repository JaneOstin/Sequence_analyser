def nucl_acid(nuc_seq: str):
    dna = set(["A", "T", "G", "C"])
    rna = set(["A", "U", "G", "C"])
    if set(nuc_seq.upper()).issubset(dna):
        return "dna"
    elif set(nuc_seq.upper()).issubset(rna):
        return "rna"
    else:
        return "other"


def transcribe_proc(trans_seq: str):
    trans_seq = trans_seq.replace('T', 'U')
    trans_seq = trans_seq.replace('t', 'u')
    return trans_seq


def reverse_proc(rev_seq: str):
    return rev_seq[::-1]


def complement_proc(comp_seq: str, nucl_acid: str):
    res_comp_seq = ''
    dict_dna = {
        'A': 'T',
        'a': 't',
        'T': 'A',
        't': 'a',
        'G': 'C',
        'g': 'c',
        'C': 'G',
        'c': 'g'
    }
    dict_rna = {
        'A': 'U',
        'a': 'u',
        'U': 'A',
        'u': 'a',
        'G': 'C',
        'g': 'c',
        'C': 'G',
        'c': 'g'
    }
    if nucl_acid == 'dna':
        res_comp_seq = ''.join([dict_dna[i] for i in comp_seq])
    else:
        res_comp_seq = ''.join([dict_rna[i] for i in comp_seq])
    return res_comp_seq


def gc_content(gc_seq: str):
    gc_seq = gc_seq.upper()
    len_seq = len(gc_seq)
    if len_seq == 0:
        return None
    else:
        count_gc = (gc_seq.count('G') + gc_seq.count('C')) * 100 / len_seq
        return count_gc


def rev_comp(seq: str, nucl: str):
    return reverse_proc(complement_proc(seq, nucl))


def proc_dna(seq: str, oper: str):
    if oper == "transcribe":
        return transcribe_proc(seq)
    elif oper == "reverse":
        return reverse_proc(seq)
    elif oper == "complement":
        return complement_proc(seq, "dna")
    elif oper == "reverse_complement":
        return rev_comp(seq, "dna")
    elif oper == "gc_content_in_percentage":
        gc_perc = str(round(gc_content(seq))) + "% GC content"
        return gc_perc
    else:
        return "Unknown procedure"


def proc_rna(seq: str, oper: str):
    if oper == "transcribe":
        return "Unknown procedure in biology for RNA"
    elif oper == "reverse":
        return reverse_proc(seq)
    elif oper == "complement":
        return complement_proc(seq, "rna")
    elif oper == "reverse_complement":
        return rev_comp(seq, "rna")
    elif oper == "gc_content_in_percentage":
        gc_perc = str(round(gc_content(seq))) + "% GC content"
        return gc_perc
    else:
        return "Unknown procedure"
    