def quality(seq: str, qual_treshhold: float = 0):
    qual_sum = 0
    len_seq = len(seq)
    for qual in seq:
        qual_sum += ord(qual) - 33
    return qual_sum/len_seq
