# Implement Prospector-like parsing of modifications. Returns the a.a.-only sequence and a list of the modifications.
__author__ = 'jakerosenberg'


def getmass(seq):

    std_aa_mass = {
        'G': 57.02146,
        'A': 71.03711,
        'S': 87.03203,
        'P': 97.05276,
        'V': 99.06841,
        'T': 101.04768,
        'C': 103.00919,
        'L': 113.08406,
        'I': 113.08406,
        'N': 114.04293,
        'D': 115.02694,
        'Q': 128.05858,
        'K': 128.09496,
        'E': 129.04259,
        'M': 131.04049,
        'H': 137.05891,
        'F': 147.06841,
        'R': 156.10111,
        'Y': 163.06333,
        'W': 186.07931,
    }

    seq_only_aa = ''.join([k for k in seq if k in std_aa_mass])
    onlymod = {n: 0 for n in range(1, len(seq_only_aa) + 1)}
    pos = 0

    for k in range(len(seq)):
        if seq[k] in std_aa_mass:
            pos += 1
        elif seq[k] == '(':  # catch the start of a modification of form (###)
            digit_position = k + 1
            modmass = []
            while seq[digit_position] != ')':
                modmass.append(seq[digit_position])
                digit_position += 1
            onlymod[pos] += float(''.join(modmass))

    return seq_only_aa, onlymod


