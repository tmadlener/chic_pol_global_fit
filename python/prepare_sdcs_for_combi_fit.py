#!/usr/bin/env python3
"""Script (one-off) to prepare all the SDCs that are then used for the fit
where NLO and NLO+LP are combined
"""

from sdc import sdc

from functools import lru_cache

# Yay for completely unnecessary optimizations in these scripts
@lru_cache
def read_sdc(filename, column):
    return sdc.read_from_file(filename, column)


# The directly usable input sdcs where Carlos didn't have to do a smoothing of
# the NLO parts
INPUT_SDCS = (
    ('1s0octetunpol', '1S0_8_total'),
    ('3s1octetlong', '3S1_8_long'),
    ('3s1octetunpol', '3S1_8_total'),
    ('3pjoctetunpol', '3PJ_8_total'),
    ('3pjoctetlong', '3PJ_8_long'),
    ('3s1singlettransverse', '3S1_1_trans'),
    ('3s1singletunpol', '3S1_1_total')
)

# The 3Pi singlets for the chic1 and chic2 where Carlos smoothed the NLO parts
SMOOTHED_SDCS = (
    ('3p2singlet', '3P2_1', 3), # 3 helicity eigenstates: h=0, |h|=1, |h|=2
    ('3p1singlet', '3P1_1', 2)  # 2 helicity eigenstates: h=0, |h|=1
)


def write_sdc(sdc, outfile, header=''):
    """Write an sdc to the desired file"""
    with open(outfile, 'w') as output:
        if header:
            output.write(header)

        for i in range(sdc.size()):
            pt, val = sdc[i]
            output.write(f'{pt} {val}\n')


def split_input_sdc(infile, outfile, sdc_col):
    """Read the SDC from an input file with multiple SDC values per pT point and
    write a new file containing only two columns"""
    input_sdc = read_sdc(infile, sdc_col)
    header = f'''# SDC file containing only one column of values created from
# original inputfile: {infile}
# using column {sdc_col} for the SDC values and column 0 for the pt values
################################################
# pT SDC
################################################
'''
    write_sdc(input_sdc, outfile, header)


def split_nlo_directly_usable():
    """Split the directly usable NLO SDCs from (i.e. the ones where no intervention
    from Carlos was necessary)
    """
    for inbase, outbase in INPUT_SDCS:
        # NLO is stored in the 4th column (0-based indexing)
        split_input_sdc(f'../SDCs_Chung/{inbase}.txt', f'./sdcs/NLO/{outbase}.txt', 3)


def split_nlo_smoothed():
    """Split the files where the SDCs have been smoothed by Carlos (@NLO)"""
    for inbase, outbase, max_hel in SMOOTHED_SDCS:
        for h in range(max_hel):
            split_input_sdc(f'../SDCs_Chung/{inbase}_NLO.txt',
                            f'./sdcs/NLO/{outbase}_h{h}.txt',
                            h + 1)


def split_lp_directly_usable():
    """Calculate the LP part as NLO+LP - NLO and write that to a separate file"""
    def get_lp_part(inbase):
        """Get the LP part from the file"""
        nlo = read_sdc(f'../SDCs_Chung/{inbase}.txt', 3)
        # NLO+LP is stored in the 2nd column
        nlo_lp = read_sdc(f'../SDCs_Chung/{inbase}.txt', 1)

        return nlo_lp - nlo

    def get_header(inbase):
        """Get the header"""
        return f'''# SDC file containing only one column of LP correction values created from
# original file: ../SDCs_Chung/{inbase}.txt
# This file contains only the LP = NLO+LP - NLO part of the SDCs
# The original values are taken from column 3 (NLO) and column 1 (NLO+LP) from the original file
################################################
# pT SDC
################################################
'''


    for inbase, outbase in INPUT_SDCS:
        write_sdc(get_lp_part(inbase), f'./sdcs/LP/{outbase}.txt', get_header(inbase))


def split_lp_smoothed():
    """Calculate the LP part as NLO+LP - nLO and write that to a separate file for
    the parts where Carlos has applied the smoothing"""
    def get_lp_part(inbase, helicity):
        """Get the LP part from the file"""
        nlo = read_sdc(f'../SDCs_Chung/{inbase}_NLO.txt', helicity + 1)
        nlo_lp = read_sdc(f'../SDCs_Chung/{inbase}_NLOLP.txt', helicity + 1)
        return nlo_lp - nlo

    def get_header(inbase, helicity):
        return f'''# SDC file containing only one column of LP correction values created from
# original files: ../SDCs_Chung/{inbase}_NLO.txt and ../SDCs_Chung/{inbase}_NLOLP.txt
# Using column {helicity + 1} for helicity eigenstate h={helicity}
# LP = NLO+LP - NLO (where the values are taken from the corresponding individual files)
################################################
# pT SDC
################################################
'''


    for inbase, outbase, max_hel in SMOOTHED_SDCS:
        for h in range(max_hel):
            write_sdc(get_lp_part(inbase, h), f'./sdcs/LP/{outbase}_h{h}.txt', get_header(inbase, h))

def main():
    """Main"""
    split_nlo_directly_usable()
    split_nlo_smoothed()
    split_lp_directly_usable()
    split_lp_smoothed()

if __name__ == '__main__':
    main()
