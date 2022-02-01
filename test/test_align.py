# Importing Dependencies
import pytest
from align import NeedlemanWunsch, read_fasta
import numpy as np

def test_nw_alignment():
    """
    TODO: Write your unit test for NW alignment
    using test_seq1.fa and test_seq2.fa by
    asserting that you have correctly filled out
    the your 3 alignment matrices.
    Use the BLOSUM62 matrix and a gap open penalty
    of -10 and a gap extension penalty of -1.
    """
    seq1, _ = read_fasta("./data/test_seq1.fa")
    seq2, _ = read_fasta("./data/test_seq2.fa")

    gap_open = -10
    gap_extend = -1

    f = NeedlemanWunsch("./substitution_matrices/BLOSUM62.mat",gap_open,gap_extend)
    alignment_score,seq1_align,seq2_align = f.align()

    #assert that the base cases for _align_matrix are correct
    assert f._align_matrix[0,0] == 0
    assert set(f._align_matrix[0,1:]) == {-np.inf} and set(f._align_matrix[1:,0]) == {-np.inf}

    #assert that the base cases for _gapA_matrix are correct
    assert set(f._gapA_matrix[:,0] == np.linspace(-10,-14,5)) == {True}
    assert set(f._gapA_matrix[0,1:]) == {-np.inf}

    #assert that the base cases for _gapB_matrix are correct
    assert set(f._gapB_matrix[0,:] == np.linspace(-10,-13,4)) == {True}
    assert set(f._gapB_matrix[1:,0]) == {-np.inf}

    #assert that a handful of representative elements are correct
    assert f._align_matrix[1,1] == 5
    assert f._gapA_matrix[1,1] == -22
    assert f._gapB_matrix[1,1] == -22
    assert f._align_matrix[1,2] == -11
    assert f._gapA_matrix[1,2] == -23
    assert f._gapB_matrix[1,2] == -6

def test_nw_backtrace():
    """
    TODO: Write your unit test for NW backtracing
    using test_seq3.fa and test_seq4.fa by
    asserting that the backtrace is correct.
    Use the BLOSUM62 matrix. Use a gap open
    penalty of -10 and a gap extension penalty of -1.
    """
    seq3, _ = read_fasta("./data/test_seq3.fa")
    seq4, _ = read_fasta("./data/test_seq4.fa")

    gap_open = -10
    gap_extend = -1

    f = NeedlemanWunsch("./substitution_matrices/BLOSUM62.mat",gap_open,gap_extend)
    alignment_score,seq3_align,seq4_align = f.align()

    #assert that backtrace returned the correct alignment and alignment score
    assert alignment_score == 17
    assert seq3_align == 'MAVHQLIRRP'
    assert seq4_align == 'M---QLIRHP'
