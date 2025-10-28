# ---------------------------------- Imports --------------------------------- #
import math
from pathlib import Path
from test_utils import Timer

import pytest
from byu_pytest_utils import with_import, test_files, tier

# -------------------------------- Test tiers -------------------------------- #
baseline = tier('baseline', 1)
core = tier('core', 2)
stretch2 = tier('stretch2', 3)

# ----------------------------- Global variables ----------------------------- #
def read_sequence(file: Path) -> str:
    return ''.join(file.read_text().splitlines())

# ------------------------------- Baseline tests ------------------------------ #
@baseline
@with_import('alignment')
def test_small_alignment(align):
    score, aseq1, aseq2 = align('polynomial', 'exponential')
    assert score == -1
    assert aseq1 == 'polyn-omial'
    assert aseq2 == 'exponential'


@baseline
@with_import('alignment')
def test_tiny_dna_alignment(align):
    score, aseq1, aseq2 = align('ATGCATGC', 'ATGGTGC')
    assert score == -12
    assert aseq1 == 'ATGCATGC'
    assert aseq2 == 'ATG-GTGC'


@baseline
@with_import('alignment')
def test_tiny_dna_alignment_gaps(align):
    score, aseq1, aseq2 = align('ATATATATAT', 'TATATATATA')
    assert score == -17
    assert aseq1 == 'ATATATATAT-'
    assert aseq2 == '-TATATATATA'


@baseline
@with_import('alignment')
def test_small_dna_alignment_not_banded(align):
    score, aseq1, aseq2 = align('GGGGTTTTAAAACCCCTTTT', 'TTTTAAAACCCCTTTTGGGG')
    assert score == -8
    assert aseq1 == 'GGGGTTTTAAAACCCCTTTT----'
    assert aseq2 == '----TTTTAAAACCCCTTTTGGGG'


@baseline
@with_import('alignment')
def test_empty_sequences(align):
    score, aseq1, aseq2 = align('', '')
    assert score == 0
    assert aseq1 == ''
    assert aseq2 == ''


@baseline
@with_import('alignment')
def test_one_empty_sequence(align):
    score, aseq1, aseq2 = align('ACGT', '')
    assert score == 5 * 4  # Assuming gap penalty of 5
    assert aseq1 == 'ACGT'
    assert aseq2 == '----'


@baseline
@with_import('alignment')
def test_identical_sequences(align):
    score, aseq1, aseq2 = align('GATTACA', 'GATTACA')
    assert score == -21  # Assuming match score of -3
    assert aseq1 == 'GATTACA'
    assert aseq2 == 'GATTACA'


@baseline
@with_import('alignment')
def test_all_mismatch_sequences(align):
    score, aseq1, aseq2 = align('AAAA', 'TTTT')
    assert score == 1 * 4  # Assuming mismatch penalty of 1
    assert aseq1 == 'AAAA'
    assert aseq2 == 'TTTT'


@baseline
@with_import('alignment')
def test_alignment_with_long_gaps(align):
    score, aseq1, aseq2 = align('ATTTATTTA', 'AAA')
    assert aseq1 == 'ATTTATTTA'
    assert aseq2 == 'A---A---A'
    assert score == 21

@baseline
@with_import('alignment')
def test_different_match_score(align):
    score_default, a1_default, a2_default = align('ACGT', 'ACGT')
    score_bonus, a1_bonus, a2_bonus = align('ACGT', 'ACGT', match_award=-10)
    assert score_default == -12
    assert score_bonus == -40
    assert a1_bonus == a2_bonus == 'ACGT'


@baseline
@with_import('alignment')
def test_different_indel_penalty(align):
    score_low_penalty, _, _ = align('ACGT', 'AGT', indel_penalty=1)
    score_default_penalty, _, _ = align('ACGT', 'AGT')
    score_high_penalty, _, _ = align('ACGT', 'AGT', indel_penalty=10)
    assert score_low_penalty < score_default_penalty < score_high_penalty


@baseline
@with_import('alignment')
def test_different_substitution_penalty(align):
    score_default_sub, _, _ = align('ACGT', 'AGT')
    score_high_sub, _, _ = align('AAAA', 'TTTT', sub_penalty=10)
    assert score_default_sub < score_high_sub


@baseline
@with_import('alignment')
def test_tie_breaking(align):
    # generate a 3 way tie by modifying the match/sub awards
    score_3_way_tie, aseq1, aseq2 = align('ABDE', 'BADE', match_award=-2, sub_penalty=4, indel_penalty=5)
    assert aseq1 == 'ABDE'
    assert aseq2 == 'BADE'
    assert score_3_way_tie == 4

    score_diag_left_tie, bseq1, bseq2 = align('AA', 'A')
    assert score_diag_left_tie == 2
    assert bseq1 == 'AA'
    assert bseq2 == '-A'

    # generate another 2 way tie between left and up
    score_left_top_tie, cseq1, cseq2 = align('AT', 'AG', sub_penalty=15)
    assert score_left_top_tie == 7
    assert cseq1 == 'AT-'
    assert cseq2 == 'A-G'


@baseline
@with_import('alignment')
def test_alignment_scoring_effects(align):
    # alignment with default parameters
    score_default, aseq1_default, aseq2_default = align('AAAA', 'CCCC')

    # alignment with different scoring
    score_diff_scoring, aseq1_diff_scoring, aseq2_diff_scoring = align(
        'AAAA', 'CCCC', match_award=-1, indel_penalty=1, sub_penalty=5
    )

    assert aseq1_default == 'AAAA' and aseq2_default == 'CCCC'
    assert score_default == 4

    assert aseq1_diff_scoring == 'AAAA----' and aseq2_diff_scoring == '----CCCC'
    assert score_diff_scoring == 8


@baseline
@with_import('alignment')
def test_medium_dna_alignment(align):
    seq1 = 'ataagagtgattggcgatatcggctccgtacgtaccctttctactctcgggctcttccccgttagtttaaatctaatctctttataaacggcacttcc'
    seq2 = 'ataagagtgattggcgtccgtacgtaccctttctactctcaaactcttgttagtttaaatctaatctaaactttataaacggcacttcctgtgtgtccat'

    score, aseq1, aseq2 = align(seq1, seq2)

    expected_align1 = 'ataagagtgattggcgatatcggctccgtacgtaccctttctactctcgggctcttccccgttagtttaaatctaatct---ctttataaacggca----c----t-tcc--'
    expected_align2 = 'ataagagtgatt-g-g----c-g-tccgtacgtaccctttctactctcaaactctt----gttagtttaaatctaatctaaactttataaacggcacttcctgtgtgtccat'

    assert score == -116
    assert aseq1 == expected_align1
    assert aseq2 == expected_align2


@baseline
@with_import('alignment')
def test_large_dna_alignment(align):
    timer = Timer()
    for N in [10, 100, 1000, 1500, 2000, 3000]:
        timer.start_lap()
        seq1 = read_sequence(test_files / 'bovine_coronavirus.txt')[:N]
        seq2 = read_sequence(test_files / 'murine_hepatitus.txt')[:N]

        score, aseq1, aseq2 = align(seq1, seq2)

        with open(test_files / f'large_bovine_murine_align_{N}.txt') as file:
            expected_score, expected_align1, expected_align2 = file.read().splitlines()

        assert score == int(expected_score)
        assert aseq1 == expected_align1
        assert aseq2 == expected_align2

        if timer.lap_time() > 120:
            pytest.fail(f"The alignment of size {N} took too long")


# -------------------------------- Core tests -------------------------------- #
@core
@with_import('alignment')
def test_small_dna_alignment_banded(align):
    score, aseq1, aseq2 = align('GGGGTTTTAAAACCCCTTTT', 'TTTTAAAACCCCTTTTGGGG', banded_width=2)
    assert score == 6
    assert aseq1 == 'GGGGTTTTAAAACCCCTT--TT'
    assert aseq2 == '--TTTTAAAACCCCTTTTGGGG'


@core
@with_import('alignment')
def test_length_discrepancy_banded(align):
    score, aseq1, aseq2 = align('AAAA', 'AAAHHHAHHHHHHHHHHHHHHHHHH', banded_width=2)
    assert score == math.inf
    assert aseq1 == None
    assert aseq2 == None


@core
@with_import('alignment')
def test_alignment_changes_with_banded_width(align):
    # banded width 1
    score_narrow, a1_narrow, a2_narrow = align('tcgctcatatatccc', 'atataccctggggtg', banded_width=1, match_award=-1,
                                               sub_penalty=1, indel_penalty=1)
    assert score_narrow == 10
    assert a1_narrow == '-tcgctcatatatccc'
    assert a2_narrow == 'atataccct-ggggtg'

    # banded width 2
    score_mid, a1_mid, a2_mid = align('tcgctcatatatccc', 'atataccctggggtg', banded_width=2, match_award=-1,
                                      sub_penalty=1, indel_penalty=1)
    assert score_mid == 8
    assert a1_mid == '-t-cgctcat-atatccc'
    assert a2_mid == 'atatac-cctggggt--g'

    # banded width 7
    score_wide, a1_wide, a2_wide = align("tcgctcatatatccc", "atataccctggggtg", banded_width=7, match_award=-1,
                                         sub_penalty=1, indel_penalty=1)
    assert score_wide == 6
    assert a1_wide == 'tcgctcatatatccc-------'
    assert a2_wide == '------atata-ccctggggtg'


@core
@with_import('alignment')
def test_medium_dna_alignment_banded(align):
    seq1 = 'ataagagtgattggcgatatcggctccgtacgtaccctttctactctcgggctcttccccgttagtttaaatctaatctctttataaacggcacttcc'
    seq2 = 'ataagagtgattggcgtccgtacgtaccctttctactctcaaactcttgttagtttaaatctaatctaaactttataaacggcacttcctgtgtgtccat'

    score, aseq1, aseq2 = align(seq1, seq2, banded_width=2)

    expected_align1 = 'ataagagtgattggcg-atatcggctccgtacgtaccctttctactctcgggctcttccccgttagtttaaatctaatctctttataaacggcacttcc--'
    expected_align2 = 'ataagagtgattggcgtccgtacgtaccctttctactc-tcaaactcttgttagtttaaatctaatctaaactttataaacggcacttcctgtgtgtccat'

    assert score == -79
    assert aseq1 == expected_align1
    assert aseq2 == expected_align2


@core
@with_import('alignment')
def test_massive_dna_alignment_banded(align):
    timer = Timer()
    for N in [10, 100, 1000, 10000, 20000, 25000, 31000]:
        timer.start_lap()
        seq1 = read_sequence(test_files / 'bovine_coronavirus.txt')[:N]
        seq2 = read_sequence(test_files / 'murine_hepatitus.txt')[:N]

        score, aseq1, aseq2 = align(seq1, seq2, banded_width=3)

        with open(test_files / f'massive_bovine_murine_align_{N}.txt') as file:
            expected_score, expected_align1, expected_align2 = file.read().splitlines()

        assert score == int(expected_score)
        assert aseq1 == expected_align1
        assert aseq2 == expected_align2

        if timer.lap_time() > 10:
            pytest.fail(f"The unrestricted alignment of size {N} took too long")


# ------------------------------ Stretch 2 tests ----------------------------- #
@stretch2
@with_import('alignment')
def test_small_dna_alignment_open_gap(align):
    score, aseq1, aseq2 = align('GGGGTTTTAAAACCCCTTTT', 'TTTTAAAACCCCTTTTGGGG', gap_open_penalty=10, indel_penalty=.5)
    assert score == -24
    assert aseq1 == 'GGGGTTTTAAAACCCCTTTT----'
    assert aseq2 == '----TTTTAAAACCCCTTTTGGGG'


@stretch2
@with_import('alignment')
def test_small_dna_alignment_gap_open(align):
    score, seq1, seq2 = align('tcgctcatatatccc', 'atataccctggggtg', gap_open_penalty=10, indel_penalty=.5)
    assert score == 7
    assert seq1 == 'tcgctcatatatcc------c'
    assert seq2 == '------atataccctggggtg'


@stretch2
@with_import('alignment')
def test_medium_dna_alignment_gap_open(align):
    seq1 = 'ataagagtgattggcgatatcggctccgtacgtaccctttctactctcgggctcttccccgttagtttaaatctaatctctttataaacggcacttcc'
    seq2 = 'ataagagtgattggcgtccgtacgtaccctttctactctcaaactcttgttagtttaaatctaatctaaactttataaacggcacttcctgtgtgtccat'

    score, aseq1, aseq2 = align(seq1, seq2, gap_open_penalty=10, indel_penalty=.5)

    expected_align1 = 'ataagagtgattggcgatatcggctccgtacgtaccctttctactctcgggctcttccccgttagtttaaatctaatct---ctttataaacggcactt-----------cc'
    expected_align2 = 'ataagagtgattgg--------cgtccgtacgtaccctttctactctcaaactctt----gttagtttaaatctaatctaaactttataaacggcacttcctgtgtgtccat'

    assert score == -177
    assert aseq1 == expected_align1
    assert aseq2 == expected_align2


@stretch2
@with_import('alignment')
def test_alignment_changes_with_open_gap_penalties(align):
    seq1 = 'ataagagtgattggcgatatcggctccgtacgtaccctttctactctcgggctcttccccgttagtttaaatctaatctctttataaacggcacttcc'
    seq2 = 'ataagagtgattggcgtccgtacgtaccctttctactctcaaactcttgttagtttaaatctaatctaaactttataaacggcacttcctgtgtgtccat'

    score_a, seq1_a, seq2_a = align(seq1, seq2, gap_open_penalty=2, indel_penalty=.5)

    expected_align1 = 'ataagagtgattggcgatatcggctccgtacgtaccctttctactctcgggctcttccccgttagtttaaatctaatct---ctttataaacggcactt---------cc--'
    expected_align2 = 'ataagagtgattggcg---t-----ccgtacgtaccctttctactctcaaactctt----gttagtttaaatctaatctaaactttataaacggcacttcctgtgtgtccat'

    assert score_a == -221
    assert seq1_a == expected_align1
    assert seq2_a == expected_align2

    score_b, seq1_b, seq2_b = align(seq1, seq2, gap_open_penalty=5, indel_penalty=.5)

    expected_align1_b = 'ataagagtgattggcgatatcggctccgtacgtaccctttctactctcgggctcttccccgttagtttaaatctaatct---ctttataaacggcactt-----------cc'
    expected_align2_b = 'ataagagtgattgg--------cgtccgtacgtaccctttctactctcaaactctt----gttagtttaaatctaatctaaactttataaacggcacttcctgtgtgtccat'

    assert score_b == -197
    assert seq1_b == expected_align1_b
    assert seq2_b == expected_align2_b


@stretch2
@with_import('alignment')
def test_large_alignment_open_gap(align):
    seq1 = read_sequence(test_files / 'bovine_coronavirus.txt')[:1000]
    seq2 = read_sequence(test_files / 'murine_hepatitus.txt')[:1000]

    score, aseq1, aseq2 = align(seq1, seq2, gap_open_penalty=10, indel_penalty=.5)

    with open(test_files / f'large_bovine_murine_open_gap_1000.txt') as file:
        expected_score, expected_align1, expected_align2 = file.read().splitlines()

    assert score == int(expected_score)
    assert aseq1 == expected_align1
    assert aseq2 == expected_align2
