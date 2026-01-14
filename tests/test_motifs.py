import pytest
from bioinf.motifs import motifs as m

def test_iupac_to_regex():
    assert m.iupac_to_regex("CCWGG") == "CC[AT]GG"

def test_iupac_to_regex_invalid():
    with pytest.raises(ValueError):
        m.iupac_to_regex("CCZGG")

def test_find_overlapping():
    assert m.find_overlapping("AAAA", "AA") == [0, 1, 2]

def test_prosite_to_regex_minimal():
    assert m.prosite_to_regex("C-x(2)-C") == "C.{2}C"

def test_prosite_empty_raises():
    with pytest.raises(ValueError):
        m.prosite_to_regex("")

def test_digest_dna_ecori():
    cuts, frags = m.digest_dna("GAATTCC", "G^AATTC")
    assert cuts == [1]
    assert frags == ["G", "AATTCC"]

def test_digest_dna_requires_caret():
    with pytest.raises(ValueError):
        m.digest_dna("GAATTCC", "GAATTC")

def test_build_pwm_simple_alphabet_AT():
    pwm = m.build_pwm(["AA", "AT"], alphabet="AT")
    assert pwm == [{"A": 1.0, "T": 0.0}, {"A": 0.5, "T": 0.5}]

def test_pssm_and_best_subsequence():
    pwm = m.build_pwm(["ACG", "ACG", "ATG"])
    pssm = m.pssm_from_pwm(pwm)
    best, pos, score = m.best_subsequence(pssm, "TTTACGAAA")
    assert best == "ACG"
    assert pos == 3
    assert score != float("-inf")

def test_best_subsequence_short_seq_raises():
    pwm = m.build_pwm(["ACG", "ACG", "ATG"])
    pssm = m.pssm_from_pwm(pwm)
    with pytest.raises(ValueError):
        m.best_subsequence(pssm, "AC")