"""Tests for mutate_genome.py. Run with: python scripts/test_mutate_genome.py"""
import random
from mutate_genome import make_mutations, apply_mutations


def test_toy_example():
    toy = "AGTGACAGTTTCAG"
    mutations = [
        {"pos": 1, "type": "SNP", "ref": "G", "alt": "C"},
        {"pos": 4, "type": "INS", "ref": "A", "alt": "ATT"},
        {"pos": 8, "type": "DEL", "ref": "TTT", "alt": "T"},
    ]
    result = apply_mutations(toy, mutations)
    assert result == "ACTGATTCAGTCAG", result


def test_mutations_match_reference():
    rng = random.Random(42)
    seq = ""
    for i in range(100000):
        seq = seq + rng.choice("ACGT")
    mutations = make_mutations(seq, 300, 10, 10, 50, seed=7)

    assert len(mutations) == 320
    for m in mutations:
        assert seq[m["pos"]:m["pos"] + len(m["ref"])] == m["ref"]
        assert m["ref"][0] == m["alt"][0] or m["type"] == "SNP"
        if m["type"] == "SNP":
            assert m["ref"] != m["alt"]
    for i in range(1, len(mutations)):
        assert mutations[i]["pos"] - mutations[i - 1]["pos"] >= 50


if __name__ == "__main__":
    test_toy_example()
    print("test_toy_example: PASS")
    test_mutations_match_reference()
    print("test_mutations_match_reference: PASS")
