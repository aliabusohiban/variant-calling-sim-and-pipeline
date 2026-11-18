"""Add random SNPs and small indels to a genome.

Writes the mutated genome as FASTA and the true mutations as a VCF
(positions are on the original reference, 1-based).

Example:
    python scripts/mutate_genome.py --ref data/EcoliK12-MG1655.fasta \
        --out_prefix sim/ecoli --seed 1
"""
import random
from argparse import ArgumentParser

BASES = "ACGT"


def read_fasta(path):
    name = None
    pieces = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                name = line[1:].split()[0]
            else:
                pieces.append(line)
    return name, "".join(pieces).upper()


def choose_positions(genome_length, n, min_gap, rng):
    chosen = []
    while len(chosen) < n:
        pos = rng.randint(min_gap, genome_length - min_gap)
        too_close = False
        for p in chosen:
            if abs(p - pos) < min_gap:
                too_close = True
                break
        if not too_close:
            chosen.append(pos)
    return chosen


def make_mutations(seq, n_snps, n_ins, n_del, min_gap, seed):
    rng = random.Random(seed)
    total = n_snps + n_ins + n_del
    positions = choose_positions(len(seq), total, min_gap, rng)
    positions.sort()

    types = ["SNP"] * n_snps + ["INS"] * n_ins + ["DEL"] * n_del
    rng.shuffle(types)

    mutations = []
    for i in range(total):
        pos = positions[i]
        mtype = types[i]
        ref_base = seq[pos]

        if mtype == "SNP":
            options = []
            for b in BASES:
                if b != ref_base:
                    options.append(b)
            ref = ref_base
            alt = rng.choice(options)
        elif mtype == "INS":
            length = rng.randint(1, 10)
            inserted = ""
            for j in range(length):
                inserted = inserted + rng.choice(BASES)
            ref = ref_base
            alt = ref_base + inserted
        else:
            length = rng.randint(1, 10)
            ref = seq[pos:pos + 1 + length]
            alt = ref_base

        mutations.append({"pos": pos, "type": mtype, "ref": ref, "alt": alt})
    return mutations


def apply_mutations(seq, mutations):
    pieces = []
    current = 0
    for m in mutations:
        pieces.append(seq[current:m["pos"]])
        pieces.append(m["alt"])
        current = m["pos"] + len(m["ref"])
    pieces.append(seq[current:])
    return "".join(pieces)


def write_fasta(name, seq, path, width=70):
    with open(path, "w") as f:
        f.write(">" + name + "\n")
        for i in range(0, len(seq), width):
            f.write(seq[i:i + width] + "\n")


def write_vcf(mutations, chrom, chrom_length, path):
    with open(path, "w") as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write("##contig=<ID=%s,length=%d>\n" % (chrom, chrom_length))
        f.write('##INFO=<ID=TYPE,Number=1,Type=String,Description="SNP, INS or DEL">\n')
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for m in mutations:
            f.write("%s\t%d\t.\t%s\t%s\t.\tPASS\tTYPE=%s\n"
                    % (chrom, m["pos"] + 1, m["ref"], m["alt"], m["type"]))


def main():
    parser = ArgumentParser(description="Add random SNPs and indels to a genome")
    parser.add_argument("--ref", required=True, help="reference genome FASTA")
    parser.add_argument("--out_prefix", required=True,
                        help="writes <prefix>_mutated.fasta and <prefix>_truth.vcf")
    parser.add_argument("--snps", type=int, default=300)
    parser.add_argument("--ins", type=int, default=10)
    parser.add_argument("--dels", type=int, default=10)
    parser.add_argument("--min_gap", type=int, default=50)
    parser.add_argument("--seed", type=int, default=1)
    args = parser.parse_args()

    name, seq = read_fasta(args.ref)
    mutations = make_mutations(seq, args.snps, args.ins, args.dels,
                               args.min_gap, args.seed)
    mutated = apply_mutations(seq, mutations)

    write_fasta(name + "_mutated", mutated, args.out_prefix + "_mutated.fasta")
    write_vcf(mutations, name, len(seq), args.out_prefix + "_truth.vcf")

    inserted = 0
    deleted = 0
    for m in mutations:
        if m["type"] == "INS":
            inserted = inserted + len(m["alt"]) - 1
        if m["type"] == "DEL":
            deleted = deleted + len(m["ref"]) - 1
    expected = len(seq) + inserted - deleted

    print("genome:", name)
    print("mutations:", len(mutations))
    print("original length:", len(seq))
    print("mutated length: ", len(mutated))
    print("expected length:", expected)
    if len(mutated) != expected:
        raise SystemExit("ERROR: mutated length is wrong")


if __name__ == "__main__":
    main()
