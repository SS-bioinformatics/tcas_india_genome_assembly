import sys

def stats(fasta_path):
    seqs = {}
    name = None
    chunks = []
    with open(fasta_path, "r") as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(chunks)
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line)
        if name is not None:
            seqs[name] = "".join(chunks)

    lengths = sorted([len(s) for s in seqs.values()], reverse=True)
    total = sum(lengths)
    n_count = sum(s.upper().count("N") for s in seqs.values())
    ungapped = total - n_count
    num_scaffolds = len(lengths)

    # N50 / L50
    half = total / 2.0
    cum = 0
    n50 = None
    l50 = None
    for i, l in enumerate(lengths, 1):
        cum += l
        if cum >= half:
            n50 = l
            l50 = i
            break

    print(f"File: {fasta_path}")
    print(f"  Num sequences (scaffolds): {num_scaffolds}")
    print(f"  Total size (Mbp): {total/1e6:.3f}")
    print(f"  Ungapped length (Mbp): {ungapped/1e6:.3f}")
    print(f"  Number of Ns (Mbp): {n_count/1e6:.3f}")
    print(f"  Scaffold N50 (Mbp): {n50/1e6:.3f}" if n50 else "  N50: n/a")
    print(f"  Scaffold L50: {l50}")
    print(f"  Per-sequence lengths:")
    for name_, l in sorted(seqs.items(), key=lambda kv: -len(kv[1])):
        print(f"    {name_}: {len(l)} bp ({l.upper().count('N')} Ns)")
    print()

if __name__ == "__main__":
    for p in sys.argv[1:]:
        stats(p)
