import argparse
import random


DNA = "ACGT"


def random_dna(length):
    """Generate a random DNA sequence."""
    return "".join(random.choices(DNA, k=length))


def make_haplotypes(length, heterozygosity):
    """Generate two haplotypes differing at heterozygous sites."""
    hap1 = random_dna(length)
    hap2 = list(hap1)

    for i, base in enumerate(hap1):
        if random.random() < heterozygosity:
            # Choose a different nucleotide
            hap2[i] = random.choice(
                [b for b in DNA if b != base]
            )

    return hap1, "".join(hap2)


def introduce_errors(read, error_rate):
    """Introduce independent substitution errors."""
    result = []

    for base in read:
        if random.random() < error_rate:
            base = random.choice(
                [b for b in DNA if b != base]
            )
        result.append(base)

    return "".join(result)


def simulate_reads(
    length,
    coverage,
    read_length,
    error_rate,
    heterozygosity,
):
    if length < read_length:
        raise ValueError(
            "Sequence length must be >= read length."
        )

    hap1, hap2 = make_haplotypes(
        length, heterozygosity
    )

    # Approximate number of reads needed for the
    # requested average coverage of the diploid genome.
    n_reads = round(coverage * length / read_length)

    for i in range(n_reads):
            # Randomly sample one of the two haplotypes
            genome = random.choice([hap1, hap2])

            # Choose a valid start position
            start = random.randint(
                0, length - read_length
            )

            read = genome[start:start + read_length]
            read = introduce_errors(read, error_rate)
            print(read)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Simulate simple diploid DNA reads."
    )

    parser.add_argument("--length", type=int, default=10000)
    parser.add_argument("--coverage", type=float, default=10)
    parser.add_argument("--read-length", type=int, default=100)
    parser.add_argument("--error-rate", type=float, default=0.01)
    parser.add_argument("--heterozygosity", type=float, default=0.001)
    parser.add_argument("--seed", type=int, default=None)

    args = parser.parse_args()

    if args.seed is not None:
        random.seed(args.seed)

    if args.length <= 0 or args.read_length <= 0:
        parser.error("Length and read length must be positive.")

    if args.coverage < 0:
        parser.error("Coverage must be non-negative.")

    if not 0 <= args.error_rate <= 1:
        parser.error("Error rate must be between 0 and 1.")

    if not 0 <= args.heterozygosity <= 1:
        parser.error("Heterozygosity must be between 0 and 1.")

    simulate_reads(
        args.length,
        args.coverage,
        args.read_length,
        args.error_rate,
        args.heterozygosity,
    )
