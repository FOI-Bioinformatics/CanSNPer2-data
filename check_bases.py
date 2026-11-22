#!/usr/bin/env python3
"""Check the SNP bases of the schemes against their reference genomes.

A SNP is placed by the name of a sequence of its reference genome (a chromosome
or plasmid, named by its GenBank accession) and a 1-based position on that
sequence. validate.py checks the scheme files without the genomes; this script
downloads the GenBank assembly of every reference genome (once, into a cache
folder) and checks that the base of the named sequence at each SNP position is
the derived or the ancestral base of that SNP, in <name>_snp.txt and in
<name>_canSNPs.txt. The order of the sequences in the assembly file does not
matter. A sequence name that is not in the assembly, or positions on another
genome, are reported as problems.

    python check_bases.py
    python check_bases.py --scheme brucella_classical --cache ~/.cache/cansnper2-data

A SNP listed under known_issues in schemes.toml is reported as a warning. Exits
non-zero if any scheme has a problem.
"""

import argparse
import gzip
import os
import shutil
import sys
import tempfile
from urllib.request import urlopen

import validate

NCBI = "https://ftp.ncbi.nlm.nih.gov/genomes/all"
TIMEOUT = 120  # seconds per connection


def genome_url(accession, assembly):
    digits = accession[4:].split(".")[0]
    prefix = accession[:3]
    stem = f"{accession}_{assembly}"
    return f"{NCBI}/{prefix}/{digits[0:3]}/{digits[3:6]}/{digits[6:9]}/{stem}/{stem}_genomic.fna.gz"


def download(url, target):
    with urlopen(url, timeout=TIMEOUT) as response, open(target, "wb") as out:
        shutil.copyfileobj(response, out)


def fetch(cache, accession, assembly, retrieve=download):
    """Path of the gzipped genome in cache, downloaded when it is not there yet.

    The download goes to a temporary file that is only renamed once gzip has
    read it to the end, so a transfer cut short is not taken for a genome by a
    later run."""
    path = os.path.join(cache, f"{accession}_{assembly}_genomic.fna.gz")
    if os.path.exists(path):
        return path
    handle, partial = tempfile.mkstemp(dir=cache, suffix=".part")
    os.close(handle)
    try:
        retrieve(genome_url(accession, assembly), partial)
        with gzip.open(partial, "rb") as stream:
            while stream.read(1 << 20):
                pass
        os.replace(partial, path)
    finally:
        if os.path.exists(partial):
            os.remove(partial)
    return path


def read_fasta(path):
    """[(name, sequence), ...] in file order, sequences in upper case."""
    records, name, parts = [], None, []
    with gzip.open(path, "rt") as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(parts).upper()))
                words = line[1:].split()
                name, parts = (words[0] if words else ""), []
            else:
                parts.append("".join(line.split()))
    if name is not None:
        records.append((name, "".join(parts).upper()))
    return records


def scheme_snps(scheme):
    """The canSNPs of a scheme as (file, line, label, node, genome, sequence, position,
    derived, ancestral). Every canSNP is checked, the alternatives of a node included."""
    snps = []
    _raw, _lines, rows = validate.read_rows(scheme.path(scheme.cansnps))
    for number, row in enumerate(rows[1:], start=2):
        if len(row) == len(validate.CANSNP_HEADER):
            node, rank, _reference, genome, sequence, position, derived, ancestral, _score, _typed = row
            label = f"canSNP {node}" if rank == "1" else f"canSNP {node} (rank {rank})"
            snps.append((scheme.cansnps, number, label, node, genome, sequence,
                         position, derived, ancestral))
    return snps


def check_scheme(scheme, cache, retrieve=download):
    """Return the problems, the warnings (known issues) and a report of the sequences
    of each reference genome with the number of SNPs on each."""
    problems, warnings, report = [], [], []
    _raw, _lines, rows = validate.read_rows(scheme.path(scheme.references))
    references = {}
    for row in rows[1:]:
        if len(row) == len(validate.REFERENCE_HEADER):
            genome, _strain, genbank, _refseq, assembly = row
            references[genome] = (genbank, assembly)

    by_genome = {}
    for snp in scheme_snps(scheme):
        by_genome.setdefault(snp[4], []).append(snp)

    for genome in sorted(by_genome):
        if genome not in references:
            continue  # validate.py reports a SNP on a genome without a reference
        genbank, assembly = references[genome]
        try:
            records = read_fasta(fetch(cache, genbank, assembly, retrieve))
        except Exception as error:  # a download or gzip failure is reported, not raised
            problems.append(f"{scheme.references}: could not get {genbank} ({assembly}) "
                            f"for genome '{genome}': {error}")
            continue
        sequences = dict(records)
        on_sequence = {name: 0 for name, _seq in records}

        for relative, number, label, snp_id, _genome, sequence, position, derived, ancestral in by_genome[genome]:
            if not validate.POSITION_RE.match(position) or int(position) < 1:
                continue  # validate.py reports it
            pos = int(position)
            report_to = warnings if snp_id in scheme.known_issues else problems
            where = f"{relative}:{number}: {label} at {sequence} position {pos} of {genome}"
            if sequence not in sequences:
                report_to.append(f"{where}: the reference {genbank} has no sequence "
                                 f"'{sequence}'; its sequences are "
                                 f"{', '.join(name for name, _seq in records)}")
                continue
            if pos > len(sequences[sequence]):
                report_to.append(f"{where}: the sequence {sequence} of {genbank} has "
                                 f"{len(sequences[sequence])} bases")
                continue
            on_sequence[sequence] += 1
            base = sequences[sequence][pos - 1]
            if base not in (derived, ancestral):
                report_to.append(f"{where}: the reference {genbank} has {base}, neither "
                                 f"the derived ({derived}) nor the ancestral ({ancestral}) base")
        listed = ", ".join(f"{name} {len(seq)} bp ({on_sequence[name]} SNP"
                           f"{'s' if on_sequence[name] != 1 else ''})" for name, seq in records)
        report.append(f"{genome}: {listed}")
    return problems, warnings, report


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--scheme", action="append", metavar="NAME",
                        help="check only this scheme; repeatable")
    parser.add_argument("--cache", default=os.path.join(validate.REPO, ".references"),
                        help="folder for the downloaded genomes (default: .references)")
    args = parser.parse_args()

    os.makedirs(args.cache, exist_ok=True)
    schemes = [s for s in validate.discover(args.scheme) if s.kind == "snp"]
    failed = warned = 0
    for scheme in schemes:
        problems, warnings, report = check_scheme(scheme, args.cache)
        if problems:
            failed += 1
            print(f"FAIL {scheme.name} ({len(problems)} problem{'s' if len(problems) != 1 else ''})")
        else:
            print(f"ok   {scheme.name}")
        for line in report:
            print(f"     {line}")
        for problem in problems:
            print(f"  - {problem}")
        for warning in warnings:
            print(f"  ! known issue: {warning}")
        warned += len(warnings)

    if warned:
        print(f"\n{warned} known issue(s) reported; see known_issues in schemes.toml.")
    if failed:
        print(f"{failed} of {len(schemes)} schemes have SNP bases that do not match their reference.")
        return 1
    print(f"All {len(schemes)} schemes match their reference genomes.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
