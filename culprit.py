#!/usr/bin/env python3
from __future__ import annotations
import argparse
import math
from pathlib import Path
from typing import Dict, Tuple

# import your align() implementation
from alignment import align

KNOWN_CODES = {
    "hg38":   "Human (Homo sapiens)",
    "panTro4":"Chimp (Pan troglodytes)",
    "rheMac3":"Rhesus macaque (Macaca mulatta)",
    "canFam3":"Dog (Canis lupus familiaris)",
    "rn5":    "Rat (Rattus norvegicus)",
    "mm10":   "Mouse (Mus musculus)",
    "unknown":"Unknown",
}

def parse_fasta_by_species(path: Path) -> Dict[str, str]:
    """
    Parses a FASTA-like file where header lines begin with '>' and contain
    a token like uc002tuu.1_hg38_8_17 ... or unknown.1_unknown_8_17
    Returns dict: species_code -> sequence (concatenated, uppercased, no spaces).
    """
    species_to_seq: Dict[str, str] = {}
    cur_species = None
    with path.open("r", encoding="utf-8") as f:
        for raw in f:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                token = line[1:].split()[0]  # take first token on header line
                # Try to find a known species code inside the token
                found = None
                for code in KNOWN_CODES.keys():
                    if code in token:
                        found = code
                        break
                if found is None:
                    # fallback: header like "unknown.1_unknown_8_17"
                    if "unknown" in token:
                        found = "unknown"
                    else:
                        raise ValueError(f"Could not identify species from header: {line}")
                cur_species = found
                species_to_seq.setdefault(cur_species, "")
            else:
                if cur_species is None:
                    raise ValueError("Sequence line encountered before any FASTA header.")
                species_to_seq[cur_species] += line.strip().upper()
    return species_to_seq

def main():
    p = argparse.ArgumentParser(description="Compare unknown DNA against known species with Needleman–Wunsch.")
    p.add_argument("fasta", type=Path, help="Path to lct_exon8.txt (FASTA-like) file.")
    p.add_argument("--band", type=int, default=-1,
                   help="Banded width d (use -1 for unrestricted). For the project core, d=3.")
    p.add_argument("--match", type=int, default=-3, help="Match award (default -3).")
    p.add_argument("--indel", type=int, default=5, help="Indel penalty (default 5).")
    p.add_argument("--sub", type=int, default=1, help="Substitution penalty (default 1).")
    args = p.parse_args()

    data = parse_fasta_by_species(args.fasta)

    if "unknown" not in data:
        raise SystemExit("No 'unknown' sequence found in the FASTA file.")

    unknown_seq = data["unknown"]

    # collect scores
    results: list[Tuple[str, float]] = []
    for code, seq in data.items():
        if code == "unknown":
            continue
        score, a1, a2 = align(
            seq1=seq,
            seq2=unknown_seq,
            match_award=args.match,
            indel_penalty=args.indel,
            sub_penalty=args.sub,
            banded_width=args.band
        )
        results.append((code, score))

    # sort by score (lower is better since we're minimizing cost)
    results.sort(key=lambda x: (math.inf if x[1] is None else x[1]))

    # pretty print
    print("\nDNA similarity to UNKNOWN (lower score = more similar):")
    print(f"Scoring: match={args.match}, indel={args.indel}, sub={args.sub}, band={args.band}")
    print("-" * 72)
    rank = 1
    for code, score in results:
        label = KNOWN_CODES.get(code, code)
        if score is math.inf:
            print(f"{rank:>2}. {label:<35}  score = inf  (not alignable with given band)")
        else:
            print(f"{rank:>2}. {label:<35}  score = {score}")
        rank += 1

    # winner (best match)
    best_code, best_score = results[0]
    best_label = KNOWN_CODES.get(best_code, best_code)
    if best_score is math.inf:
        print("\nBest: none — unknown is outside the band for all species.")
    else:
        print(f"\nBest match: {best_label}  (score {best_score})")

if __name__ == "__main__":
    main()
