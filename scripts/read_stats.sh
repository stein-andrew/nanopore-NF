#!/usr/bin/env bash
# Read-level accuracy of one aligned direct-RNA BAM, as in Fig. 1 of Ghohabi Esfahani et al.
# (Genome Biology 2026): per-read identity and indel rates, and the base-substitution matrix.
#
#   read_stats.sh <bam> <reference.fa> <out prefix> [fraction 0.02] [region chr22] [threads 4]
#
#   <prefix>.read_identity.tsv.gz  one row per sampled primary read with >= 200 aligned bases:
#                                  aligned_len (M+D, reference bases), matches, mismatches,
#                                  insertions, deletions, identity, ins_rate, del_rate
#                                  identity = matches / (matches + mismatches + insertions + deletions)
#   <prefix>.substitutions.tsv     ref_base read_base count, in RNA sense (U for T), from a pileup
#                                  of <region> at positions with >= 20 reads on that strand
#
# Reads are a random <fraction> of primary alignments with MAPQ >= MIN_MAPQ (default 40, the same
# reads Stage 2 uses everywhere; low-MAPQ reads include rRNA placed on rRNA-like loci at ~80%
# identity, which would measure misplacement, not basecalling). Mismatches come
# from NM (edit distance) minus inserted and deleted bases; NM is written by minimap2. The read's
# strand is the RNA's strand in direct RNA, so reverse-strand bases are complemented to RNA sense.
set -euo pipefail
die() { echo "read_stats.sh: $*" >&2; exit 1; }
(( $# >= 3 )) || die "usage: $0 <bam> <reference.fa> <out prefix> [fraction] [region] [threads]"
bam=$1 ref=$2 prefix=$3 frac=${4:-0.02} region=${5:-chr22} threads=${6:-4} minq=${MIN_MAPQ:-40}
[[ -s $bam && -s $ref ]] || die "missing $bam or $ref"
seed_frac="42${frac#0}"                       # samtools -s SEED.FRACTION

samtools view -@ "$threads" -F 2308 -q "$minq" -s "$seed_frac" "$bam" | awk -F'\t' -v OFS='\t' '
    BEGIN { print "aligned_len", "matches", "mismatches", "insertions", "deletions", "identity", "ins_rate", "del_rate" }
    { m = 0; ins = 0; del = 0; nm = -1
      n = split($6, L, /[MIDNSHP=X]/); split($6, O, /[0-9]+/)
      for (i = 1; i < n; i++) { op = O[i + 1]
          if (op == "M" || op == "=" || op == "X") m += L[i]; else if (op == "I") ins += L[i]; else if (op == "D") del += L[i] }
      for (i = 12; i <= NF; i++) if (substr($i, 1, 5) == "NM:i:") { nm = substr($i, 6) + 0; break }
      if (nm < 0 || m + del < 200) next
      mm = nm - ins - del; if (mm < 0) mm = 0
      tot = m + ins + del
      printf "%d\t%d\t%d\t%d\t%d\t%.5f\t%.5f\t%.5f\n", m + del, m - mm, mm, ins, del, (m - mm) / tot, ins / tot, del / tot }
' | gzip > "$prefix.read_identity.tsv.gz.partial"
mv "$prefix.read_identity.tsv.gz.partial" "$prefix.read_identity.tsv.gz"
echo "read_stats.sh: $(( $(zcat "$prefix.read_identity.tsv.gz" | wc -l) - 1 )) reads in $prefix.read_identity.tsv.gz" >&2

samtools mpileup -r "$region" -f "$ref" -q "$minq" -Q 0 -B -d 0 --ff 2308 "$bam" 2> /dev/null | python3 -c '
import re, sys
comp = {"A": "U", "C": "G", "G": "C", "T": "A"}          # reverse-strand reference base, in RNA sense
sense = {"A": "A", "C": "C", "G": "G", "T": "U"}
rev_read = {"a": "U", "c": "G", "g": "C", "t": "A"}
fwd_read = {"A": "A", "C": "C", "G": "G", "T": "U"}
bases = "ACGU"
counts = {(r, b): 0 for r in bases for b in bases}
skip_mark = re.compile(r"\^.|\$")
indel = re.compile(r"[+-](\d+)")
for line in sys.stdin:
    f = line.rstrip("\n").split("\t")
    if len(f) < 5 or f[2].upper() not in sense:
        continue
    refb, s = f[2].upper(), skip_mark.sub("", f[4])
    out, i = [], 0                                         # drop inserted / deleted sequences
    for m in indel.finditer(s):
        if m.start() < i:
            continue
        out.append(s[i:m.start()]); i = m.end() + int(m.group(1))
    out.append(s[i:]); s = "".join(out)
    fwd = {"A": s.count("A"), "C": s.count("C"), "G": s.count("G"), "T": s.count("T")}
    rev = {"a": s.count("a"), "c": s.count("c"), "g": s.count("g"), "t": s.count("t")}
    nf, nr = s.count("."), s.count(",")
    if nf + sum(fwd.values()) >= 20:                       # forward strand: RNA sense as written
        r = sense[refb]; counts[(r, r)] += nf
        for b, n in fwd.items(): counts[(r, fwd_read[b])] += n
    if nr + sum(rev.values()) >= 20:                       # reverse strand: complement both
        r = comp[refb]; counts[(r, r)] += nr
        for b, n in rev.items(): counts[(r, rev_read[b])] += n
print("ref_base\tread_base\tcount")
for r in bases:
    for b in bases:
        print(f"{r}\t{b}\t{counts[(r, b)]}")
' > "$prefix.substitutions.tsv.partial"
mv "$prefix.substitutions.tsv.partial" "$prefix.substitutions.tsv"
echo "read_stats.sh: substitutions over $region in $prefix.substitutions.tsv" >&2