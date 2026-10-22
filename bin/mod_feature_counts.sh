#!/usr/bin/env bash
# Count `modkit extract calls` rows per site and call code: over all reads, per isoform and
# per gene of each read. Reads the calls from stdin as they stream out of modkit, so the
# per-read table (sites x depth rows) is never written anywhere.
#
#   modkit extract calls S.bam - --include-bed positions.bed ... | mod_feature_counts.sh S.read_features.tsv > counts.tsv
#
# read_features.tsv  read_id, isoform, gene ('.' = none), from bin/read_features.py
# counts.tsv         level, feature, chrom, pos, strand, base, call_code, n   (passing calls only)
#   base             the read's base at the site (canonical_base): a read whose base differs
#                    from the reference's is a different-base read in modkit pileup (Ndiff),
#                    not valid coverage, so downstream only calls on the mod's own base count
#   level 'all'      every passing call, feature '.': checked against the pileup downstream
#   level 'isoform'  reads with a unique isoform assignment
#   level 'gene'     reads with a single gene
# call_code '-' is a canonical call; it counts toward valid coverage, like in modkit pileup.
#
# Columns are found by header name. The strand is ref_mod_strand (the modified base's strand
# on the reference) when present, else ref_strand. A fail value other than true/false, or a
# missing column, is an error rather than a guess.
set -euo pipefail

feats=${1:?usage: modkit extract calls ... | $0 <read_features.tsv> > counts.tsv}
[[ -r $feats ]] || { echo "mod_feature_counts.sh: cannot read $feats" >&2; exit 1; }

# POSIX awk only: the image's awk is Ubuntu's mawk.
awk -F'\t' -v feats="$feats" '
BEGIN {
    OFS = "\t"
    while ((getline line < feats) > 0) {
        if (split(line, f, "\t") < 3) continue
        if (f[2] != ".") iso[f[1]] = f[2]
        if (f[3] != ".") gene[f[1]] = f[3]
    }
    close(feats)
}
NR == 1 {
    for (i = 1; i <= NF; i++) col[$i] = i
    split("read_id chrom ref_position call_code fail canonical_base", need, " ")
    for (k = 1; k <= 6; k++)
        if (!(need[k] in col)) { msg = "extract calls output has no " need[k] " column"; exit }
    if ("ref_mod_strand" in col) sc = col["ref_mod_strand"]
    else if ("ref_strand" in col) sc = col["ref_strand"]
    else { msg = "extract calls output has no ref_mod_strand or ref_strand column"; exit }
    r = col["read_id"]; ch = col["chrom"]; p = col["ref_position"]; cc = col["call_code"]; fl = col["fail"]
    cb = col["canonical_base"]
    next
}
{
    if ($fl == "true") { failed++; next }
    if ($fl != "false") { msg = "unrecognised fail value \"" $fl "\" on line " NR; exit }
    k = $ch OFS $p OFS $sc OFS $cb OFS $cc
    all[k]++
    if ($r in iso)  n_iso[iso[$r] OFS k]++
    if ($r in gene) n_gene[gene[$r] OFS k]++
    passed++
}
END {
    if (msg != "") { print "mod_feature_counts.sh: " msg > "/dev/stderr"; exit 1 }
    if (NR == 0)   { print "mod_feature_counts.sh: no input at all (not even a header)" > "/dev/stderr"; exit 1 }
    for (k in all)    print "all", ".", k, all[k]
    for (k in n_iso)  print "isoform", k, n_iso[k]
    for (k in n_gene) print "gene", k, n_gene[k]
    printf "mod_feature_counts.sh: %d passing calls, %d failing\n", passed, failed > "/dev/stderr"
}'