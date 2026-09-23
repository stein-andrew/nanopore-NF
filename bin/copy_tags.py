#!/usr/bin/env python3
"""
Transfer per-read tags (base modifications, poly-A, move tables) from a dorado
uBAM onto an aligned BAM.

Why this exists
---------------
pbrun minimap2 drops MM/ML when reading a uBAM (--in-bam), and segfaults on the
--copy-comment route when the tags are large (8-mod RNA004 produces very long
comment lines). Since the GPU alignment is byte-identical to the CPU one, the
cheapest correct fix is to align without tags and graft them back afterwards.

Unlike `modkit repair`, this can carry any tag, which matters for RNA: Clair3-RNA
wants the move table (mv) and IsoQuant wants the poly-A estimate (pt).

Tags are copied verbatim from the source BAM -- no FASTQ text round-trip, so
nothing can be silently truncated or reformatted on the way through.

Method
------
Strip the uBAM down to the wanted tags, name-sort both files, then walk them in
lockstep. Constant memory: nothing is held in RAM but the current record.
Aligned read names must be a subset of the uBAM's, which is true whenever the
aligned file came from that uBAM.

Hard-clipped records (supplementary alignments without minimap2 -Y) carry a
truncated SEQ, so MM/ML/MN/mv copied verbatim would point past the end of it.
Those tags are skipped on such records and counted; everything else is copied.
Use `modkit repair` instead when supplementary alignments matter (the DNA arm).

Usage
-----
  copy_tags.py \\
      --ubam    sample.ubam.bam \\
      --aligned sample.aligned.bam \\
      --out     sample.tagged.bam \\
      --tags    MM,ML,MN,pt,qs,mv \\
      --threads 16
"""

import argparse
import os
import shutil
import subprocess
import sys
import tempfile

DEFAULT_TAGS = "MM,ML,MN,pt,qs"
# RG is deliberately absent: the aligner already wrote a read group into the
# aligned BAM, and a second RG tag on the record would conflict with the header.
# ts is absent too: dorado's ts:i (trimmed samples) collides with minimap2's
# ts:A (transcript strand), which IsoQuant reads.

# Tags whose values index into SEQ. Wrong on a hard-clipped record.
POSITIONAL = {"MM", "ML", "MN", "mv"}


def die(msg):
    sys.exit(f"FATAL: {msg}")


def sh(cmd):
    """Run a command, failing loudly. Used for the samtools sort passes."""
    proc = subprocess.run(cmd)
    if proc.returncode != 0:
        die(f"command failed ({proc.returncode}): {' '.join(cmd)}")


def name_sort(samtools, src, dst, threads, tmpdir, keep_tags=None, mem_per_thread="1G"):
    print(f"  name-sorting {os.path.basename(src)}", flush=True)
    sort = [samtools, "sort", "-n", "-@", str(threads), "-m", mem_per_thread,
            "-T", os.path.join(tmpdir, os.path.basename(dst) + ".tmp"), "-o", dst]
    if not keep_tags:
        sh(sort + [src])
        return
    # Drop everything but the wanted tags first: on an emit-moves uBAM the move
    # table is most of the file, and there is no reason to sort it if unwanted.
    view = subprocess.Popen([samtools, "view", "-u", "-@", str(threads),
                             "--keep-tag", keep_tags, src], stdout=subprocess.PIPE)
    srt = subprocess.run(sort + ["-"], stdin=view.stdout)
    view.stdout.close()
    if view.wait() != 0 or srt.returncode != 0:
        die(f"name sort of {src} failed (view {view.returncode}, sort {srt.returncode})")


def count_tagged(samtools, bam, tag, n=2000):
    """How many of the first `n` records carry `tag`. Reads only those records."""
    p = subprocess.Popen([samtools, "view", bam], stdout=subprocess.PIPE, text=True)
    found = 0
    for i, line in enumerate(p.stdout):
        if i >= n:
            break
        if f"\t{tag}:" in line:
            found += 1
    p.stdout.close()
    p.kill()
    p.wait()
    return found


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ubam", required=True, help="source uBAM carrying the tags")
    ap.add_argument("--aligned", required=True, help="aligned BAM missing the tags")
    ap.add_argument("--out", required=True, help="output BAM (coordinate sorted, indexed)")
    ap.add_argument("--tags", default=DEFAULT_TAGS,
                    help=f"comma-separated tags to copy (default: {DEFAULT_TAGS})")
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--samtools", default="samtools")
    ap.add_argument("--tmpdir", default=os.environ.get("TMPDIR", "/tmp"))
    args = ap.parse_args()

    samtools = shutil.which(args.samtools)
    if not samtools:
        die(f"samtools not found on PATH (looked for '{args.samtools}')")
    for f in (args.ubam, args.aligned):
        if not os.path.isfile(f):
            die(f"missing input: {f}")

    names = [t.strip() for t in args.tags.split(",") if t.strip()]
    if not names:
        die("--tags is empty")
    wanted = tuple(t + ":" for t in names)
    print(f"copying tags: {', '.join(names)}")

    workdir = tempfile.mkdtemp(prefix="tagxfer_", dir=args.tmpdir)
    try:
        ns_ubam = os.path.join(workdir, "ubam.ns.bam")
        ns_aln = os.path.join(workdir, "aln.ns.bam")
        name_sort(samtools, args.ubam, ns_ubam, args.threads, workdir, keep_tags=",".join(names))
        name_sort(samtools, args.aligned, ns_aln, args.threads, workdir)

        # Aligned side keeps its header (-h); tag side is records only.
        p_aln = subprocess.Popen([samtools, "view", "-h", "-@", str(args.threads), ns_aln],
                                 stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
        p_tag = subprocess.Popen([samtools, "view", "-@", str(args.threads), ns_ubam],
                                 stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
        # Stream SAM straight into the coordinate sort: no intermediate BAM on disk.
        p_out = subprocess.Popen([samtools, "sort", "-@", str(args.threads), "-m", "1G",
                                  "-T", os.path.join(workdir, "cs"), "--write-index",
                                  "-o", f"{args.out}##idx##{args.out}.bai", "-"],
                                 stdin=subprocess.PIPE, text=True, bufsize=1 << 20)

        cur_name = None
        cur_tags = []
        n_records = n_tagged = n_missing = n_clipped = 0

        def advance_to(name):
            """Pull uBAM records until we reach `name`. Returns its tags, or None."""
            nonlocal cur_name, cur_tags
            while cur_name != name:
                line = p_tag.stdout.readline()
                if not line:
                    return None
                f = line.rstrip("\n").split("\t")
                cur_name = f[0]
                cur_tags = [t for t in f[11:] if t.startswith(wanted)]
            return cur_tags

        for line in p_aln.stdout:
            if line.startswith("@"):
                p_out.stdin.write(line)
                continue
            n_records += 1
            fields = line.rstrip("\n").split("\t")
            tags = advance_to(fields[0])
            if tags is None:
                die(f"read '{fields[0]}' not found in uBAM -- are these the same sample? "
                    "(name sort order must come from the same samtools version)")
            if "H" in fields[5]:
                kept = [t for t in tags if t[:2] not in POSITIONAL]
                if len(kept) != len(tags):
                    n_clipped += 1
                tags = kept
            # Don't duplicate a tag the aligner already emitted.
            have = {t[:2] for t in fields[11:]}
            add = [t for t in tags if t[:2] not in have]
            if add:
                fields.extend(add)
                n_tagged += 1
            else:
                n_missing += 1
            p_out.stdin.write("\t".join(fields) + "\n")

        print(f"  records: {n_records}  tagged: {n_tagged}  without tags: {n_missing}  "
              f"hard-clipped (positional tags skipped): {n_clipped}", flush=True)
        print("  coordinate sorting", flush=True)
        p_out.stdin.close()
        for name, p in (("aligned reader", p_aln), ("tag reader", p_tag), ("sorter", p_out)):
            if p.wait() != 0:
                die(f"{name} exited {p.returncode}")
        if n_tagged == 0:
            die("nothing was tagged -- check --tags matches what is in the uBAM")
    finally:
        shutil.rmtree(workdir, ignore_errors=True)

    # Verify rather than assume. Cheap, and the whole point of the exercise.
    sh([samtools, "quickcheck", "-v", args.out])
    for tag in ("MM", "mv"):
        if tag in names:
            found = count_tagged(samtools, args.out, tag)
            print(f"  {tag} present in {found} of the first 2000 output records")
            if found == 0:
                die(f"output has no {tag} tags")

    print(f"done: {args.out}")


if __name__ == "__main__":
    main()
