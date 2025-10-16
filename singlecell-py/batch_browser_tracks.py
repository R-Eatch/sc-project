#!/usr/bin/env python3
"""
Batch browser‑track drawing with automatic BW / narrowPeak pairing (v2.4)
--------------------------------------------------------------------------
Fixes & improvements
-------------------
* **Gene name parsing**: always prefer `gene_name` if present; fallback to
  `Name`, then `gene_id` / `ID`. This fixes missing EPCAM/KRT8/LALBA etc.
* Support `--gene_case_insensitive` and `--debug_gene` for troubleshooting.
* Other features unchanged: automatic BW/peak pairing, optional gene track,
  customizable heights, hide peak labels, etc.

Dependencies
------------
* Python ≥3.8
* pyGenomeTracks ≥3.7
"""
import argparse, glob, subprocess, sys, re
from collections import defaultdict
from pathlib import Path

# ---------- GTF parsing helpers ---------- #
_ATTR_RE = re.compile(r'\s*([^\s]+)\s+"?([^;\"]+)"?\s*')

PRIORITY_KEYS = ("gene_name", "name", "gene")  # preferred
FALLBACK_KEYS = ("gene_id", "id")                 # fallback if no name


def extract_attr(attrs: str, case_insensitive: bool = False):
    """Return chosen gene key according to priority list"""
    chosen = {}
    for field in attrs.split(';'):
        field = field.strip()
        if not field:
            continue
        m = _ATTR_RE.match(field)
        if not m:
            continue
        key, val = m.groups()
        key_lc = key.lower()
        if case_insensitive:
            val = val.upper()
        if key_lc in PRIORITY_KEYS:
            # first preferred key wins
            return val
        if key_lc in FALLBACK_KEYS and 'fallback' not in chosen:
            chosen['fallback'] = val
    return chosen.get('fallback')  # may be None


def build_gene_index(gtf_file: Path, case_insensitive: bool = False):
    """Return dict: gene_name -> (chrom,start,end,strand)"""
    idx = {}
    with gtf_file.open() as fh:
        for ln in fh:
            if ln.startswith('#'):
                continue
            chrom, _, feat, s, e, _, strand, _, attrs = ln.rstrip().split('\t')
            if feat != 'gene':
                continue
            name = extract_attr(attrs, case_insensitive)
            if name:
                idx[name] = (chrom, int(s), int(e), strand)
    return idx


def tss_region(pos, flank):
    chrom, start, end, strand = pos
    tss = start if strand == '+' else end
    return f"{chrom}:{max(0, tss - flank)}-{tss + flank}"


def run(cmd):
    print(' '.join(map(str, cmd)))
    subprocess.run(cmd, check=True)


# ---------- main ---------- #

def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument('--gtf', required=True, help='Genome GTF file')
    ap.add_argument('--bw', nargs='+', required=True, help='bigWig file globs')
    ap.add_argument('--peaks', nargs='+', help='narrowPeak/BED file globs')
    ap.add_argument('--genes', required=True, help='Gene list file or comma list')
    ap.add_argument('--gene_case_insensitive', action='store_true')
    ap.add_argument('--debug_gene', help='Print GTF line(s) for this gene and exit')
    ap.add_argument('--gtf_track', help='Add gene annotation track')
    ap.add_argument('--flank', type=int, default=20000)
    ap.add_argument('--bw_height', type=float, default=2.0)
    ap.add_argument('--peak_height', type=float, default=0.6)
    ap.add_argument('--hide_peak_labels', action='store_true')
    ap.add_argument('--dpi', type=int, default=300)
    ap.add_argument('--outdir', default='tracks_out')
    args = ap.parse_args(argv)

    outdir = Path(args.outdir).expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    # Collect files
    bw_files = [Path(f).resolve() for pat in args.bw for f in glob.glob(pat)]
    peak_files = [Path(f).resolve() for pat in (args.peaks or []) for f in glob.glob(pat)]
    if not bw_files:
        sys.exit('[ERROR] No bigWig files matched')

    # Group by stem label
    label_map = defaultdict(dict)
    for f in bw_files:
        label_map[f.stem]['bw'] = f
    for f in peak_files:
        label_map[f.stem]['pk'] = f

    # Compose tracks.ini
    ini_lines = []
    for label in sorted(label_map):
        ent = label_map[label]
        if 'bw' in ent:
            ini_lines += [
                f'[{label}_bw]',
                f'file = {ent["bw"]}',
                f'title = {label}',
                f'height = {args.bw_height}',
                'color = grey',
                'plot_type = fill',
                ''
            ]
        if 'pk' in ent:
            ini_lines += [
                f'[{label}_pk]',
                f'file = {ent["pk"]}',
                'track_type = narrowPeak',
                f'height = {args.peak_height}',
            ]
            if args.hide_peak_labels:
                ini_lines.append('display_labels = no')
            ini_lines.append('')

    if args.gtf_track:
        ini_lines += [
            '[Genes]',
            f'file = {args.gtf_track}',
            'track_type = gene',
            'height = 2',
            'color = darkblue',
            'display = stacked',
            'gene_rows = 3',
            ''
        ]

    ini_path = outdir / 'tracks.ini'
    ini_path.write_text('\n'.join(ini_lines))
    print(f'[INFO] tracks.ini written to {ini_path}')

    # Build gene index
    gene_idx = build_gene_index(Path(args.gtf).resolve(), args.gene_case_insensitive)
    print(f'[INFO] Parsed {len(gene_idx):,} genes from GTF')

    if args.debug_gene:
        key = args.debug_gene.upper() if args.gene_case_insensitive else args.debug_gene
        if key in gene_idx:
            print('[DEBUG] Gene', args.debug_gene, 'found at', gene_idx[key])
        else:
            print('[DEBUG] Gene', args.debug_gene, 'NOT found')
        sys.exit(0)

    # Prepare gene list
    if Path(args.genes).is_file():
        genes = [l.strip() for l in Path(args.genes).read_text().splitlines() if l.strip()]
    else:
        genes = [g.strip() for g in args.genes.split(',') if g.strip()]
    if args.gene_case_insensitive:
        genes = [g.upper() for g in genes]

    # Plot per gene
    for g in genes:
        if g not in gene_idx:
            print(f'[WARN] {g} not in GTF')
            continue
        region = tss_region(gene_idx[g], args.flank)
        out_pdf = outdir / f'{g}.pdf'
        run(['pyGenomeTracks', '--tracks', str(ini_path), '--region', region,
             '--outFileName', str(out_pdf), '--dpi', str(args.dpi), '--fontSize', '10'])

    print('[DONE] All browser tracks saved to', outdir)


if __name__ == '__main__':
    main()

