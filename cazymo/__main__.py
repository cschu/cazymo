# pylint: disable=C0103,C0301

""" module docstring """

import logging
import os
import pathlib
import subprocess
import sys

# pylint: disable=W0611
from gq.profilers import RegionQuantifier
from . import __version__
from .handle_args import handle_args


logger = logging.getLogger(__name__)


def main():

    args = handle_args(sys.argv[1:])

    logger.info("Version: %s", __version__)
    logger.info("Command: %s %s", os.path.basename(sys.argv[0]), " ".join(sys.argv[1:]))

    print(args)
    if args.input_files != "-" and not all(os.path.exists(f) for f in args.input_files):
        raise ValueError("bam file does not exist", args.bam_file)
    if not os.path.exists(args.annotation_db):
        raise ValueError("annotation database does not exist", args.annotation_db)

    if os.path.dirname(args.out_prefix):
        pathlib.Path(os.path.dirname(args.out_prefix)).mkdir(
            exist_ok=True, parents=True
        )

    fq = RegionQuantifier(
        db=args.annotation_db,
        out_prefix=args.out_prefix,
        ambig_mode=args.ambig_mode,
        strand_specific=args.strand_specific,
        calc_coverage=True,
        paired_end_count=args.paired_end_count,
        unmarked_orphans=args.unmarked_orphans,
        reference_type="domain",
    )

    ncpus = 8
    with subprocess.Popen(
        ("bwa", "mem", "-a", "-t", f"{ncpus}", "-K", "10000000", args.bwa_index, *args.input_files),
        stdout=subprocess.PIPE
    ) as bwa_proc:
        with subprocess.Popen(
            ("samtools", "view", "-F", "4", "-buSh", "-"),
            stdin=bwa_proc.stdout, stdout=subprocess.PIPE
        ) as samtools_proc:
            fq.process_bamfile(
                samtools_proc.stdout,
                aln_format=args.format,
                min_identity=args.min_identity, min_seqlen=args.min_seqlen
            )
    
    # "bwa mem -a -t ${align_cpus} ${blocksize} \$(readlink ${reference}) ${sample.id}_R1.sorted.fastq.gz ${reads2} | samtools view -F 4 -buSh - | ${sort_cmd}"

    # elif source.endswith(".bz2"):
	# 	bz2_pr = subprocess.Popen(("bzip2", "-dc", resolved_src), stdout=subprocess.PIPE)
	# 	with open(dest, "wt") as _out:
	# 		subprocess.run(("gzip", "-c", "-"), stdin=bz2_pr.stdout, stdout=_out)
	# else:

    # fq.process_bamfile(
    #     args.bam_file, aln_format=args.format, min_identity=args.min_identity, min_seqlen=args.min_seqlen
    # )


if __name__ == "__main__":
    main()
