# pylint: disable=C0301,C0103
""" module docstring """

import argparse
import logging
import textwrap

from . import __version__


def handle_args(args):

    log_ap = argparse.ArgumentParser(prog="cazymo", add_help=False)
    log_ap.add_argument("-l", "--log_level", type=int, choices=range(1, 5), default=logging.INFO)
    log_args, _ = log_ap.parse_known_args(args)

    try:
        logging.basicConfig(
            level=log_args.log_level,
            format='[%(asctime)s] %(message)s'
        )
    except ValueError as invalid_loglevel_err:
        raise ValueError(f"Invalid log level: {log_args.log_level}") from invalid_loglevel_err

    ap = argparse.ArgumentParser(
        prog="cazymo",
        formatter_class=argparse.RawTextHelpFormatter,
        parents=(log_ap,),
    )
    ap.add_argument(
        "annotation_db",
        type=str,
        help=textwrap.dedent(
            """\
            Path to an sqlite3 database containing the reference annotation.
			"""
        ),
    )
    ap.add_argument(
        "bwa_index",
        type=str, help="",
    )
    ap.add_argument(
        "input_files",
        type=str,
        nargs="*",
        help=textwrap.dedent(
            """\
            Path to metagenomic reads in fastq format or alignments in sam or bam format.
            Fastq files can be supplied as a single unpaired file or two paired-end files.
            Sam or bam files need to be name-sorted and need to be specified via the --format flag.
            Ambiguous alignments need to be flagged as secondary
            alignments with the same read id as their primary alignment.
            (e.g. output from BWA mem -a).
            Input from STDIN can be used with '-'."""
        ),
    )
    ap.add_argument(
        "--out_prefix",
        "-o",
        type=str,
        default="cazymo",
        help="Prefix for output files.",
    )
    # ap.add_argument(
    #     "--ambig_mode",
    #     type=str,
    #     choices=("unique_only", "all1", "primary_only", "1overN"),
    #     default="unique_only",
    #     help=textwrap.dedent(
    #         """\
    #         Determines how ambiguous alignments should be treated. This setting mimics NGLess' behaviour.
    #         - 'unique_only' ignores any alignment flagged as ambiguous (MAPQ=0). This is the default setting.
    #         - 'all1' treats each alignment as unique (each ambiguous alignment contributes 1 count to features it aligns to.)
    #         - 'primary_only' takes the unique alignments and the primary and alignment of each ambiguous read group.
    #         - '1overN' each alignment contributes 1/(n=number of ambiguous alignments of the same read) counts to features it aligns to."""
    #     ),
    # )

    ap.add_argument(
        "--strand_specific",
        action="store_true",
        help="Perform strand-specific counting for RNAseq reads. "
        "This flag is currently ignored for paired-end data.",
    )

    ap.add_argument(
        "--min_identity",
        type=float,
        default=0.97,
        help="Minimum sequence identity [n_match/length] for an alignment to be considered.",
    )

    ap.add_argument(
        "--min_seqlen",
        type=int,
        default=45,
        help="Minimum read length [bp] for an alignment to be considered.",
    )

    ap.add_argument(
        "--format",
        type=str,
        choices=("fastq", "fq", "FASTQ", "FQ", "sam", "bam", "SAM", "BAM"),
        default="fastq",
        help="Format of the alignment input. Supported: fastq, sam, bam. Fastq files can be gzipped.",
    )

    # ap.add_argument(
    #     "--paired_end_count",
    #     type=int,
    #     choices=(1, 2),
    #     default=1,
    #     help="Paired-end count contribution: 0.5 / mate (1) or 1 / mate (2) [1]",
    # )

    # orphan reads will not have flag 0x1 set
    ap.add_argument(
        "--unmarked_orphans",
        action="store_true",
        help="Ensure that alignments from unmarked orphan reads (from preprocessing) are properly accounted for.",
    )

    ap.add_argument(
        "--no_prefilter",
        action="store_true",
        help="",
    )

    ap.add_argument(
        "--cpus_for_alignment", "-t",
        type=int, default=1,
        help="",
    )

    ap.add_argument(
        "--version", "-v", action="version", version="%(prog)s " + __version__
    )
    ap.add_argument("--debug", action="store_true")

    return ap.parse_args(args)
