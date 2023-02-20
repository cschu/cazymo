# pylint: disable=C0103,C0301

""" module docstring """

import logging
import os
import pathlib
import subprocess
import sys

# pylint: disable=W0611
from gq.profilers import RegionQuantifier
from gq.db.db_import import DomainBedDatabaseImporter

from . import __version__
from .handle_args import handle_args


logger = logging.getLogger(__name__)


def check_bwa_index(prefix):
    suffixes = (".amb", ".ann", ".bwt", ".pac", ".sa")
    return all(os.path.isfile(prefix + suffix) for suffix in suffixes)


# pylint: disable=R0913
def run_alignment(
    profiler,
    input_files,
    bwa_index,
    annotation_db,
    out_prefix,
    cpus_for_alignment=1,
    no_prefilter=False,
    min_identity=None,
    min_seqlen=None,
    unmarked_orphans=False,
):
    samtools_io_flags = "-buSh" if no_prefilter else "-Sh"

    commands = [
        f"bwa mem -v 1 -a -t {cpus_for_alignment} -K 10000000 {bwa_index} {' '.join(input_files)}",
        f"read_count {out_prefix} --all",
        f"samtools view -F 4 {samtools_io_flags} -",
    ]

    if not no_prefilter:
        logging.info("Prefiltering activated.")
        commands += [
            f"read_count {out_prefix}",
            "samtools view -buSh -",
            f"bedtools intersect -u -ubam -a stdin -b {annotation_db}",
        ]

    logger.info("Used command: %s", " | ".join(commands))

    try:
        with subprocess.Popen(" | ".join(commands), shell=True, stdout=subprocess.PIPE) as read_processing_proc:
            profiler.count_alignments(
                read_processing_proc.stdout,
                aln_format="bam",
                min_identity=min_identity,
                min_seqlen=min_seqlen,
                external_readcounts=None if no_prefilter else (out_prefix + ".readcount.json"),
                unmarked_orphans=unmarked_orphans,
            )
    except Exception as err:
        logger.error("Caught some exception:")
        logger.error("%s", err)
        raise Exception from err


def check_input_reads(fwd=None, rev=None, singles=None, orphans=None):
    fwd_reads = fwd.split(",") if fwd else None
    rev_reads = rev.split(",") if rev else None
    single_reads = singles.split(",") if singles else None
    orphan_reads = orphans.split(",") if orphans else None

    all_readsets = []

    if fwd_reads and rev_reads:
        if len(fwd_reads) == len(rev_reads):
            all_readsets += zip((["paired"] * len(fwd_reads)), fwd_reads, rev_reads)
        else:
            raise ValueError(f"Found different numbers of forward/R1 {len(fwd_reads)} and reverse/R2 {len(rev_reads)} reads.")
    elif fwd_reads:
            logger.warning("Found -1 forward/R1 reads but no -2 reverse/R2 reads. Treating these as single-end reads.")
            all_readsets += zip((["single"] * len(fwd_reads)), fwd_reads)
    elif rev_reads:
        raise ValueError(f"Found -2 reverse/R2 reads but no -1 forward/R1 reads.")
    
    
    if single_reads:
        all_readsets += zip((["single"] * len(single_reads)), single_reads)
    if orphan_reads:
        all_readsets += zip((["orphan"] * len(orphan_reads)), orphan_reads)

    if not all_readsets:
        raise ValueError(f"No input reads specified.")

    for _, *reads in all_readsets:
        for r in reads:
            if not os.path.isfile(r):
                raise ValueError(f"{r} does not seem to be a valid read file.")

    return all_readsets


def main():

    args = handle_args(sys.argv[1:])

    logger.info("Version: %s", __version__)
    logger.info("Command: %s %s", os.path.basename(sys.argv[0]), " ".join(sys.argv[1:]))

    print(args)

    input_data = check_input_reads(args.reads1, args.reads2, args.singles, args.orphans)

    # if args.input_files != "-" and not all(os.path.exists(f) for f in args.input_files):
    #     input_files_str = "\n".join(args.input_files)
    #     raise ValueError(f"There is an issue with your input files. Please check.\n{input_files_str}")
    if not os.path.exists(args.annotation_db):
        raise ValueError(f"{args.annotation_db} is not a valid annotation database")
    if not check_bwa_index(args.bwa_index):
        raise ValueError(f"{args.bwa_index} is not a valid bwa index.")

    if os.path.dirname(args.out_prefix):
        pathlib.Path(os.path.dirname(args.out_prefix)).mkdir(
            exist_ok=True, parents=True
        )

    

    db_importer = DomainBedDatabaseImporter(logger, args.annotation_db, single_category="cazy")
    logger.info("Finished loading database.")

    fq = RegionQuantifier(
        db=db_importer,
        out_prefix=args.out_prefix,
        ambig_mode="1overN",
        strand_specific=args.strand_specific,
        reference_type="domain",
    )

    for input_type, *reads in input_data:

        logger.info("Running %s alignment: %s", input_type, ",".join(reads))

        run_alignment(
            fq,
            # args.input_files,
            reads,
            args.bwa_index,
            args.annotation_db,
            args.out_prefix,
            cpus_for_alignment=args.cpus_for_alignment,
            no_prefilter=args.no_prefilter,
            min_identity=args.min_identity,
            min_seqlen=args.min_seqlen,
            unmarked_orphans=input_type=="orphan",
        )

    fq.finalise(restrict_reports=("rpkm",))


if __name__ == "__main__":
    main()
