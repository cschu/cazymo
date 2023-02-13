# pylint: disable=C0103,C0301

""" module docstring """

import contextlib
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


def main():

    args = handle_args(sys.argv[1:])

    logger.info("Version: %s", __version__)
    logger.info("Command: %s %s", os.path.basename(sys.argv[0]), " ".join(sys.argv[1:]))

    print(args)
    if args.input_files != "-" and not all(os.path.exists(f) for f in args.input_files):
        input_files_str = "\n".join(args.input_files)
        raise ValueError(f"There is an issue with your input files. Please check.\n{input_files_str}")
    if not os.path.exists(args.annotation_db):
        raise ValueError(f"{args.annotation_db} is not a valid annotation database")
    if not check_bwa_index(args.bwa_index):
        raise ValueError(f"{args.bwa_index} is not a valid bwa index.")

    if os.path.dirname(args.out_prefix):
        pathlib.Path(os.path.dirname(args.out_prefix)).mkdir(
            exist_ok=True, parents=True
        )

    db_importer = DomainBedDatabaseImporter(logger, args.annotation_db)
    # db_importer.gather_category_and_feature_data(args.annotation_db)
    # db_importer.process_annotations(args.annotation_db)
    logger.info("Finished loading database.")

    fq = RegionQuantifier(
        db=db_importer,
        out_prefix=args.out_prefix,
        ambig_mode="1overN",
        strand_specific=args.strand_specific,
        unmarked_orphans=args.unmarked_orphans,
        reference_type="domain",
    )

    try:
        # pylint: disable=E1124,R1732
        logging.info("Initiated read alignment process.")
        with subprocess.Popen(
            ("bwa", "mem", "-a", "-t", f"{args.cpus_for_alignment}", "-K", "10000000", args.bwa_index, *args.input_files),
            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        ) as bwa_proc:
            if args.no_prefilter:
                samtools_filter_proc = subprocess.Popen(
                    ("samtools", "view", "-F 4", "-buSh", "-",), stdin=bwa_proc.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                )
                read_count_proc = contextlib.nullcontext()
                samtools_convert_proc = contextlib.nullcontext()
                bedtools_proc = contextlib.nullcontext()

                align_stream = samtools_filter_proc.stdout
            else:
                logging.info("Prefiltering activated.")
                samtools_filter_proc = subprocess.Popen(
                    ("samtools", "view", "-F 4", "-Sh", "-",), stdin=bwa_proc.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                )
                read_count_proc = subprocess.Popen(
                    ("read_count", args.out_prefix), stdin=samtools_filter_proc.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                )
                samtools_convert_proc = subprocess.Popen(
                    ("samtools", "view", "-buSh", "-",), stdin=read_count_proc.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                )
                bedtools_proc = subprocess.Popen(
                    ("bedtools", "intersect", "-u", "-ubam", "-a", "stdin", "-b", f"{args.annotation_db}",),
                    stdin=samtools_convert_proc.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                )
                align_stream = bedtools_proc.stdout

            # with samtools_filter_proc, read_count_proc, samtools_convert_proc, bedtools_proc:
            with bedtools_proc:
                fq.process_bamfile(
                    align_stream,
                    aln_format="bam",
                    min_identity=args.min_identity, min_seqlen=args.min_seqlen,
                    external_readcounts=None if args.no_prefilter else (args.out_prefix + ".readcount.json")
                )

    except Exception as err:
        logger.error("Caught exception:")
        logger.error("%s", err)
        raise Exception from err

        # with subprocess.Popen(
        #     ("samtools", "view", "-F", "4", "-buSh", "-"),
        #     stdin=bwa_proc.stdout, stdout=subprocess.PIPE
        # ) as samtools_proc:
        #     if args.no_prefilter:
        #         bedtools_proc = contextlib.nullcontext()
        #         fq_input = samtools_proc.stdout
        #     else:
        #         logging.info("Prefiltering activated.")
        #         bedtools_proc = subprocess.Popen(
        #             ("bedtools", "intersect", "-u", "-ubam", "-a", "stdin", "-b", f"{args.annotation_db}"),
        #             stdin=samtools_proc.stdout, stdout=subprocess.PIPE
        #         )
        #         # bedtools intersect -u -ubam -a ${bam} -b ${db_bedfile} > filtered_bam/${sample}.bam

        #         fq_input = bedtools_proc.stdout

        #     with bedtools_proc:
        #         fq.process_bamfile(
        #             fq_input,
        #             aln_format="bam",
        #             min_identity=args.min_identity, min_seqlen=args.min_seqlen
        #         )


if __name__ == "__main__":
    main()
