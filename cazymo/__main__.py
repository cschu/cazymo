# pylint: disable=C0103,C0301

""" module docstring """

import contextlib
import json
import logging
import os
import pathlib
import subprocess
import sys

# pylint: disable=W0611
from gq.profilers import RegionQuantifier
from gq.db import get_writable_database
from gq.bin.build_bed_database import gather_category_and_feature_data, process_annotations
from gq.db.models import db

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

    _, db_session = get_writable_database()

    with db_session:
        code_map, nseqs = gather_category_and_feature_data(args.annotation_db, db_session=db_session)
        process_annotations(args.annotation_db, db_session, code_map, nseqs)

        print("QUERY", db_session.query(db.Feature).filter(db.Feature.id == 1).join(db.Category, db.Feature.category_id == db.Category.id).one_or_none().name)
        logging.info("Finished loading database.")

        fq = RegionQuantifier(
            db=db_session,
            out_prefix=args.out_prefix,
            ambig_mode=args.ambig_mode,
            strand_specific=args.strand_specific,
            calc_coverage=False,
            paired_end_count=args.paired_end_count,
            unmarked_orphans=args.unmarked_orphans,
            reference_type="domain",
        )

        # cmd = [f"bwa mem -a -t {args.cpus} -K 10000000 {args.bwa_index} {' '.join(args.input_files)}"]
        # if args.no_prefilter:
        #     cmd += ["samtools view -F 4 -buSh -"]
        # else:
        #     cmd += ["samtools view -F 4 -Sh -"]
        #     cmd += ["count_reads"]
        #     cmd += ["samtools view -buSh -"]
        #     cmd += ["bedtools intersect -u -ubam -a stdin -b {args.annotation_db}"]

        # with subprocess.Popen(
        #     "|".join(cmd), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,    
        # ) as align_flow:
        #     fq.process_bamfile(

        #     )

        # pylint: disable=E1124
        with subprocess.Popen(
            ("bwa", "mem", "-a", "-t", f"{args.cpus}", "-K", "10000000", args.bwa_index, *args.input_files),
            stdout=subprocess.PIPE
        ) as bwa_proc:
            if args.no_prefilter:
                samtools_filter_proc = subprocess.Popen(
                    ("samtools", "view", "-F 4", "-buSh", "-",), stdin=bwa_proc.stdout, stdout=subprocess.PIPE
                )
                read_count_proc = contextlib.nullcontext()
                samtools_convert_proc = contextlib.nullcontext()
                bedtools_proc = contextlib.nullcontext()

                align_stream = samtools_filter_proc.stdout
            else:
                logging.info("Prefiltering activated.")
                samtools_filter_proc = subprocess.Popen(
                    ("samtools", "view", "-F 4", "-Sh", "-",), stdin=bwa_proc.stdout, stdout=subprocess.PIPE
                )
                read_count_proc = subprocess.Popen(
                    ("read_count", args.out_prefix), stdin=samtools_filter_proc.stdout, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                )
                samtools_convert_proc = subprocess.Popen(
                    ("samtools", "view", "-buSh", "-",), stdin=read_count_proc.stdout, stdout=subprocess.PIPE
                )
                bedtools_proc = subprocess.Popen(
                    ("bedtools", "intersect", "-u", "-ubam", "-a", "stdin", "-b", f"{args.annotation_db}",),
                    stdin=samtools_convert_proc.stdout, stdout=subprocess.PIPE
                )
                align_stream = bedtools_proc.stdout

            # with samtools_filter_proc, read_count_proc, samtools_convert_proc, bedtools_proc:
            with bedtools_proc:
                # readcounts = None
                # if not args.no_prefilter:
                #     try:
                #         # readcounts = json.loads(open(args.out_prefix + ".readcount.json", "rt")).get("n_reads", 0)
                #         readcounts = json.loads(read_count_proc.stderr).get("n_reads")
                #     except:
                #         logger.warn("Could not access pre-filter readcounts. Using post-filter readcounts.")
                        
                fq.process_bamfile(
                    align_stream,
                    aln_format="bam",
                    min_identity=args.min_identity, min_seqlen=args.min_seqlen,
                    external_readcounts=None if args.no_prefilter else read_count_proc.stderr,
                )


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
