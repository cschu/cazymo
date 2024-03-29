# pylint: disable=C0103
""" module docstring """

import argparse
import json
import sys


def get_lines_from_chunks(_in, bufsize=400000000):
    tail = ""
    while True:
        chunk = "".join((tail, _in.read(bufsize)))
        if not chunk:
            break
        chunk = chunk.split("\n")
        *chunk, tail = chunk
        for line in chunk:
            yield line
    if tail:
        yield tail


def main():

    ap = argparse.ArgumentParser()
    ap.add_argument("out_prefix")
    ap.add_argument("--all", action="store_true")
    args = ap.parse_args()

    nreads, nalign, nlines = 0, 0, 0
    lastread = None

    nlines = -1
    for nlines, line in enumerate(get_lines_from_chunks(sys.stdin)):
        if line[0] != "@":
            nalign += 1
            rname = line[:line.find("\t")]
            if lastread is None or lastread != rname:
                nreads += 1
                lastread = rname

        print(line, file=sys.stdout)

    counts = {
        "n_lines": nlines, "n_align": nalign, "n_reads": nreads
    }
    # print(json.dumps(counts), file=sys.stderr)

    suffix = ".all" if args.all else ""
    with open(f"{args.out_prefix}{suffix}.readcount.json", "wt", encoding="UTF-8") as json_out:
        print(json.dumps(counts), file=json_out)


if __name__ == "__main__":
    main()
