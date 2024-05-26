from pathlib import Path

import click

from scitrack import CachingLogger

from cogent3 import make_unaligned_seqs
from cogent3.core.alignment import SequenceCollection

__author__ = "Zongjing Han"
__copyright__ = "Copyright 2016-2021, Zongjing Han"
__credits__ = ["Zongjing Han"]
__license__ = "BSD"
__version__ = "2020.6.5"  # A DATE BASED VERSION
__maintainer__ = "Zongjing Han"
__email__ = "u7627363@anu.edu.au"
__status__ = "alpha"


LOGGER = CachingLogger()


@click.group()
@click.version_option(__version__)  # add version option
def main():
    """docstring explains what your cl app does"""
    pass

def unique_kmers(seqs: SequenceCollection, k: int = 2) -> dict[str, set[str]]:
    """returns a {seqname: set(str, ...)}"""
    # Dictionary to store k-mers and the sequences they appear in
    kmer_dict = {}
    
    # Iterate over sequences and extract k-mers
    for s in seqs.seqs:
        name = s.name
        sequence = str(s)  # Access the sequence data correctly
        for i in range(len(sequence) - k + 1):
            kmer = sequence[i:i + k]
            if kmer not in kmer_dict:
                kmer_dict[kmer] = set()
            kmer_dict[kmer].add(name)

    unique_kmers_dict = {}
    for kmer, names in kmer_dict.items():
        if len(names) == 1:  
            name = next(iter(names))
            if name not in unique_kmers_dict:
                unique_kmers_dict[name] = set()
            unique_kmers_dict[name].add(kmer)
    
    return unique_kmers_dict

@main.command(no_args_is_help=True)
@click.option(
    "-i",
    "--infile",
    required=True,
    type=click.Path(exists=True),
    help="Path to the input file containing sequences"
)
@click.option(
    "-k",
    "--kmer_length",
    required=True,
    type=int,
    help="Length of the k-mers"
)
@click.option(
    "-o",
    "--outfile",
    required=True,
    type=click.Path(),
    help="Path to the output file"
)
def find_unique_kmers(infile, kmer_length, outfile):
    """Identify unique k-mers from sequences"""
    with open(infile, 'r') as f:
        seqs_dict = {}
        seq_id = None
        sequence = []
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if seq_id is not None:
                    seqs_dict[seq_id] = "".join(sequence)
                seq_id = line[1:]  # Get rid of '>'
                sequence = []
            else:
                sequence.append(line)
        if seq_id is not None:  # For the last sequence
            seqs_dict[seq_id] = "".join(sequence)

    seqs = SequenceCollection(seqs_dict, moltype="dna")
    unique_kmers_dict = unique_kmers(seqs, kmer_length)

    with open(outfile, 'w') as f:
        for seq_id, kmers in unique_kmers_dict.items():
            f.write(f"{seq_id}\t{','.join(kmers)}\n")

    click.echo(f"Unique k-mers written to {outfile}")



# note that I often define reusable options in a central place
_verbose = click.option(
    "-v",
    "--verbose",
    count=True,
    help="is an integer indicating number of cl occurrences",
)

# you can define custom parsers / validators
def _parse_csv_arg(*args) -> list:
    return args[-1].split(",")


_names = click.option(
    "--names",
    callback=_parse_csv_arg,
    help="converts comma separated values",
)

_outpath = click.option(
    "-o", "--outpath", type=Path, help="the input string will be cast to Path instance"
)

# the no_args_is_help=True means help is displayed if a
# user doesn't provide any arguments to a subcommand.
# Should be a click default I think!
@main.command(no_args_is_help=True)
@click.option(
    "-i",
    "--infile",
    required=True,
    type=click.Path(exists=True),
    help="fails if provided value is non-existent path",
)
@_outpath
@click.option(
    "--achoice",
    type=click.Choice(["choice1", "choice2"]),
    default="choice1",
    help="make a choice",
)
@_names
@click.option(
    "-O",
    "--overwrite",
    is_flag=True,
    help="overwrite results",
)
@click.option(
    "--ensembl_account",
    envvar="ENSEMBL_ACCOUNT",
    help="shell variable with MySQL account "
    "details, e.g. export "
    "ENSEMBL_ACCOUNT='myhost.com jill jills_pass'",
)
@_verbose
def demo_log(infile, outpath, achoice, names, overwrite, verbose, ensembl_account):
    """what demo_log subcommand does"""
    # capture the local variables, at this point just provided arguments
    LOGGER.log_args()
    LOGGER.log_versions("numpy")
    LOGGER.input_file(infile)

    LOGGER.log_file_path = outpath / "some_path.log"


def _parse_csv_arg(*args):
    return args[-1].split(",")


@main.command(no_args_is_help=True)
@click.argument("message", required=True, type=str)
@click.option("-t", "--test", is_flag=True, help="test run")
@_verbose
def demo_echo(message, test, verbose):
    """what demo_echo subcommand does"""
    for _ in range(verbose):
        click.secho(message, fg="blue")


if __name__ == "__main__":
    main()
