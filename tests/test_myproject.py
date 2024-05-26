import pytest
from click.testing import CliRunner
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))

from myproject.cli import demo_echo, demo_log, find_unique_kmers



DATADIR = Path(__file__).parent / "data"


@pytest.fixture(scope="session")
def tmp_dir(tmpdir_factory):
    return tmpdir_factory.mktemp("demo")


@pytest.fixture(scope="session")
def runner():
    """exportrc works correctly."""
    return CliRunner()


def test_demo_echo(tmp_dir, runner):
    args = "'hello!' -vv".split()
    r = runner.invoke(demo_echo, args, catch_exceptions=False)
    assert r.exit_code == 0, r.output


def test_demo_log(tmp_dir, runner):
    args = f"-i {DATADIR / 'demo.fasta'} --names human,chimp -o {tmp_dir}".split()
    r = runner.invoke(demo_log, args, catch_exceptions=False)
    assert r.exit_code == 0, r.output

def test_find_unique_kmers(tmp_dir, runner):
    input_file = DATADIR / 'demo.fasta'
    output_file = tmp_dir / 'output.txt'
    args = f"-i {input_file} -k 2 -o {output_file}".split()
    r = runner.invoke(find_unique_kmers, args, catch_exceptions=False)
    assert r.exit_code == 0, r.output
    with open(output_file, 'r') as f:
        lines = f.readlines()
    assert len(lines) > 0  
    expected_output = {
        "s1": {"AA"},
        "s2": {"GA", "TG"},
        "s3": {"AC", "CT"}
    }
    for line in lines:
        seq_id, kmers = line.strip().split("\t")
        kmer_set = set(kmers.split(","))
        assert kmer_set == expected_output[seq_id]