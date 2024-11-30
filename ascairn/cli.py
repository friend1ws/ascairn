import click
from ascairn.commands import type
from ascairn.commands import parse_marker
from ascairn.commands import kmer_count

@click.group()
def main():
    """ASCairn: Alpha Satellite Centromere Analysis."""
    pass

main.add_command(type.type_command, name="type")
main.add_command(kmer_count.kmer_count_command, name = "kmer_count")
main.add_command(parse_marker.parse_marker_command, name="parse_marker")

