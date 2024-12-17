import click
from ascairn.commands import quick_type_all
from ascairn.commands import parse_marker
from ascairn.commands import kmer_count
from ascairn.commands import type
from ascairn.commands import check_depth

@click.group()
def main():
    """ASCairn: Alpha Satellite Centromere Analysis."""
    pass

# main.add_command(quick_type_all.quick_type_all_command, name="quick_type_all")
main.add_command(kmer_count.kmer_count_command, name = "kmer_count")
main.add_command(parse_marker.parse_marker_command, name="parse_marker")
main.add_command(type.type_command, name="type")
main.add_command(check_depth.check_depth_command, name="check_depth")


