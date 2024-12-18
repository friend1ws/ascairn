import click
from ascairn.commands import parse_marker
from ascairn.commands import kmer_count
from ascairn.commands import cen_type
from ascairn.commands import check_depth

@click.group()
def main():
    """ASCairn: Alpha Satellite Centromere Analysis."""
    pass

main.add_command(kmer_count.kmer_count_command, name = "kmer_count")
main.add_command(parse_marker.parse_marker_command, name="parse_marker")
main.add_command(cen_type.cen_type_command, name="cen_type")
main.add_command(check_depth.check_depth_command, name="check_depth")


