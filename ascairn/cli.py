import click
from ascairn.commands import type

@click.group()
def main():
    """ASCairn: Alpha Satellite Centromere Analysis."""
    pass

main.add_command(type.type_command, name="type")

