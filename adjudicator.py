#!/usr/bin/env python3
"""Adjudicator: A CLI tool for processing GFF3 and GFA files."""

import csv
import os

import click

from tools.waointersect import WAOIntersect


def _validate_sample_row(row: list[str], line_num: int) -> tuple[str, str, str]:
    if len(row) != 3:
        raise click.BadParameter(
            f"Line {line_num}: expected 3 columns, got {len(row)}.",
            param_hint="--input-tsv",
        )

    label, gff3_path, gfa_path = row[0].strip(), row[1].strip, row[2].strip

    if not label:
        raise click.BadParameter(
            f"Line {line_num}: column 1 (label) must not be empty.",
            param_hint="--input-tsv",
        )

    if not gff3_path.endswith(".gff3"):
        raise click.BadParameter(
            f"Line {line_num}: column 2 must end in '.gff3', got '{gff3_path}'.",
            param_hint="--input-tsv",
        )

    if not gfa_path.endswith(".gfa"):
        raise click.BadParameter(
            f"Line {line_num}: column 3 must end in '.gfa', got '{gfa_path}'.",
            param_hint="--input-tsv",
        )

    return label, gff3_path, gfa_path


def _check_files_exist(rows: list[tuple[str, str, str]], strict: bool) -> None:
    missing = [
        f"  [{label}] {path}"
        for label, gff3_path, gfa_path in rows
        for path in (gff3_path, gfa_path)
        if not os.path.exists(path)
    ]

    if missing:
        msg = "The following files were not found:\n" + "\n".join(missing)
        if strict:
            raise click.ClickException(msg)
        click.echo(f"Warning: {msg}", err=True)


def _parse_tsv(input_tsv: str) -> list[tuple[str, str, str]]:
    rows = []
    with open(input_tsv, newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        for line_num, raw_row in enumerate(reader, start=1):
            if not raw_row or raw_row[0].startswith("#"):
                continue
            rows.append(_validate_sample_row(raw_row, line_num))

    if not rows:
        raise click.ClickException(f"No data rows found in '{input_tsv}'.")

    return rows


@click.group()
@click.version_option(version="0.1.0", prog_name="adjudicator")
def cli():
    """
    The Adjudicator collapses multiple structural annotations using
    best hit domain hmm scores for gene family assignments as a basis
    to compare overlapping models and their respective exons.

    It will also attempt to insert unique annotations from each set
    and remove (or assign) calls against repeat annotations if available.
    """


@cli.command("collapse")
@click.option(
    "--input-tsv",
    "-i",
    required=True,
    type=click.Path(exists=True, dir_okay=False, readable=True, path_type=str),
    help=(
        """
       Tab-separated file with three columns:

          COL 1  Label / sample name (string)
          COL 2  Path to a .gff3 annotation file
          COL 3  Path to a .gfa graph file

        The top-down order is the order of precedence when models have highly
        similar scores.
    """
    ),
)
@click.option(
    "--output-dir",
    "-o",
    default=".",
    show_default=True,
    type=click.Path(file_okay=False, writable=True, path_type=str),
    help="Directory where output files will be written.",
)
@click.option(
    "--strict/--no-strict",
    default=False,
    show_default=True,
    help="Exit with an error if any referenced input file does not exist on disk.",
)
@click.option(
    "--verbose",
    "-v",
    is_flag=True,
    default=False,
    help="Enable verbose logging.",
)
def collapse(input_tsv: str, output_dir: str, strict: bool, verbose: bool):
    """Collapse structural annotations in top-down order.

    \b
    The TSV must contain exactly three tab-separated columns (no header):

      COLUMN 1  A unique label or sample name  (string)
      COLUMN 2  Path to a GFF3 annotation file (.gff3)
      COLUMN 3  Path to a GFA graph file       (.gfa)

    \b
    Example TSV:
      sample_A  /data/sample_A.gff3  /data/sample_A.gfa
      sample_B  /data/sample_B.gff3  /data/sample_B.gfa

    \b
    Usage:
      adjudicator collapse --input-tsv samples.tsv [options] [--help]
    """
    rows = _parse_tsv(input_tsv)
    click.echo(
        f"Parsed {len(rows)} entr{'y' if len(rows) == 1 else 'ies'} from '{input_tsv}'."
    )

    _check_files_exist(rows, strict=strict)
    os.makedirs(output_dir, exist_ok=True)

    for label, gff3_path, gfa_path in rows:
        if verbose:
            click.echo(f"Processing '{label}':")
            click.echo(f"GFF3 : {gff3_path}")
            click.echo(f"GFA  : {gfa_path}")

    click.echo(f"Done. Results written to '{output_dir}'.")


if __name__ == "__main__":
    cli()
