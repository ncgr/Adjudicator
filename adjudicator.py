#!/usr/bin/env python3
"""Adjudicator: A CLI tool for processing GFF3 and GFA files."""

import csv
import os

import click

from tools.waointersect import WAOIntersect
from tools.adjudicate_model import AdjudicateModel
from tools.gff3_in_memory import GFF3InMemory


def _validate_sample_row(row: list[str], line_num: int) -> tuple[str, str, str]:
    if len(row) != 3:
        raise click.BadParameter(
            f"Line {line_num}: expected 3 columns, got {len(row)}.",
            param_hint="--input-tsv",
        )

    label, gff3_path, gfa_path = row[0].strip(), row[1].strip(), row[2].strip()

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


def _parse_tsv(input_tsv: str, strict: bool) -> list[tuple[str, str, str]]:
    rows = []
    with open(input_tsv, newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        for line_num, row in enumerate(reader, start=1):
            if not row or row[0].startswith("#"):
                continue
            label, gff3_path, gfa_path = row[0].strip(), row[1].strip(), row[2].strip()
            if strict:  # placeholder, method is a bit silly
                label, gff3_path, gfa_path = _validate_sample_row(row, line_num)
            rows.append([label, gff3_path, gfa_path])

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


@cli.command("repeat-filter")
@click.option(
    "--input-tsv",
    "-i",
    required=True,
    type=click.Path(exists=True, dir_okay=False, readable=True, path_type=str),
    help=(
        """Tab-separated file with three columns:

          COL 1  Label / sample name (string)\n
          COL 2  Path to a .gff3 annotation file\n
          COL 3  Path to a .gfa LIS gene family assignment file.

        The top-down order is the order of precedence when models have highly
        similar scores."""
    ),
)
@click.option(
    "--annotation",
    "-a",
    required=True,
    type=click.Path(exists=True, dir_okay=False, readable=True, path_type=str),
    help=("""GFF3 file of annotated regions to exclude."""),
)
@click.option(
    "--max-coverage",
    "-m",
    default=0.4,
    show_default=True,
    type=float,
    help=(
        """Maximum overlap coverage between features and annotations.

           The default feature is 'exon'. Adjust this with '--feature'
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
def repeat_filter(
    input_tsv: str,
    annotation: str,
    max_coverage: float,
    output_dir: str,
    strict: bool,
    verbose: bool
):
    """
    Filter the gff3 structural annotations from the input TSV file against
    the provided 'annotation' file.

    Usage:

        adjudicator repeat-filter <--input-tsv samples.tsv> <--annotation repeats.gff3> [options] [--help]
    """
    rows = _parse_tsv(input_tsv, strict)
    click.echo(
        f"Parsed {len(rows)} entr{'y' if len(rows) == 1 else 'ies'} from '{input_tsv}'."
    )
    os.makedirs(output_dir, exist_ok=True)

    count = 0
    compare = list()
    for label, gff3_path, gfa_path in rows:
        if verbose:
            click.echo(f"Processing '{label}':")
            click.echo(f"GFF3 : {gff3_path}")
            click.echo(f"GFA  : {gfa_path}")
        compare.append([label, gff3_path, gfa_path])
        label_a = label
        label_b = "repeat_filter"
        gff3_a = gff3_path
        gff3_b = annotation
        gfa_a = gfa_path
        name = f"{label_a}_{label_b}"
        compare_out = f"{output_dir}/{name}"
        overlaps = f"{compare_out}.wao.gff3"
        intersection = WAOIntersect(gff3_a, gff3_b, True, 0.4)
        intersection.write(compare_out)
        adjudicator = AdjudicateModel(
            overlaps, gfa_a, None, 0.4, False, compare_out
        )
        blacklist = adjudicator.repeat_filter()
        filtered_blacklist = []
        for gid in blacklist:
            filtered_blacklist.append(".".join(gid.split(".")[:-1]))
        filtered_gff3 = GFF3InMemory(gff3_a, filtered_blacklist)
        filtered_records_out = f"{compare_out}.final.wsubfeatures.gff3"
        with open(filtered_records_out, "w") as fh:
            for gid in filtered_gff3.gene_ids():
                filtered_gff3.write_gene(gid, fh)
    #        print(all_gffs.get_gene(gid))
    click.echo(f"Done. Results written to '{output_dir}'.")



@cli.command("collapse")
@click.option(
    "--input-tsv",
    "-i",
    required=True,
    type=click.Path(exists=True, dir_okay=False, readable=True, path_type=str),
    help=(
        """Tab-separated file with three columns:

          COL 1  Label / sample name (string)\n
          COL 2  Path to a .gff3 annotation file\n
          COL 3  Path to a .gfa LIS gene family assignment file.

        The top-down order is the order of precedence when models have highly
        similar scores."""
    ),
)
@click.option(
    "--min-overlap",
    "-m",
    default=0.00001,
    show_default=True,
    type=float,
    help="Minimum overlap coverage between A and B wrt to A.",
)
@click.option(
    "--no-orphans",
    "-n",
    is_flag=True,
    default=False,
    show_default=True,
    help="Do not include genes with no gene family assignment.",
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
def collapse(
    input_tsv: str,
    min_overlap: float,
    no_orphans: bool,
    output_dir: str,
    strict: bool,
    verbose: bool,
):
    """
    Collapse overlapping structural annotations using gene family assignments from the Legume Information System:

    https://data.legumeinfo.org/LEGUMES/Fabaceae/genefamilies/legume.fam3.VLMQ/

    Usage:

        adjudicator collapse <--input-tsv samples.tsv> [options] [--help]
    """
    rows = _parse_tsv(input_tsv, strict)
    click.echo(
        f"Parsed {len(rows)} entr{'y' if len(rows) == 1 else 'ies'} from '{input_tsv}'."
    )

    _check_files_exist(rows, strict=strict)
    os.makedirs(output_dir, exist_ok=True)

    count = 0
    compare = list()
    compare_out = ""
    last_final = ""
    last_comparison = ""
    last_gfa = ""
    name = ""
    for label, gff3_path, gfa_path in rows:
        if verbose:
            click.echo(f"Processing '{label}':")
            click.echo(f"GFF3 : {gff3_path}")
            click.echo(f"GFA  : {gfa_path}")
        compare.append([label, gff3_path, gfa_path])
        if count < 2:
            if count == 1:
                label_a = compare[0][0]
                label_b = compare[1][0]
                gff3_a = compare[0][1]
                gff3_b = compare[1][1]
                gfa_a = compare[0][2]
                gfa_b = compare[1][2]
                name = f"{label_a}_{label_b}"
                compare_out = f"{output_dir}/{name}"
                intersection = WAOIntersect(gff3_a, gff3_b, False, None)
                intersection.write(compare_out)
                overlaps = f"{compare_out}.wao.gff3"
                unique_b = f"{compare_out}.unique_b.gff3"
                adjudicator = AdjudicateModel(
                    overlaps, gfa_a, gfa_b, min_overlap, no_orphans, compare_out
                )
                adjudicator.choose_model()
        else:
            name = f"{last_comparison}_{label}"
            compare_out = f"{output_dir}/{name}"
            intersection = WAOIntersect(last_final, gff3_path, False, None)
            intersection.write(compare_out)
            overlaps = f"{compare_out}.wao.gff3"
            unique_b = f"{compare_out}.unique_b.gff3"
            adjudicator = AdjudicateModel(
                overlaps, last_gfa, gfa_path, min_overlap, no_orphans, compare_out
            )
            adjudicator.choose_model()

        last_comparison = name
        last_final = f"{compare_out}.final.gff3"
        last_gfa = f"{compare_out}.gfa"
        count += 1

    input_gffs = [sample[1] for sample in compare]
    all_gffs = GFF3InMemory(input_gffs, [])
    final_gff3 = GFF3InMemory(last_final, [])
    all_records_out = f"{compare_out}.final.wsubfeatures.gff3"
    with open(all_records_out, "w") as fh:
        for gid in final_gff3.gene_ids():
            all_gffs.write_gene(gid, fh)
    #        print(all_gffs.get_gene(gid))
    click.echo(f"Done. Results written to '{output_dir}'.")


if __name__ == "__main__":
    cli()
