"""
GFF3InMemory: Parses a GFF3 file into a hierarchical in-memory data structure.

Hierarchy: gene -> transcript -> exon -> CDS
Only gene, transcript, exon, and CDS features are parsed.
"""

from collections import defaultdict
from typing import Optional, Union


class GFF3InMemory:
    """
    Reads a GFF3 file and builds an in-memory hierarchical data structure.

    Structure of gff3 dict:
    {
        "<gene_id>": {
            "feature": "gene",
            "seqid": str,
            "source": str,
            "start": int,
            "end": int,
            "score": str,
            "strand": str,
            "phase": str,
            "attributes": dict,
            "children": {
                "<transcript_id>": {
                    "feature": "transcript",
                    ...same fields as gene...,
                    "children": {
                        "<exon_id>": {
                            "feature": "exon",
                            ...same fields...,
                            "children": {}
                        },
                        "<cds_id>": {
                            "feature": "CDS",
                            ...same fields...,
                            "children": {}
                        }
                    }
                }
            }
        }
    }
    """

    FEATURE_ORDER = ["gene", "mRNA", "exon", "CDS"]
    FEATURE_RANK = {f: i for i, f in enumerate(FEATURE_ORDER)}

    def __init__(self, input_gff3: Union[str, list]):
        """
        Parameters
        ----------
        input_gff3 : str or list of str
            Path to a single GFF3 file, or a list of paths to multiple GFF3
            files.  All files are merged into one shared hierarchy; if the
            same feature ID appears in more than one file the last file's
            record wins.
        """
        if isinstance(input_gff3, str):
            self.input_gff3 = [input_gff3]
        else:
            self.input_gff3 = list(input_gff3)

        self.gff3: dict = {}
        self._id_index: dict = {}
        self._parse()

    def _parse_attributes(self, attr_str: str) -> dict:
        """Parse the GFF3 column-9 attribute string into a dict."""
        attrs = {}
        for field in attr_str.split(";"):
            field = field.strip()
            if not field or "=" not in field:
                continue
            key, _, value = field.partition("=")
            attrs[key.strip()] = value.strip()
        return attrs

    def _make_record(self, parts: list) -> dict:
        """Build a record dict from the 9 GFF3 columns."""
        attributes = self._parse_attributes(parts[8]) if len(parts) > 8 else {}
        return {
            "feature": parts[2],
            "seqid": parts[0],
            "source": parts[1],
            "start": int(parts[3]),
            "end": int(parts[4]),
            "score": parts[5],
            "strand": parts[6],
            "phase": parts[7],
            "attributes": attributes,
            "children": {},
        }

    def _parse(self):
        """
        Read all GFF3 files, keep only the four target features, sort them
        by feature rank (gene first, CDS last), then wire the hierarchy.
        Records from all files are pooled before sorting so the hierarchy
        is built correctly regardless of file order.
        """
        target_features = set(self.FEATURE_ORDER)
        raw_records = []

        for gff3_file in self.input_gff3:
            with open(gff3_file, "r") as fh:
                for line in fh:
                    line = line.rstrip()

                    if not line or line.startswith("#"):
                        continue

                    parts = line.split("\t")
                    if len(parts) < 9:
                        continue

                    feature_type = parts[2]
                    if feature_type not in target_features:
                        continue

                    record = self._make_record(parts)
                    raw_records.append(record)

        raw_records.sort(key=lambda r: self.FEATURE_RANK[r["feature"]])

        for record in raw_records:
            attrs = record["attributes"]
            feature_id = attrs.get("ID", "")
            parent_id = attrs.get("Parent", "")

            if record["feature"] == "gene":
                self.gff3[feature_id] = record
                self._id_index[feature_id] = record

            else:
                if parent_id and parent_id in self._id_index:
                    parent_record = self._id_index[parent_id]
                    parent_record["children"][feature_id] = record
                    self._id_index[feature_id] = record
                else:
                    self._id_index[feature_id] = record

    def get_gene(self, gene_id: str) -> Optional[dict]:
        """Return the gene record for *gene_id*, or None if not found."""
        return self.gff3.get(gene_id)

    def get_transcript(self, transcript_id: str) -> Optional[dict]:
        """Return any transcript record by its ID (searches flat index)."""
        record = self._id_index.get(transcript_id)
        if record and record["feature"] == "transcript":
            return record
        return None

    def gene_ids(self) -> list:
        """Return a list of all gene IDs."""
        return list(self.gff3.keys())

    def transcripts_for_gene(self, gene_id: str) -> dict:
        """Return the children dict (transcripts) of a gene."""
        gene = self.gff3.get(gene_id, {})
        return gene.get("children", {})

    def exons_for_transcript(self, transcript_id: str) -> dict:
        """Return exon children for a transcript."""
        tx = self.get_transcript(transcript_id)
        if tx is None:
            return {}
        return {k: v for k, v in tx["children"].items() if v["feature"] == "exon"}

    def cds_for_transcript(self, transcript_id: str) -> dict:
        """Return CDS children for a transcript."""
        tx = self.get_transcript(transcript_id)
        if tx is None:
            return {}
        return {k: v for k, v in tx["children"].items() if v["feature"] == "CDS"}

    def _record_to_gff3_line(self, feature_id: str, record: dict) -> str:
        """Serialise a single record back to a GFF3-formatted line."""
        attrs = record["attributes"]
        # Rebuild attribute string; ensure ID comes first, then Parent if present
        attr_parts = []
        if "ID" in attrs:
            attr_parts.append("ID={}".format(attrs["ID"]))
        if "Parent" in attrs:
            attr_parts.append("Parent={}".format(attrs["Parent"]))
        for k, v in attrs.items():
            if k not in ("ID", "Parent"):
                attr_parts.append("{}={}".format(k, v))
        attr_str = ";".join(attr_parts)

        return "\t".join(
            [
                record["seqid"],
                record["source"],
                record["feature"],
                str(record["start"]),
                str(record["end"]),
                record["score"],
                record["strand"],
                record["phase"],
                attr_str,
            ]
        )

    def _flatten_record(self, record: dict) -> list:
        """
        Recursively walk a record and all its children in hierarchy order
        (gene -> transcript -> exon -> CDS) and return a flat list of
        GFF3-formatted lines.
        """
        feature_id = record["attributes"].get("ID", "")
        lines = [self._record_to_gff3_line(feature_id, record)]
        # Sort children by feature rank so output order is deterministic
        children = sorted(
            record["children"].values(),
            key=lambda r: self.FEATURE_RANK.get(r["feature"], 99),
        )
        for child in children:
            lines.extend(self._flatten_record(child))
        return lines

    def flatten_gene(self, gene_id: str) -> list:
        """
        Return a flat list of GFF3-formatted strings for *gene_id* and all
        of its descendant features (transcripts, exons, CDS).

        Returns an empty list if *gene_id* is not found.
        """
        gene = self.gff3.get(gene_id)
        if gene is None:
            return []
        return self._flatten_record(gene)

    def print_gene(self, gene_id: str) -> None:
        """Print the GFF3 lines for *gene_id* and all its children to stdout."""
        for line in self.flatten_gene(gene_id):
            print(line)

    def write_gene(self, gene_id: str, file_handle) -> None:
        """
        Write the GFF3 lines for *gene_id* and all its children to an open
        file handle.

        Example
        -------
        with open("out.gff3", "w") as fh:
            gff.write_gene("gene:ENSG00000001", fh)
        """
        for line in self.flatten_gene(gene_id):
            file_handle.write(line + "\n")

    def print_all(self) -> None:
        """Print every gene and all its children to stdout."""
        for gene_id in self.gff3:
            self.print_gene(gene_id)

    def __len__(self) -> int:
        """Number of top-level gene records."""
        return len(self.gff3)

    def __repr__(self) -> str:
        n = len(self.input_gff3)
        files = self.input_gff3[0] if n == 1 else f"{n} files"
        return f"GFF3InMemory(files={files!r}, genes={len(self)})"
