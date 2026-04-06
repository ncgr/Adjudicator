import os

import intervaltree


class WAOIntersect:
    """Reproduce the output of ``bedtools intersect -wao`` for two GFF3 files.

    GFF3 coordinates are 1-based fully-closed.  Coordinate conversion to
    0-based half-open intervals is performed internally for IntervalTree
    arithmetic only.

    Output column layout (19 columns, native GFF3 order throughout):
        [0-8]  A record: seqname, source, feature, start, end,
                         score, strand, frame, attributes
        [9-17] B record: same nine columns  (``"."`` fields when no overlap)
        [18]   overlap base-pairs  (``"0"`` when no overlap)
    """

    _NULL_B = [".", ".", ".", "-1", "-1", ".", ".", ".", "."]

    def __init__(self, gff3_a: str, gff3_b: str):
        """Initialise, parse both GFF3 files, and build the B interval trees.

        Args:
            gff3_a: Path to the query GFF3 file (the ``-a`` operand).
            gff3_b: Path to the database GFF3 file (the ``-b`` operand).
        """
        self.gff3_a = os.path.abspath(gff3_a)
        self.gff3_b = os.path.abspath(gff3_b)
        self.trees_b: dict[str, intervaltree.IntervalTree] = {}
        self.records_a: list[tuple[str, int, int, list[str]]] = []
        self._load_b()
        self._load_a()

    def _parse_gff3(self, path: str) -> list[tuple[str, int, int, list[str]]]:
        """Parse a GFF3 file into a list of interval records.

        Comment lines (``#``) and blank lines are skipped.  Lines with fewer
        than nine tab-delimited fields are skipped silently.

        The original nine fields are preserved unchanged for output.
        Coordinates are also stored in 0-based half-open form for internal
        IntervalTree arithmetic only.

        Args:
            path: Absolute path to the GFF3 file.

        Returns:
            A list of ``(seqname, start_0, end, fields)`` tuples where
            ``start_0`` is 0-based, ``end`` is the GFF3 end (used as the
            half-open upper bound), and ``fields``
        """
        records = []
        with open(path) as fh:
            for line in fh:
                if line.startswith("#") or not line.strip():
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 9:
                    continue
                seqname = fields[0]
                start_0 = int(fields[3]) - 1
                end = int(fields[4])
                records.append((seqname, start_0, end, fields))
        return records

    def _load_b(self):
        """Parse ``gff3_b`` and insert every record into a per-sequence IntervalTree.

        Each interval's data payload is the original nine GFF3 fields so that
        they can be appended verbatim to an A row in ``intersect``.
        """
        for seqname, start_0, end, fields in self._parse_gff3(self.gff3_b):
            if seqname not in self.trees_b:
                self.trees_b[seqname] = intervaltree.IntervalTree()
            self.trees_b[seqname][start_0:end] = fields

    def _load_a(self):
        """Parse ``gff3_a`` and store the records for iteration in ``intersect``."""
        self.records_a = self._parse_gff3(self.gff3_a)

    def intersect(self) -> list[list[str]]:
        """Run the ``-wao`` intersection and return all output rows.

        For every A record all overlapping B intervals are found via the
        IntervalTree.

        ``overlap = min(end_a, iv.end) - max(start_a, iv.begin)``.

        A records with no matching B sequence name, or with no overlapping B
        interval, are emitted once with ``_NULL_B`` fields and an overlap of
        ``"0"``.
        """
        rows = []
        self.non_overlapping_b: set[tuple] = set()

        for seqname, start_0, end, fields in self._parse_gff3(self.gff3_b):
            self.non_overlapping_b.add(tuple(fields))

        for seqname, start_a, end_a, fields_a in self.records_a:
            tree = self.trees_b.get(seqname)
            if tree is None:
                rows.append(fields_a + self._NULL_B + ["0"])
                continue

            overlaps = sorted(
                tree.overlap(start_a, end_a),
                key=lambda iv: (iv.begin, iv.end),
            )
            if not overlaps:
                rows.append(fields_a + self._NULL_B + ["0"])
                continue

            for iv in overlaps:
                overlap_bp = min(end_a, iv.end) - max(start_a, iv.begin)
                self.non_overlapping_b.discard(tuple(iv.data))  # mark as matched
                rows.append(fields_a + list(iv.data) + [str(overlap_bp)])

        return rows

    def write(self, output_prefix: str):
        """Write the intersections to two gff3(ish) files to mimic wao
        and report non-overlapping features in b.
        """
        wao_out = f"{output_prefix}.wao.gff3"
        unique_b_out = f"{output_prefix}.unique_b.gff3"
        with open(wao_out, "w") as fh:
            for row in self.intersect():
                fh.write("\t".join(row) + "\n")

        with open(unique_b_out, "w") as fh:
            for row in self.non_overlapping_b:
                fh.write("\t".join(row) + "\n")


if __name__ == "__main__":
    test_a = "./test_a.gff3"
    test_b = "./test_b.gff3"
    wao = WAOIntersect(test_a, test_b)
    wao.write("./test.wao.gff3")
