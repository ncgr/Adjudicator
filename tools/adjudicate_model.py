#!/usr/bin/env python3
"""Adjudicate overlapping gene model predictions using gene family scores."""

import sys
from typing import Optional


class AdjudicateModel:
    """Choose the best gene model for overlapping features by family score.

    Compares two sets of gene predictions at overlapping loci and retains
    the feature with the higher-scoring gene family assignment.  Results are
    written to stdout (retained) and a TSV file (eliminated).

    Args:
        overlaps:      Path to BEDTools intersect file of overlapping features.
        genefamiliesa: Path to gene-family assignments for annotation set A.
        genefamiliesb: Path to gene-family assignments for annotation set B.
        minoverlap:    Minimum fractional overlap to trigger adjudication.
        no_orphans:    Skip features with no gene-family assignment when True.
        name:          Output prefix for the eliminated-features TSV.
    """

    def __init__(
        self,
        overlaps: str,
        genefamiliesa: str,
        genefamiliesb: str,
        minoverlap: float,
        no_orphans: bool,
        name: str,
    ) -> None:
        self.overlaps = overlaps
        self.genefamiliesa = genefamiliesa
        self.genefamiliesb = genefamiliesb
        self.minoverlap = minoverlap
        self.no_orphans = no_orphans
        self.name = name

        self.families: dict = {}
        self.my_ids: dict = {}
        self.alt_ids: dict = {}
        self.eliminated_ids: dict = {}

    def __repr__(self) -> str:
        return (
            f"AdjudicateModel(overlaps={self.overlaps!r}, "
            f"genefamiliesa={self.genefamiliesa!r}, "
            f"genefamiliesb={self.genefamiliesb!r}, "
            f"minoverlap={self.minoverlap}, "
            f"no_orphans={self.no_orphans}, "
            f"name={self.name!r})"
        )

    def get_families(self, genefamilies: str) -> dict:
        """Read a gene-families file and return a feature-to-family mapping.

        Args:
            genefamilies: Path to the tab-delimited gene-families file.

        Returns:
            Mapping of feature ID to gene-family, score, and model.
        """
        families: dict = {}
        with open(genefamilies) as fh:
            for raw in fh:
                line = raw.rstrip()
                if not line or line.startswith("#"):
                    continue
                fields = line.split("\t")
                feature = fields[0]
                genefamily = fields[1]
                model = fields[2]
                bscore = float(fields[5])
                score = float(fields[4])
                # Skip features with inflated scores from multi-copy / fusions
                if score / bscore >= 1.5:
                    continue
                families[feature] = {
                    "gene_family": genefamily,
                    "score": bscore,
                    "model": model,
                }
        return families

    def load_families(self) -> None:
        """Populate self.families for both annotation sets."""
        self.families = {
            "a": self.get_families(self.genefamiliesa),
            "b": self.get_families(self.genefamiliesb),
        }

    def choose_model(self) -> None:
        """Adjudicate overlapping features and write results."""
        self.load_families()

        with open(self.overlaps) as fh:
            for raw in fh:
                line = raw.rstrip()
                if not line:
                    continue
                self._process_line(line)

        self._write_results()

    def _process_line(self, line: str) -> None:
        """Parse one intersect line and update adjudication state."""
        fields = line.split("\t")
        overlap = int(fields[-1])
        featureid1 = fields[8].split(";")[0].split("=")[1]
        record1 = "\t".join(fields[:9])

        if not overlap:  # no overlap — keep unconditionally
            self._init_my_id(featureid1, record1)
            return

        record2 = "\t".join(fields[9:-1])
        featureid2 = fields[-2].split(";")[0].split("=")[1]
        start1, stop1 = int(fields[3]), int(fields[4])

        if overlap / (stop1 - start1) <= self.minoverlap:
            return

        self._init_my_id(featureid1, record1)

        family1, score1, model1 = self._lookup("a", featureid1)
        family2, score2, model2 = self._lookup("b", featureid2)

        if not family1 and not family2:
            self.eliminated_ids[featureid2] = self._elim_entry(
                featureid1, featureid2, score1, score2, "no gene families"
            )
            return

        if family1 == family2:
            self._adjudicate_same_family(
                featureid1,
                record1,
                family1,
                score1,
                model1,
                featureid2,
                record2,
                family2,
                score2,
                model2,
            )
        else:
            self._adjudicate_diff_family(
                featureid1,
                record1,
                family1,
                score1,
                model1,
                featureid2,
                record2,
                family2,
                score2,
                model2,
            )

    def _lookup(self, key: str, featureid: str) -> tuple:
        """Return (gene_family, score, model) for a feature, or defaults."""
        entry = self.families[key].get(featureid)
        if entry:
            return entry["gene_family"], entry["score"], entry["model"]
        return None, 0.0, None

    def _adjudicate_same_family(
        self,
        featureid1: str,
        record1: str,
        family1: Optional[str],
        score1: float,
        model1: Optional[str],
        featureid2: str,
        record2: str,
        family2: Optional[str],
        score2: float,
        model2: Optional[str],
    ) -> None:
        """Handle the case where both features share the same gene family."""
        current = self.my_ids[featureid1]["score"]

        if score1 >= score2 and current <= score1:
            self._update_my_id(
                featureid1, family1, score1, model1, featureid1, record1
            )
            reason = (
                "score1 = score2" if score1 == score2 else "score1 > score2"
            )
            self.eliminated_ids[featureid2] = self._elim_entry(
                featureid1, featureid2, score1, score2, reason
            )
        elif score2 > score1 and current < score2:
            self._update_my_id(
                featureid1, family2, score2, model2, featureid2, record2
            )
            self.eliminated_ids[featureid1] = self._elim_entry(
                featureid1, featureid2, score1, score2, "score1 < score2"
            )

    def _adjudicate_diff_family(
        self,
        featureid1: str,
        record1: str,
        family1: Optional[str],
        score1: float,
        model1: Optional[str],
        featureid2: str,
        record2: str,
        family2: Optional[str],
        score2: float,
        model2: Optional[str],
    ) -> None:
        """Handle the case where features belong to different gene families."""
        current = self.my_ids[featureid1]["score"]

        if not family1:  # only family2 has an assignment
            if current < score2:
                self._update_my_id(
                    featureid1, family2, score2, model2, featureid2, record2
                )
                self.eliminated_ids[featureid1] = self._elim_entry(
                    featureid1, featureid2, score1, score2, "score1 NA"
                )
        elif not family2:  # only family1 has an assignment
            if current <= score1:
                self._update_my_id(
                    featureid1, family1, score1, model1, featureid1, record1
                )
                self.eliminated_ids[featureid2] = self._elim_entry(
                    featureid1, featureid2, score1, score2, "score2 NA"
                )
        else:  # both features have assignments in different families
            if featureid2 not in self.eliminated_ids:
                if not family2 and self.no_orphans:
                    return
                self.my_ids[featureid2] = {
                    "family": family2,
                    "score": score2,
                    "model": model2,
                    "parent": featureid2,
                    "record": record2,
                }

    def _update_my_id(
        self,
        featureid: str,
        family: Optional[str],
        score: float,
        model: Optional[str],
        parent: str,
        record: str,
    ) -> None:
        """Overwrite all fields of an existing my_ids entry."""
        entry = self.my_ids[featureid]
        entry["score"] = score
        entry["family"] = family
        entry["model"] = model
        entry["parent"] = parent
        entry["record"] = record

    def _init_my_id(
        self,
        featureid: str,
        record: str,
        family: Optional[str] = None,
        score: float = 0.0,
        model: Optional[str] = None,
    ) -> None:
        """Add featureid to my_ids only if it is not already present."""
        if featureid not in self.my_ids:
            self.my_ids[featureid] = {
                "family": family,
                "score": score,
                "model": model,
                "parent": featureid,
                "record": record,
            }

    @staticmethod
    def _elim_entry(
        featureid1: str,
        featureid2: str,
        score1: float,
        score2: float,
        reason: str,
    ) -> str:
        """Format a tab-delimited eliminated-features log entry."""
        return f"{featureid1}\t{featureid2}\t{score1}\t{score2}\t{reason}"

    def _write_results(self) -> None:
        """Print retained records to stdout; write eliminated ones to TSV."""
        for data in self.my_ids.values():
            if data["parent"] not in self.eliminated_ids:
                print(data["record"])

        out_path = f"./{self.name}.eliminated.tsv"
        with open(out_path, "w") as fh:
            for entry in self.eliminated_ids.values():
                fh.write(f"{entry}\n")


if __name__ == "__main__":
    _intersect = sys.argv[1]
    _genefamiliesa = sys.argv[2]
    _genefamiliesb = sys.argv[3]
    _minoverlap = float(sys.argv[4])
    _no_orphans = sys.argv[5]
    _name = sys.argv[6]

    adjudicator = AdjudicateModel(
        _intersect,
        _genefamiliesa,
        _genefamiliesb,
        _minoverlap,
        _no_orphans,
        _name,
    )
    adjudicator.choose_model()
