#!/usr/bin/env python3
from __future__ import annotations

import csv
from pathlib import Path


def load_drb1_dqb1_map(path: Path) -> tuple[dict[str, str], dict[str, str]]:
    drb1_to_dqb1: dict[str, str] = {}
    dqb1_to_drb1: dict[str, str] = {}
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"drb1", "dqb1"}
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise ValueError(f"{path} must contain tab-delimited columns: drb1, dqb1")
        for row in reader:
            drb1 = (row.get("drb1") or "").strip()
            dqb1 = (row.get("dqb1") or "").strip()
            if not drb1 or not dqb1:
                continue
            drb1_to_dqb1[drb1] = dqb1
            dqb1_to_drb1[dqb1] = drb1
    return drb1_to_dqb1, dqb1_to_drb1