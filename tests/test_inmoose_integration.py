"""Integration tests for PCPFM functionality backed by InMoose."""

from __future__ import annotations

import importlib.util
import sys
import types
from pathlib import Path

import numpy as np
import pandas as pd
import pytest


@pytest.fixture
def feature_table_module(monkeypatch: pytest.MonkeyPatch):
    """Load FeatureTable with real InMoose and minimal unrelated imports."""
    package_path = Path(__file__).resolve().parents[1] / "pcpfm"
    package = types.ModuleType("pcpfm")
    package.__path__ = [str(package_path)]
    monkeypatch.setitem(sys.modules, "pcpfm", package)
    monkeypatch.setitem(sys.modules, "pcpfm.utils", types.ModuleType("pcpfm.utils"))

    feature_table_spec = importlib.util.spec_from_file_location(
        "pcpfm.FeatureTable", package_path / "FeatureTable.py"
    )
    assert feature_table_spec is not None and feature_table_spec.loader is not None
    feature_table_module = importlib.util.module_from_spec(feature_table_spec)
    monkeypatch.setitem(sys.modules, "pcpfm.FeatureTable", feature_table_module)
    feature_table_spec.loader.exec_module(feature_table_module)
    return feature_table_module


def test_batch_correct_reduces_a_deterministic_batch_shift(
    feature_table_module,
) -> None:
    sample_names = ["batch1_a", "batch1_b", "batch2_a", "batch2_b"]

    class Acquisition:
        def __init__(self, name: str) -> None:
            self.name = name

    class Experiment:
        def __init__(self) -> None:
            self.acquisitions = [Acquisition(name) for name in sample_names]

        @staticmethod
        def batches(field: str) -> dict[str, list[str]]:
            assert field == "Batch"
            return {
                "batch1": ["batch1_a", "batch1_b"],
                "batch2": ["batch2_a", "batch2_b"],
            }

    feature_index = np.arange(20, dtype=float)
    baseline = 100.0 + 10.0 * feature_index
    batch_shift = 25.0 + 3.0 * (feature_index % 5)
    batch1_spread = feature_index % 3 + 1.0
    batch2_spread = 1.5 * (feature_index % 4 + 1.0)
    feature_metadata = {
        "id_number": np.arange(1, 21),
        **{f"metadata_{index}": np.full(20, index) for index in range(1, 11)},
    }
    table = feature_table_module.FeatureTable(
        pd.DataFrame(
            {
                **feature_metadata,
                "batch1_a": baseline - batch1_spread,
                "batch1_b": baseline + batch1_spread,
                "batch2_a": baseline + batch_shift - batch2_spread,
                "batch2_b": baseline + batch_shift + batch2_spread,
            }
        ),
        Experiment(),
        "input",
    )
    original_samples = table.feature_table[sample_names].copy()
    original_metadata = table.feature_table[list(feature_metadata)].copy()
    batch_gap_before = (
        original_samples[["batch2_a", "batch2_b"]].mean(axis=1)
        - original_samples[["batch1_a", "batch1_b"]].mean(axis=1)
    ).abs()
    assert batch_gap_before.tolist() == [25.0, 28.0, 31.0, 34.0, 37.0] * 4

    table.batch_correct("Batch")

    corrected_samples = table.feature_table[sample_names]
    batch_gap_after = (
        corrected_samples[["batch2_a", "batch2_b"]].mean(axis=1)
        - corrected_samples[["batch1_a", "batch1_b"]].mean(axis=1)
    ).abs()
    assert np.isfinite(corrected_samples.to_numpy()).all()
    assert not corrected_samples.equals(original_samples)
    pd.testing.assert_frame_equal(
        table.feature_table[list(feature_metadata)], original_metadata
    )
    assert (batch_gap_after < batch_gap_before * 0.10).all()
    np.testing.assert_allclose(
        corrected_samples.mean(axis=1),
        original_samples.mean(axis=1),
        rtol=0,
        atol=0.25,
    )
