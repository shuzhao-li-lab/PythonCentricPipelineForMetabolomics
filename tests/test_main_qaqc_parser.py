"""Focused tests for QAQC command-line selector parsing."""

import copy
import importlib.util
import sys
import types
from collections import OrderedDict
from pathlib import Path
from unittest.mock import Mock

import pytest


@pytest.fixture
def main_module(monkeypatch: pytest.MonkeyPatch):
    """Load ``pcpfm.main`` without requiring optional pipeline dependencies."""

    package_path = Path(__file__).resolve().parents[1] / "pcpfm"
    package = types.ModuleType("pcpfm")
    package.__path__ = [str(package_path)]
    monkeypatch.setitem(sys.modules, "pcpfm", package)

    monkeypatch.setitem(sys.modules, "gdown", types.ModuleType("gdown"))
    for dependency_name in ("Experiment", "EmpCpds", "Report"):
        dependency = types.ModuleType(f"pcpfm.{dependency_name}")
        monkeypatch.setitem(sys.modules, dependency.__name__, dependency)
        setattr(package, dependency_name, dependency)

    defaults_spec = importlib.util.spec_from_file_location(
        "pcpfm.default_parameters", package_path / "default_parameters.py"
    )
    assert defaults_spec is not None and defaults_spec.loader is not None
    defaults = importlib.util.module_from_spec(defaults_spec)
    monkeypatch.setitem(sys.modules, "pcpfm.default_parameters", defaults)
    defaults_spec.loader.exec_module(defaults)
    setattr(package, "default_parameters", defaults)

    main_spec = importlib.util.spec_from_file_location(
        "pcpfm.main", package_path / "main.py"
    )
    assert main_spec is not None and main_spec.loader is not None
    main_module = importlib.util.module_from_spec(main_spec)
    monkeypatch.setitem(sys.modules, "pcpfm.main", main_module)
    main_spec.loader.exec_module(main_module)
    return main_module


@pytest.fixture
def feature_table_module(monkeypatch: pytest.MonkeyPatch):
    """Load ``pcpfm.FeatureTable`` while avoiding unrelated package imports."""

    package_path = Path(__file__).resolve().parents[1] / "pcpfm"
    package = types.ModuleType("pcpfm")
    package.__path__ = [str(package_path)]
    monkeypatch.setitem(sys.modules, "pcpfm", package)
    monkeypatch.setitem(sys.modules, "pcpfm.utils", types.ModuleType("pcpfm.utils"))

    intervaltree = types.ModuleType("intervaltree")
    intervaltree.IntervalTree = type("IntervalTree", (), {})
    monkeypatch.setitem(sys.modules, "intervaltree", intervaltree)
    combat = types.ModuleType("combat")
    pycombat = types.ModuleType("combat.pycombat")
    pycombat.pycombat = lambda *args, **kwargs: None
    monkeypatch.setitem(sys.modules, "combat", combat)
    monkeypatch.setitem(sys.modules, "combat.pycombat", pycombat)

    feature_table_spec = importlib.util.spec_from_file_location(
        "pcpfm.FeatureTable", package_path / "FeatureTable.py"
    )
    assert feature_table_spec is not None and feature_table_spec.loader is not None
    feature_table_module = importlib.util.module_from_spec(feature_table_spec)
    monkeypatch.setitem(sys.modules, "pcpfm.FeatureTable", feature_table_module)
    feature_table_spec.loader.exec_module(feature_table_module)
    return feature_table_module


def test_pca_selector_merges_with_real_qaqc_defaults(
    monkeypatch: pytest.MonkeyPatch, main_module
) -> None:
    monkeypatch.setattr(
        main_module.default_parameters,
        "PARAMETERS",
        copy.deepcopy(main_module.default_parameters.PARAMETERS),
    )
    monkeypatch.setattr(sys, "argv", ["pcpfm", "QAQC", "--pca"])

    params = main_module.Main.process_params()

    assert params["pca"] is True
    assert params["all"] is False
    assert params["tsne"] is False
    assert params["pearson"] is False
    assert params["spearman"] is False
    assert params["kendall"] is False
    assert params["missing_feature_distribution"] is False
    assert params["missing_feature_percentiles"] is False
    assert params["median_correlation_outlier_detection"] is False
    assert params["missing_feature_outlier_detection"] is False
    assert params["intensity_analysis"] is False
    assert params["feature_distribution"] is False
    assert params["feature_outlier_detection"] is False
    assert params["interactive_plots"] is False
    assert params["save_plots"] is True


def test_qaqc_without_selector_retains_run_all_default(
    monkeypatch: pytest.MonkeyPatch, main_module
) -> None:
    monkeypatch.setattr(
        main_module.default_parameters,
        "PARAMETERS",
        copy.deepcopy(main_module.default_parameters.PARAMETERS),
    )
    monkeypatch.setattr(sys, "argv", ["pcpfm", "QAQC"])

    params = main_module.Main.process_params()

    assert params["all"] is True


def test_qaqc_all_selector_retains_run_all_default(
    monkeypatch: pytest.MonkeyPatch, main_module
) -> None:
    monkeypatch.setattr(
        main_module.default_parameters,
        "PARAMETERS",
        copy.deepcopy(main_module.default_parameters.PARAMETERS),
    )
    monkeypatch.setattr(sys, "argv", ["pcpfm", "QAQC", "--all"])

    params = main_module.Main.process_params()

    assert params["all"] is True


def test_qaqc_pca_dispatches_only_pca(
    monkeypatch: pytest.MonkeyPatch, main_module, feature_table_module
) -> None:
    monkeypatch.setattr(
        main_module.default_parameters,
        "PARAMETERS",
        copy.deepcopy(main_module.default_parameters.PARAMETERS),
    )
    monkeypatch.setattr(sys, "argv", ["pcpfm", "QAQC", "--pca"])
    params = main_module.Main.process_params()

    feature_table = feature_table_module.FeatureTable(None, None, "test")
    feature_table.generate_figure_params = Mock()
    feature_table.method_map = OrderedDict(
        (name, Mock(return_value={"Type": name})) for name in feature_table.method_map
    )

    results = feature_table.QAQC(params)

    assert results == [{"Type": "pca"}]
    assert feature_table.method_map["pca"].call_count == 1
    assert all(
        method.call_count == 0
        for name, method in feature_table.method_map.items()
        if name != "pca"
    )


def test_pca_selector_rejects_a_following_value(
    monkeypatch: pytest.MonkeyPatch, main_module
) -> None:
    monkeypatch.setattr(main_module.default_parameters, "PARAMETERS", {"multicores": 1})
    monkeypatch.setattr(sys, "argv", ["pcpfm", "QAQC", "--pca", "unexpected"])

    with pytest.raises(SystemExit) as exc_info:
        main_module.Main.process_params()

    assert exc_info.value.code == 2
