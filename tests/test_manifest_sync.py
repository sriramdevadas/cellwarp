"""
Test that environment.yml and pyproject.toml [lock] declare identical
version pins for every package present in both.

Drift between the two manifests would mean conda installers and pip
installers reproduce against different dependency versions — which
defeats the whole point of having a [lock] group.
"""
from __future__ import annotations

import re
import tomllib
from pathlib import Path

import pytest
import yaml

REPO_ROOT = Path(__file__).resolve().parent.parent
PYPROJECT = REPO_ROOT / "pyproject.toml"
ENVIRONMENT_YML = REPO_ROOT / "environment.yml"


def _parse_pyproject_lock() -> dict[str, str]:
    """Return {package_name_normalized: version} from pyproject.toml [lock]."""
    with open(PYPROJECT, "rb") as f:
        data = tomllib.load(f)
    lock = data["project"]["optional-dependencies"]["lock"]
    pins: dict[str, str] = {}
    for entry in lock:
        # entry like "numpy==2.4.3"
        m = re.match(r"^([A-Za-z0-9_.\-]+)==([A-Za-z0-9_.\-+]+)$", entry.strip())
        assert m, f"pyproject.toml [lock] entry not in 'name==version' form: {entry!r}"
        name = m.group(1).lower().replace("_", "-")
        version = m.group(2)
        pins[name] = version
    return pins


def _parse_environment_yml() -> dict[str, str]:
    """
    Return {package_name_normalized: version} from environment.yml.
    Combines conda deps and pip deps. Skips the `python` entry and
    skips any local install directives (e.g. `-e .`).
    """
    with open(ENVIRONMENT_YML) as f:
        env = yaml.safe_load(f)
    pins: dict[str, str] = {}
    for dep in env.get("dependencies", []):
        if isinstance(dep, str):
            # conda dep, e.g. "numpy=2.4.3" or "python=3.12.12"
            if "=" not in dep:
                continue
            name, _, version = dep.partition("=")
            name = name.lower().replace("_", "-")
            if name == "python" or name == "pip":
                continue
            pins[name] = version
        elif isinstance(dep, dict) and "pip" in dep:
            for pip_entry in dep["pip"]:
                pip_entry = pip_entry.strip()
                if pip_entry.startswith("-e ") or pip_entry == "-e .":
                    continue
                m = re.match(
                    r"^([A-Za-z0-9_.\-]+)==([A-Za-z0-9_.\-+]+)$", pip_entry
                )
                assert m, (
                    f"environment.yml pip entry not in 'name==version' form: "
                    f"{pip_entry!r}"
                )
                name = m.group(1).lower().replace("_", "-")
                version = m.group(2)
                pins[name] = version
    return pins


def test_manifests_exist() -> None:
    assert PYPROJECT.exists(), f"missing: {PYPROJECT}"
    assert ENVIRONMENT_YML.exists(), f"missing: {ENVIRONMENT_YML}"


def test_lock_and_environment_versions_agree() -> None:
    """
    For every package present in BOTH pyproject.toml [lock] and
    environment.yml, the pinned version must be identical.
    """
    lock_pins = _parse_pyproject_lock()
    env_pins = _parse_environment_yml()

    common = set(lock_pins) & set(env_pins)
    assert common, (
        "no overlap between pyproject.toml [lock] and environment.yml — "
        "one of the manifests is empty or malformed"
    )

    mismatches: list[str] = []
    for pkg in sorted(common):
        if lock_pins[pkg] != env_pins[pkg]:
            mismatches.append(
                f"{pkg}: pyproject.toml=={lock_pins[pkg]} vs "
                f"environment.yml={env_pins[pkg]}"
            )
    assert not mismatches, (
        "version drift between pyproject.toml [lock] and environment.yml:\n  "
        + "\n  ".join(mismatches)
    )


def test_every_lock_pkg_appears_in_environment() -> None:
    """
    Every package in pyproject.toml [lock] must appear somewhere in
    environment.yml (either conda deps or pip deps). This catches the
    failure mode where a package is added to [lock] but the conda
    install path silently omits it.
    """
    lock_pins = _parse_pyproject_lock()
    env_pins = _parse_environment_yml()
    missing = sorted(set(lock_pins) - set(env_pins))
    assert not missing, (
        "packages in pyproject.toml [lock] but absent from environment.yml: "
        + ", ".join(missing)
    )


def test_requirements_txt_pins_match_lock():
    """Every package in requirements.txt that has an exact == pin must match
    the corresponding pin in pyproject.toml [lock] extras. requirements.txt is
    allowed to contain extra transitive packages not in [lock]; the test only
    checks the overlap.
    """
    import re
    from pathlib import Path

    try:
        import tomllib
    except ImportError:
        import tomli as tomllib

    project_root = Path(__file__).resolve().parent.parent
    requirements_txt = project_root / "requirements.txt"
    pyproject_toml = project_root / "pyproject.toml"

    assert requirements_txt.exists(), "requirements.txt missing at repo root"
    assert pyproject_toml.exists(), "pyproject.toml missing at repo root"

    # Parse [lock] extras
    with open(pyproject_toml, "rb") as f:
        pyproject = tomllib.load(f)
    lock_extras = (
        pyproject.get("project", {})
        .get("optional-dependencies", {})
        .get("lock", [])
    )
    lock_pins = {}
    for entry in lock_extras:
        m = re.match(r"^([A-Za-z0-9_.\-]+)==([A-Za-z0-9_.\-+]+)$", entry.strip())
        if m:
            lock_pins[m.group(1).lower()] = m.group(2)

    # Parse requirements.txt (skip comments, blank lines, editable installs)
    req_pins = {}
    for line in requirements_txt.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#") or line.startswith("-e"):
            continue
        m = re.match(r"^([A-Za-z0-9_.\-]+)==([A-Za-z0-9_.\-+]+)", line)
        if m:
            req_pins[m.group(1).lower()] = m.group(2)

    # Compare overlap
    mismatches = []
    for pkg, lock_version in lock_pins.items():
        if pkg in req_pins and req_pins[pkg] != lock_version:
            mismatches.append(
                f"{pkg}: lock=={lock_version} vs requirements.txt=={req_pins[pkg]}"
            )

    assert not mismatches, "Pin drift between [lock] and requirements.txt:\n" + "\n".join(mismatches)


def test_every_lock_pkg_appears_in_requirements_txt():
    """Every package in pyproject.toml [lock] extras must appear in
    requirements.txt (any version specifier). Catches the failure mode where
    a new pin is added to [lock] but not propagated to requirements.txt.
    """
    import re
    from pathlib import Path

    try:
        import tomllib
    except ImportError:
        import tomli as tomllib

    project_root = Path(__file__).resolve().parent.parent
    requirements_txt = project_root / "requirements.txt"
    pyproject_toml = project_root / "pyproject.toml"

    with open(pyproject_toml, "rb") as f:
        pyproject = tomllib.load(f)
    lock_extras = (
        pyproject.get("project", {})
        .get("optional-dependencies", {})
        .get("lock", [])
    )
    lock_pkgs = set()
    for entry in lock_extras:
        m = re.match(r"^([A-Za-z0-9_.\-]+)==", entry.strip())
        if m:
            lock_pkgs.add(m.group(1).lower())

    req_content = requirements_txt.read_text().lower()

    missing = sorted(p for p in lock_pkgs if p not in req_content)
    assert not missing, f"[lock] packages absent from requirements.txt: {missing}"
