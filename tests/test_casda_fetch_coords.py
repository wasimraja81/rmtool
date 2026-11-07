#!/usr/bin/env python3
"""Value-level tests for scripts/casda_fetch.py's coordinate
notation/frame support (--ra/--dec/--lon/--lat/--frame) -- the pure,
network-free parts of the script (resolve_target's own parsing,
describe_coord, parse_args' validation). Everything else in
casda_fetch.py talks to CASDA's live service and has no test coverage
for that reason; this covers what can be tested without it.

Run directly: python3 tests/test_casda_fetch_coords.py
"""
from __future__ import annotations
import argparse
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
SCRIPT = REPO_ROOT / "scripts" / "casda_fetch.py"
sys.path.insert(0, str(REPO_ROOT / "scripts"))
import casda_fetch as cf  # noqa: E402

ok = True


def report(flag_ok: bool, msg: str):
    global ok
    print(f"[{'OK' if flag_ok else 'FAIL'}] {msg}")
    if not flag_ok:
        ok = False


def make_args(**kw):
    base = dict(username="u", ra=None, dec=None, lon=None, lat=None,
                frame="icrs", target=None, abort=None, run_mode="dry",
                select=None)
    base.update(kw)
    return argparse.Namespace(**base)


def close(a: float, b: float, tol: float = 1e-6) -> bool:
    return abs(a - b) < tol


def check_resolve_target():
    # Decimal degrees unchanged from the script's pre-existing behaviour.
    c = cf.resolve_target(make_args(ra="323.56667", dec="-53.61528"))
    report(c.frame.name == "icrs" and close(c.ra.deg, 323.56667)
           and close(c.dec.deg, -53.61528),
           "decimal degrees --ra/--dec (default icrs) unchanged")

    # Colon-separated sexagesimal RA/Dec -> same position as above.
    c = cf.resolve_target(make_args(ra="21:34:16", dec="-53:36:55"))
    report(close(c.ra.deg, 323.56667, tol=1e-3) and close(c.dec.deg, -53.61528, tol=1e-3),
           "colon-sexagesimal --ra/--dec matches the decimal-degree equivalent")

    # Unit-lettered sexagesimal -> same position again.
    c = cf.resolve_target(make_args(ra="21h34m16s", dec="-53d36m55s"))
    report(close(c.ra.deg, 323.56667, tol=1e-3) and close(c.dec.deg, -53.61528, tol=1e-3),
           "unit-lettered sexagesimal --ra/--dec matches the decimal-degree equivalent")

    # fk5/fk4 pass-through, no implicit conversion to icrs at this stage
    # (astroquery's own Casda class does that later, not this script).
    c = cf.resolve_target(make_args(ra="180.5", dec="-45.2", frame="fk5"))
    report(c.frame.name == "fk5" and close(c.ra.deg, 180.5) and close(c.dec.deg, -45.2),
           "--frame=fk5 stays fk5 (no premature frame conversion)")
    c = cf.resolve_target(make_args(ra="180.5", dec="-45.2", frame="fk4"))
    report(c.frame.name == "fk4" and close(c.ra.deg, 180.5) and close(c.dec.deg, -45.2),
           "--frame=fk4 stays fk4")

    # Galactic via --lon/--lat, bare decimal.
    c = cf.resolve_target(make_args(lon="308.5", lat="-40.2", frame="galactic"))
    report(c.frame.name == "galactic" and close(c.l.deg, 308.5) and close(c.b.deg, -40.2),
           "galactic --lon/--lat, decimal degrees")

    # Galactic via --lon/--lat, colon-sexagesimal -- must resolve as
    # DEGREES, not hour angle (the one case this whole feature exists
    # to get right: Galactic longitude has no hour-angle convention).
    c = cf.resolve_target(make_args(lon="308:30:00", lat="-40:12:00", frame="galactic"))
    report(close(c.l.deg, 308.5, tol=1e-6) and close(c.b.deg, -40.2, tol=1e-6),
           "galactic colon-sexagesimal --lon/--lat resolves as degrees, not hour angle")

    # Malformed coordinate string -- fails loudly, not silently.
    try:
        cf.resolve_target(make_args(ra="not-a-coordinate", dec="-45"))
        report(False, "malformed --ra string should raise, not resolve silently")
    except (ValueError, Exception):
        report(True, "malformed --ra string raises instead of resolving silently")


def check_describe_coord():
    for frame, kw in (("icrs", dict(ra="323.56667", dec="-53.61528")),
                       ("fk5", dict(ra="180.5", dec="-45.2", frame="fk5")),
                       ("fk4", dict(ra="180.5", dec="-45.2", frame="fk4")),
                       ("galactic", dict(lon="308.5", lat="-40.2", frame="galactic"))):
        c = cf.resolve_target(make_args(**kw))
        desc = cf.describe_coord(c)
        report(bool(desc) and "ICRS" in desc,
               f"describe_coord does not crash for --frame={frame} and names ICRS")


def run_cli(*cli_args: str) -> subprocess.CompletedProcess:
    return subprocess.run(
        [sys.executable, str(SCRIPT), *cli_args],
        capture_output=True, text=True, timeout=30)


def check_cli_validation():
    cases = [
        (["--username", "u", "--ra", "180"], "--ra given without --dec"),
        (["--username", "u", "--dec", "-45"], "--dec given without --ra"),
        (["--username", "u", "--lon", "180"], "--lon given without --lat"),
        (["--username", "u", "--lat", "-45"], "--lat given without --lon"),
        (["--username", "u", "--ra", "180", "--dec", "-45", "--frame", "galactic"],
         "--ra/--dec are equatorial-only"),
        (["--username", "u", "--ra", "180", "--dec", "-45", "--lat", "-40"],
         "give either --ra/--dec or --lon/--lat, not both"),
        (["--username", "u"],
         "one of --target, --ra/--dec, or --lon/--lat is required"),
        (["--username", "u", "--target", "foo", "--ra", "180"],
         "not allowed with argument --target"),
        (["--username", "u", "--target", "foo", "--lon", "180"],
         "not allowed with argument --target"),
    ]
    for cli_args, expect_substr in cases:
        proc = run_cli(*cli_args)
        report(proc.returncode != 0 and expect_substr in proc.stderr,
               f"CLI rejects {cli_args} with '{expect_substr}'")

    # Positive cases: each should parse successfully (argparse exit code
    # 0 from parse_args alone is not observable without running main(),
    # which would hit the network -- instead confirm parse_args() itself
    # succeeds in-process, same function the CLI subprocess runs).
    positive_cases = [
        dict(ra="323.56667", dec="-53.61528"),
        dict(ra="21:34:16", dec="-53:36:55"),
        dict(ra="21h34m16s", dec="-53d36m55s"),
        dict(lon="308.5", lat="-40.2", frame="galactic"),
        dict(ra="180.5", dec="-45.2", frame="fk5"),
    ]
    for kw in positive_cases:
        argv = ["casda_fetch.py", "--username", "u"]
        for key, val in kw.items():
            argv.append(f"--{key.replace('_', '-')}")
            argv.append(str(val))
        old_argv = sys.argv
        sys.argv = argv
        try:
            cf.parse_args()
            report(True, f"parse_args() accepts {kw}")
        except SystemExit:
            report(False, f"parse_args() unexpectedly rejected {kw}")
        finally:
            sys.argv = old_argv


def main():
    check_resolve_target()
    check_describe_coord()
    check_cli_validation()
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
