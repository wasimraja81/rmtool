"""ICRS<->FK4 position conversion, independent of astropy's own FK4/FK5
transform graph (which differs from AST by up to ~47 mas -- astropy's
frame classes omit SOFA's eraHfk5z epoch-dependent ICRS<->FK5 term,
documented in docs/dev/REPROJECT_EXACT_PLAN.md sections 5d/5e). This is
the single, shared implementation: a from-source reimplementation of
AST's own published PAL/SOFA algorithm (palFk524.c/palFk45z.c, AST's
bundled pal/ source -- not a call to AST, not astropy's transform_to),
position-only (zero proper motion/parallax/radial velocity, matching
AST's HFK5Z/FK54Z/FK45Z path). Verified against AST's compiled output
on 480 whole-sky points: 0.0000 mas for ICRS->FK4, sub-0.02 mas
roundtrip (REPROJECT_EXACT_PLAN.md section 5f).

Used by tests/check_target_projection.py (B1, the committed FK4
coverage) and scratch/exact/largeimage.py (the UR-24 matrix harness) --
one implementation, not two, per the project's standing no-duplication
rule. scratch/ is gitignored, so this module lives here, in the
committed tree, and the harness imports it from here rather than
keeping its own copy.
"""
import numpy as np
import erfa

_FK524_EMI = np.array([[0.9999256795, 0.0111814828, 0.0048590039],
                       [-0.0111814828, 0.9999374849, -0.0000271771],
                       [-0.0048590040, -0.0000271557, 0.9999881946]])
_FK524_A = np.array([-1.62557e-6, -0.31919e-6, -0.13843e-6])
_FK45Z_EM = np.array([[0.9999256782, -0.0111820611, -0.0048579477],
                      [0.0111820610, 0.9999374784, -0.0000271765],
                      [0.0048579479, -0.0000271474, 0.9999881997],
                      [-0.000551, -0.238565, 0.435739],
                      [0.238514, -0.002667, -0.008541],
                      [-0.435623, 0.012254, 0.002117]])
_FK45Z_AD = np.array([1.245e-3, -1.580e-3, -0.659e-3])
_PAL_PMF = 100.0 * 60.0 * 60.0 * 360.0 / (2.0 * np.pi)


def fk524_pos(r2000, d2000):
    """palFk524.c, position-only: FK5 J2000.0 -> FK4 B1950.0. r2000/d2000
    in radians, any broadcastable shape."""
    x = np.cos(r2000) * np.cos(d2000)
    y = np.sin(r2000) * np.cos(d2000)
    z = np.sin(d2000)
    v2 = np.stack([x, y, z], axis=-1) @ _FK524_EMI.T
    rxyz = np.linalg.norm(v2, axis=-1)
    w = (v2 * _FK524_A).sum(axis=-1)
    pos = v2 + _FK524_A * rxyz[..., None] - w[..., None] * v2
    rxy = np.hypot(pos[..., 0], pos[..., 1])
    r1950 = np.mod(np.arctan2(pos[..., 1], pos[..., 0]), 2 * np.pi)
    d1950 = np.arctan2(pos[..., 2], rxy)
    return r1950, d1950


def fk45z_pos(r1950, d1950, bepoch):
    """palFk45z.c, zero proper motion in FK5: FK4 B1950.0 at Besselian
    epoch bepoch -> FK5 J2000.0. r1950/d1950 in radians."""
    r0 = np.stack([np.cos(r1950) * np.cos(d1950), np.sin(r1950) * np.cos(d1950),
                   np.sin(d1950)], axis=-1)
    a1 = _FK524_A + (bepoch - 1950.0) / _PAL_PMF * _FK45Z_AD
    wdot = (r0 * a1).sum(axis=-1)
    v1 = r0 - a1 + wdot[..., None] * r0
    v2_full = v1 @ _FK45Z_EM.T
    epj = erfa.epj(*erfa.epb2jd(bepoch))
    v2 = v2_full[..., 0:3] + (epj - 2000.0) / _PAL_PMF * v2_full[..., 3:6]
    rxy = np.hypot(v2[..., 0], v2[..., 1])
    r2000 = np.mod(np.arctan2(v2[..., 1], v2[..., 0]), 2 * np.pi)
    d2000 = np.arctan2(v2[..., 2], rxy)
    return r2000, d2000


def icrs_to_fk4(ra_icrs_deg, dec_icrs_deg, equinox_byear=1950.0):
    """ICRS (degrees) -> FK4 at the given Besselian equinox (degrees),
    via eraHfk5z then fk524_pos -- AST's own ICRS<->FK4 path."""
    djm0, djm = erfa.epb2jd(equinox_byear)
    r5, d5, _, _ = erfa.hfk5z(np.radians(ra_icrs_deg), np.radians(dec_icrs_deg), djm0, djm)
    r4, d4 = fk524_pos(r5, d5)
    return np.degrees(r4) % 360.0, np.degrees(d4)


def fk4_to_icrs(ra_fk4_deg, dec_fk4_deg, equinox_byear=1950.0):
    """FK4 at the given Besselian equinox (degrees) -> ICRS (degrees),
    via fk45z_pos then eraFk5hz -- the inverse of icrs_to_fk4."""
    djm0, djm = erfa.epb2jd(equinox_byear)
    r5, d5 = fk45z_pos(np.radians(ra_fk4_deg), np.radians(dec_fk4_deg), equinox_byear)
    rh, dh = erfa.fk5hz(r5, d5, djm0, djm)
    return np.degrees(rh) % 360.0, np.degrees(dh)
