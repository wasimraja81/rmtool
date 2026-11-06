# Demo data

Real, public cutouts used by `scripts/demo_reproject_gallery.py` and
[docs/user/EXAMPLES.md](../EXAMPLES.md#8-visual-gallery-what-reprojection-does-to-the-sky).
Both are small float32 crops, not the full surveys, fetched via CDS's
`hips2fits` cutout service (`http://alasky.cds.unistra.fr/hips-image-services/hips2fits`).
`coordsys=galactic` only relabels the output `CTYPE` axes -- it does
not reinterpret `ra=`/`dec=` as galactic `l`/`b` -- so both fetches
below pass the ICRS coordinates of `l=0,b=0`
(`ra=266.4049882865447,dec=-28.936177761791473`, from
`SkyCoord(l=0*u.deg, b=0*u.deg, frame='galactic').icrs`), not `ra=0,dec=0`.

## `haslam_408MHz_galactic.fits`

Haslam 408 MHz all-sky map, Remazeilles et al. 2014 reprocessing
(`CDS/P/Haslam408`), whole sky, galactic/AIT:

```bash
curl -s -o haslam_408MHz_galactic.fits \
  "http://alasky.cds.unistra.fr/hips-image-services/hips2fits?hips=CDS/P/Haslam408&width=1800&height=900&projection=AIT&coordsys=galactic&ra=266.4049882865447&dec=-28.936177761791473&fov=360&format=fits"
```

## `meerkat_gc_1284MHz_galactic.fits`

MeerKAT 1.28 GHz Galactic Centre mosaic, Heywood et al. 2022 initial
data release (`CDS/P/MeerKAT/Galactic-Centre-1284MHz-StokesI`), a 1.5
deg cutout around Sgr A*, galactic/SIN:

```bash
curl -s -o meerkat_gc_1284MHz_galactic.fits \
  "http://alasky.cds.unistra.fr/hips-image-services/hips2fits?hips=CDS/P/MeerKAT/Galactic-Centre-1284MHz-StokesI&width=1500&height=1000&projection=SIN&coordsys=galactic&ra=266.4049882865447&dec=-28.936177761791473&fov=1.5&format=fits"
```
