# Databases

These SQLite files are generated from the text files under `schemes/` by
`build.py`. They are not edited by hand: a change belongs in the source files,
after which the database is rebuilt.

```sh
python build.py           # rebuild and compare, changing nothing
python build.py --write   # install the rebuilt databases
```

A rebuild only replaces a database whose content has actually changed. The
SQLite bytes differ between two builds of the same scheme, because
`snp_annotation.date` records the build date, so replacing an unchanged database
would produce a diff that means nothing.

## Manifests

Each database has a `<scheme>.manifest.json` beside it recording what produced
it: the commit the sources were at, the SHA-256 of each source file, the
installed CanSNPer2 and Python versions, and the row count of every table. The manifest is text, so its diff is reviewable where the database's is
not, and the recorded hashes let `build.py --changed-only` tell which schemes
need rebuilding without building them. `python build.py` fails when a manifest
does not record the current sources, which catches a source edit that does not
change the database, such as the strain column of a SNP file.

`build.py --write` leaves the manifest of a database it did not replace as it
is, so `built`, `versions` and `counts` keep describing the build that produced
the committed file. If only the source hashes are out of date, those are updated.

`uncommitted_sources` lists any of that scheme's source files that had
uncommitted changes when the database was built. A non-empty list means
`source_commit` does not fully describe the build.

## Archive

`archive/` holds the database each rebuild replaced, named
`<scheme>_<YYMMDD>.db`. The Francisella files dated 2020 predate this
convention being automated but follow the same pattern.

The databases built before 2026-11-21 count SNP positions along the reference
sequences joined in file order. They are read by CanSNPer2 2.2 and earlier;
CanSNPer2 2.3.0 refuses them unless they are converted with
`CanSNPer2-database --convert --database <file> --refdir <folder> -o <outdir>`.

## Using a database

```sh
CanSNPer2-download --database database/francisella_tularensis.db
CanSNPer2 --database database/francisella_tularensis.db sample.fasta --summary
```

The `genomes` table is empty in these databases, which is expected:
`CanSNPer2-download` reads reference genomes from `snp_references`.

## Changing a database directly

`CanSNPer2-database --mod_file` can change the tree of an existing database
without rebuilding it. This repository does not use it: a database changed that
way no longer matches the scheme files it was built from, and `build.py` reports
the difference. Change the scheme files and run `python build.py --write`
instead.
