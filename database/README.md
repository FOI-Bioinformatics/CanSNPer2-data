# Databases

These SQLite files are generated from the text files under `source/` by
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
need rebuilding without building them.

`uncommitted_sources` lists any of that scheme's source files that had
uncommitted changes when the database was built. A non-empty list means
`source_commit` does not fully describe the build.

## Archive

`archive/` holds the database each rebuild replaced, named
`<scheme>_<YYMMDD>.db`. The three dated Francisella files predate this
convention being automated but follow the same pattern.

## Using a database

```sh
CanSNPer2-download --database database/francisella_tularensis.db
CanSNPer2 --database database/francisella_tularensis.db sample.fasta --summary
```

The `genomes` table is empty in these databases, which is expected:
`CanSNPer2-download` reads reference genomes from `snp_references`.
