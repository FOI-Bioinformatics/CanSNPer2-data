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
it: `source_commit`, the last commit that changed the source files of the
scheme; the SHA-256 of each source file; `database_sha256`, the SHA-256 of the
database file; `format`, the database format CanSNPer2 records; the installed
CanSNPer2 and Python versions; and the row count of every table. The manifest is
text, so its diff is reviewable where the database's is not, and the recorded
hashes let `build.py --changed-only` tell which schemes need rebuilding without
building them. `python build.py` fails when a manifest does not record the
current sources, which catches a source edit that does not change the database,
such as the strain column of a SNP file, and when `database_sha256` does not
match the committed database.

`build.py --write` refuses a scheme whose source files, README or `schemes.toml`
have uncommitted changes, since the database would then not be described by a
commit. Commit the sources first, then build, then commit the database and
manifest.

`build.py --write` leaves the manifest of a database it did not replace as it
is, so `built`, `versions` and `counts` keep describing the build that produced
the committed file. If only the record of the sources or of the database file is
out of date, that is updated.

## Scheme information

Each database built for CanSNPer2 2.3.0 and later holds a `scheme_info` table, which
`CanSNPer2-database --database <file> --info` prints: the scheme name, its
status from the Status row of the scheme README, the source
(`FOI-Bioinformatics/CanSNPer2-data`), the commit, the known issues listed in
`schemes.toml`, the database format, and the CanSNPer2 version and date of the
build. No `ref` is written, since the branch a database is built on goes stale
once it is merged.

`commit` is the same as `source_commit` in the manifest: the last commit that
changed the tree, SNP and references files of the scheme. It does not cover the
README or `schemes.toml`, so that an edit to the text of a README does not
require a rebuild; the status and known issues taken from them are compared by
value instead. `CanSNPer2-download` records instead the commit the requested
branch pointed to; both identify the same scheme files.

When `build.py` compares a rebuild with the committed database, the build date
and the CanSNPer2 version are left out; the other keys must match. `--write`
refuses a scheme for which no commit is found, such as outside a git checkout.

## Archive

`archive/` holds the database each rebuild replaced, named
`<scheme>_<YYMMDD>.db`. The Francisella files dated 2020 predate this
convention being automated but follow the same pattern.

The databases built before 2026-11-21 count SNP positions along the reference
sequences joined in file order. They are read by CanSNPer2 2.2 and earlier;
CanSNPer2 2.3.0 and later refuse them unless they are converted with
`cansnper2 db convert -d <file> --ref-dir <folder> -o <outdir>`.

The databases built on 2026-11-21 hold one canSNP per node, without the rank,
score and typed columns. They are read by CanSNPer2 2.3.0; CanSNPer2 2.4.0
refuses them unless they are upgraded with
`cansnper2 db convert -d <file> -o <outdir>`, which needs no
reference files, because no position moves.

## Using a database

```sh
cansnper2 ref get -d database/francisella_tularensis.db
cansnper2 type -d database/francisella_tularensis.db sample.fasta --summary
```

`cansnper2 ref get` downloads the reference genomes listed in `snp_references`.
The `genomes` table holds the genomes file of the scheme; `cansnper2 ref get
--local-references` also downloads its representative genomes, which read typing
with `--local-references` uses. CanSNPer2 2.6.0 or later is needed for these
commands; `cansnper2 scheme get <name>` builds the same database from the scheme
files.

## Changing a database directly

`cansnper2 db modify --mod-file` can change the tree of an existing database
without rebuilding it. This repository does not use it: a database changed that
way no longer matches the scheme files it was built from, and `build.py` reports
the difference. Change the scheme files and run `python build.py --write`
instead.
