# CanSNPer2-data

Data repository for CanSNPer2 https://github.com/FOI-Bioinformatics/CanSNPer2

## Repository layout

| Path | Contents |
| ---- | -------- |
| `schemes/<name>/` | One folder per canSNP scheme, holding its files as tab-separated text. This is what you edit. See `schemes/README.md` for the formats. |
| `database/` | The SQLite databases generated from `schemes/`, with a provenance manifest beside each one. |
| `database/archive/` | The database each rebuild replaced. |
| `updates/` | Tree modification files, for changing a database without rebuilding it. |
| `work_in_progress/` | Schemes under development. These are not supported and not published. |
| `schemes.toml` | Declares schemes that are not plain SNP schemes, and records known data issues. |

## Changing a scheme

Edit the text files under `schemes/<name>/`, then check and rebuild:

```sh
python validate.py          # check the scheme files
python build.py             # rebuild and compare, changing nothing
python build.py --write     # install the rebuilt databases and manifests
```

Commit the changed scheme files together with the database and manifest that
`--write` produced. Continuous integration runs the first two commands on every
pull request, so a scheme file that breaks an invariant, or a database that was
not regenerated after its sources changed, fails the build.

`validate.py` needs only Python 3.11 or later. `build.py` also needs CanSNPer2
2.2.0 or later, installed from the `master` branch; it does not depend on
flextaxd:

```sh
pip install git+https://github.com/FOI-Bioinformatics/CanSNPer2.git@master
```

progressiveMauve is not needed to build a database, only to type a sample.

The two scripts have their own tests, which CI also runs:

```sh
python -m pytest tests
```

If you want to reserve numbers for your publication (without making the snp and tree public) contact us and we will change the NEXT FREE SNP.

**NEXT FREE SNP _Francisella_**   B.317

**NEXT FREE SNP _Bacillus_**   A.Br.177 (B.Br.012, C.Br.003) Start with unified naming?

_Francisella_ - Published but excluded cansnps
|cansnp| publication   | comment                          |
| -----| ------------- | -------------------------------- |
|B.14  | Vogler2009    | Identifies same branch as B.24   |
|B.15  | Svensson2009  | Identifies same branch as B.1 and Ftind47  |
|B.17  | Svensson2009  | Identifies same branch as B.4 and Ftind38   |
|B.18  | Svensson2009  | Identifies same branch as B.11      |
|B.19  | Svensson2009  | Identifies same branch as B.12 and Ftind33      |
|B.20  | Svensson2009  | Identifies same branch as B.42     |
|B.40  | Karlsson2013  |  Identifies same branch as B.7    |
|B.41  | Karlsson2013  |  Identifies same branch as B.10     |
|B.43  | Karlsson2013  |  Identifies same branch as B.27     |
|B.69  | Sissonen2015  | Identifies same branch as B.11      |
|B.310 | Kevin2020     | Identifies same branch as B.215     |
|B.311 | Kevin2020     | Identifies same branch as B.168     |
|B.312 | Kevin2020     | Identifies same branch as B.66      |
|B.313 | Kevin2020     | Identifies same branch as B.190     |
|A.2   | Svensson2009 | |
|A.3   | Svensson2009 | |
|A.4   | Svensson2009 | |
|A.5   | Svensson2009 | |
|A.6   | Svensson2009 | |
|A.II.5| Vogler2009 | |
|A.II.7| Vogler2009 | |

