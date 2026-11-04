# francisellaceae (work in progress)

This scheme is under development. It is not supported, it is not published with
releases, and it should not be used for typing.

| File | Contents |
| ---- | -------- |
| `francisellaceae_tree.txt` | Tree, one root-to-node path per line, tab separated. |
| `francisellaceae_genomes.txt` | Two columns, no header: a genome file name and the tree node it belongs to. |
| `francisellaceae.db` | Database built from the two files above. |
| `francisellaceae.manifest.json` | Provenance of that database. |

## Status

- CanSNPer2 cannot type against this database: it holds a taxonomy and genome
  assignments but no SNP annotations, and CanSNPer2 rejects it for lacking the
  `snp_annotation` and `snp_references` tables.
- No database has been published for this scheme, so a rebuild has no independent
  baseline. `build.py` reports a difference rather than failing on one.
- The genome assignments have been checked only for naming a node in the tree.

## Building

`validate.py` checks these files on every run. `build.py` skips the scheme unless
it is requested:

```sh
python build.py --scheme francisellaceae
```

CanSNPer2-database has no option for a genome-to-node file, so this scheme is
built with the flextaxd command. flextaxd master stops with `NameError: name 'cnp'
is not defined` when reading a CanSNPer tree, and the committed database was built
with a local, unpublished correction to `flextaxd/modules/ReadTaxonomyCanSNPer.py`.
The flextaxd version in the manifest does not reflect that correction.
