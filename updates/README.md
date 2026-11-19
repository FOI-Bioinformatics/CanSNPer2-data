# Tree modifications

This directory holds tree modification files, for changing the shape of an
existing database rather than rebuilding it from source.

A modification file lists parent and child pairs, one per line, and is applied
to a branch of a database:

```sh
CanSNPer2-database --database database/francisella_tularensis.db \
    --mod_file updates/<file>.txt --parent <node> --replace
```

For most changes this is not the route to take. Editing the files under
`source/` and running `python build.py --write` keeps the sources and the
database in step, and is what the validation in CI checks. A modification file
applied directly to a database leaves the sources describing something else.
