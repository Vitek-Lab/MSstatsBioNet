# Network backend classes

A backend is a source of prior-knowledge networks, such as INDRA.
[`get_network()`](https://vitek-lab.github.io/MSstatsBioNet/reference/get_network.md)
dispatches on the backend and the query, so each (backend, query) pair
has its own method. `NetworkBackend` is virtual: create a backend with a
constructor such as
[`indra_backend()`](https://vitek-lab.github.io/MSstatsBioNet/reference/indra_backend.md).
Other packages can add a backend by extending `NetworkBackend` and
writing methods for the generics.

## Details

`IndraBackend` queries INDRA CoGEx for networks and grounds names with
Gilda, INDRA's grounding service. Curations come from the INDRA
database.

## Slots

- `cogex_url`:

  base URL of INDRA CoGEx

- `grounding_url`:

  base URL of Gilda

- `curation_url`:

  base URL of the INDRA database, which holds the curations

## See also

[`indra_backend()`](https://vitek-lab.github.io/MSstatsBioNet/reference/indra_backend.md),
[`backend_capabilities()`](https://vitek-lab.github.io/MSstatsBioNet/reference/backend_capabilities.md)
