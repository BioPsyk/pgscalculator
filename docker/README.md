# Container image versioning

The Singularity/Docker image has its **own** version in `docker/VERSION`, separate from the v2 pipeline release in `VERSION` at the repo root.

| File | Meaning | When to bump |
|------|---------|----------------|
| `VERSION` | pgscalculator v2 CLI / wrapper release (e.g. `2.2.0`) | Pipeline features, shell/R scripts under `bin/` |
| `docker/VERSION` | Container base image (e.g. `0.7.0`) | Only when `docker/Dockerfile` changes (new tools, R packages, system libs) |

**Pending 0.8.0** — `docker/Dockerfile` no longer bundles Nextflow (v1 was removed; v2 is
pure bash/R/plink). The runtime user was renamed `nextflow` → `pgsuser`. This is a Dockerfile
change, so the **next** published image should be built and tagged `0.8.0` and the
`0.7.0` references in the root `README.md` updated to match. The already-published `0.7.0`
image still runs the v2 pipeline (the toolchain — R/bigsnpr, gctb, plink/plink2, PRScs — is
unchanged), so a rebuild is optional until you want the slimmer Nextflow-free image.

**0.7.0** — added multi-stage `r_builder` and LDpred2 R stack (`bigsnpr`, etc.) for the LDpred2 integration.

**0.6.0** — previous image; sufficient for sBayesR-only v2 runs without LDpred2.

Rebuild/publish the image only when `docker/VERSION` changes:

```bash
./scripts/docker-build.sh
singularity build sif/ibp-pgscalculator-base_version-$(cat docker/VERSION).sif \
  docker-daemon://ibp-pgscalculator-base:$(cat docker/VERSION)
```
