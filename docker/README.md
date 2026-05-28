# Container image versioning

The Singularity/Docker image has its **own** version in `docker/VERSION`, separate from the v2 pipeline release in `VERSION.v2` at the repo root.

| File | Meaning | When to bump |
|------|---------|----------------|
| `VERSION.v2` | pgscalculator v2 CLI / wrapper release (e.g. `2.2.0`) | Pipeline features, shell/R scripts under `bin/` |
| `docker/VERSION` | Container base image (e.g. `0.7.0`) | Only when `docker/Dockerfile` changes (new tools, R packages, system libs) |

**0.7.0** — added multi-stage `r_builder` and LDpred2 R stack (`bigsnpr`, etc.) for the LDpred2 integration.

**0.6.0** — previous image; sufficient for sBayesR-only v2 runs without LDpred2.

Rebuild/publish the image only when `docker/VERSION` changes:

```bash
./scripts/docker-build.sh
singularity build sif/ibp-pgscalculator-base_version-$(cat docker/VERSION).sif \
  docker-daemon://ibp-pgscalculator-base:$(cat docker/VERSION)
```
