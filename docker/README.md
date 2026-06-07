# Container image versioning

The Singularity/Docker image has its **own** version in `docker/VERSION`, separate from the v2 pipeline release in `VERSION` at the repo root.

| File | Meaning | When to bump |
|------|---------|----------------|
| `VERSION` | pgscalculator v2 CLI / wrapper release (e.g. `2.2.0`) | Pipeline features, shell/R scripts under `bin/` |
| `docker/VERSION` | Container base image (e.g. `0.8.0`) | Only when `docker/Dockerfile` changes (new tools, R packages, system libs) |

**0.8.0** — `docker/Dockerfile` no longer bundles Nextflow (v1 was retired; v2 is pure
bash/R/plink); the runtime user was renamed `nextflow` → `pgsuser`. Functionally equivalent
to `0.7.0` for running v2 (same toolchain: R/bigsnpr, gctb, plink/plink2, PRScs), just
slimmer. Build + push with `scripts/docker-deploy-build.sh` / `docker-deploy-push.sh`.

**0.7.0** — added multi-stage `r_builder` and LDpred2 R stack (`bigsnpr`, etc.) for the LDpred2 integration.

**0.6.0** — previous image; sufficient for sBayesR-only v2 runs without LDpred2.

Rebuild/publish the image only when `docker/VERSION` changes:

```bash
./scripts/docker-build.sh
singularity build sif/ibp-pgscalculator-base_version-$(cat docker/VERSION).sif \
  docker-daemon://ibp-pgscalculator-base:$(cat docker/VERSION)
```
