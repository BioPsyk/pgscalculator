#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "${SCRIPT_DIR}/../.." && pwd)"
source "${PROJECT_DIR}/bin/lib/common.sh"

VERBOSE=0
tmpdir=""
failures=0

cleanup() {
    [[ -n "$tmpdir" && -d "$tmpdir" ]] && rm -rf "$tmpdir"
}
trap cleanup EXIT

assert_eq() {
    local label="$1"
    local got="$2"
    local want="$3"
    if [[ "$got" != "$want" ]]; then
        echo "- [FAIL] ${label}: got '${got}', want '${want}'"
        failures=$((failures + 1))
    else
        echo "- [OK] ${label}"
    fi
}

assert_exit_fail() {
    local label="$1"
    shift
    if "$@" >/dev/null 2>&1; then
        echo "- [FAIL] ${label}: expected failure"
        failures=$((failures + 1))
    else
        echo "- [OK] ${label}"
    fi
}

new_config() {
    local name="$1"
    shift
    cat > "${tmpdir}/${name}.yaml" "$@"
    echo "${tmpdir}/${name}.yaml"
}

echo ">> Test methods config helpers"

tmpdir=$(mktemp -d)

cfg_default=$(new_config default <<'EOF'
outdir: /tmp/out
lddir: /tmp/ld
EOF
)
parse_config "$cfg_default"
load_active_methods "" "$cfg_default"
assert_eq "default methods" "$CFG_METHODS" "sbayesr"

cfg_inline=$(new_config inline <<'EOF'
methods: [sbayesr, ldpred2]
ldpred2:
  ld_dir: /data/ldpred2/hm3_plus
  ld_meta_file: /data/ldpred2/hm3_plus/map_hm3_plus.rds
EOF
)
parse_config "$cfg_inline"
load_active_methods "" "$cfg_inline"
assert_eq "inline methods list" "$CFG_METHODS" "sbayesr ldpred2"
has_method sbayesr "$CFG_METHODS"
assert_eq "has_method sbayesr" "$?" "0"
has_method ldpred2 "$CFG_METHODS"
assert_eq "has_method ldpred2" "$?" "0"

cfg_block=$(new_config block <<'EOF'
methods:
  - ldpred2
ldpred2:
  ld_dir: /data/ldpred2/hm3_plus
EOF
)
parse_config "$cfg_block"
load_active_methods "" "$cfg_block"
assert_eq "block methods list" "$CFG_METHODS" "ldpred2"
assert_eq "resolved ld_meta default" "${CFG_LDPRED2_LD_META_FILE_RESOLVED:-}" "/data/ldpred2/hm3_plus/map_hm3_plus.rds"

cfg_cli=$(new_config cli <<'EOF'
methods: [ldpred2]
ldpred2:
  ld_dir: /data/ldpred2/hm3_plus
EOF
)
parse_config "$cfg_cli"
load_active_methods "sbayesr" "$cfg_cli"
assert_eq "CLI override wins" "$CFG_METHODS" "sbayesr"

cfg_missing_ld=$(new_config missing_ld <<'EOF'
methods: [ldpred2]
EOF
)
unset CFG_LDPRED2_LD_DIR CFG_LDPRED2_LD_META_FILE CFG_LDPRED2_LD_META_FILE_RESOLVED CFG_METHODS
parse_config "$cfg_missing_ld"
if ( load_active_methods "" "$cfg_missing_ld" ) >/dev/null 2>&1; then
    echo "- [FAIL] ldpred2 without ld_dir fails: expected failure"
    failures=$((failures + 1))
else
    echo "- [OK] ldpred2 without ld_dir fails"
fi

cfg_slurm=$(new_config slurm <<'EOF'
slurm:
  score: { mem: 8g, cpus: 2, time: '1:00:00', max_parallel: 11 }
EOF
)
assert_eq "score_sbayesr inherits score" "$(resolve_slurm_step_settings score_sbayesr "$cfg_slurm")" "mem: 8g, cpus: 2, time: '1:00:00', max_parallel: 11"
assert_eq "score_ldpred2 inherits score" "$(resolve_slurm_step_settings score_ldpred2 "$cfg_slurm")" "mem: 8g, cpus: 2, time: '1:00:00', max_parallel: 11"

cfg_slurm_override=$(new_config slurm_override <<'EOF'
slurm:
  score: { mem: 8g, cpus: 2, time: '1:00:00', max_parallel: 11 }
  score_ldpred2: { mem: 12g, cpus: 8, time: '2:00:00', max_parallel: 22 }
EOF
)
assert_eq "score_ldpred2 explicit profile" "$(resolve_slurm_step_settings score_ldpred2 "$cfg_slurm_override")" "mem: 12g, cpus: 8, time: '2:00:00', max_parallel: 22"

if [[ "$failures" -gt 0 ]]; then
    echo "FAILED: ${failures} assertion(s)"
    exit 1
fi

echo "All methods config tests passed."
