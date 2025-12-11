#!/bin/bash
# Test script to verify mounted directories are accessible inside container
# Run this inside an interactive container session

echo "=== Testing pgscalculator v2.0.0 Mount Points ==="
echo ""

# Check if we're in container
if [[ -d "/pgscalculator" ]]; then
  echo "[OK] Container detected"
else
  echo "[WARN] Not in container - some paths may not exist"
fi

echo ""
echo "=== Checking CLI ==="
if command -v pgscalculator &> /dev/null; then
  echo "[OK] pgscalculator CLI found"
  pgscalculator --version
else
  echo "[ERROR] pgscalculator CLI not found in PATH"
fi

echo ""
echo "=== Checking Mount Points ==="

# Check standard mount points
mount_points=(
  "/pgscalculator/input"
  "/pgscalculator/outdir"
  "/pgscalculator/genodir"
  "/pgscalculator/genodir2"
  "/pgscalculator/confdir"
  "/pgscalculator/work"
  "/tmp"
)

for mp in "${mount_points[@]}"; do
  if [[ -d "$mp" ]] || [[ -f "$mp" ]]; then
    echo "[OK] $mp exists"
    if [[ -d "$mp" ]]; then
      count=$(find "$mp" -maxdepth 1 | wc -l)
      echo "      Contains $count items"
    fi
  else
    echo "[INFO] $mp not mounted (this is OK if not used)"
  fi
done

# Check LD directory (name varies)
echo ""
echo "=== Checking LD Reference ==="
ld_dirs=$(find /pgscalculator -maxdepth 1 -type d -name "*ld*" -o -name "*band*" 2>/dev/null)
if [[ -n "$ld_dirs" ]]; then
  echo "[OK] LD directories found:"
  echo "$ld_dirs" | while read dir; do
    echo "      $dir"
    if [[ -d "$dir" ]]; then
      count=$(find "$dir" -maxdepth 1 -name "*.info" -o -name "*.bin" 2>/dev/null | wc -l)
      echo "        Contains $count LD files"
    fi
  done
else
  echo "[INFO] No LD directories found (will be mounted at runtime)"
fi

echo ""
echo "=== Testing Config File Access ==="
if [[ -f "/pgscalculator/outdir/config.yaml" ]]; then
  echo "[OK] Config file found"
  echo "First few lines:"
  head -10 /pgscalculator/outdir/config.yaml
else
  echo "[INFO] Config file not found (will be created at runtime)"
fi

echo ""
echo "=== Testing CLI Help ==="
pgscalculator --help 2>&1 | head -20

echo ""
echo "=== Testing Status Command ==="
if [[ -d "/pgscalculator/outdir" ]]; then
  pgscalculator status --config /pgscalculator/outdir/config.yaml 2>&1 || echo "Status check completed"
fi

echo ""
echo "=== Test Complete ==="


