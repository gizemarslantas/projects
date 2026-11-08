#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
# Eski komut adı aynı güncel 10k / 1 iş parçacıklı akışa yönlenir.
exec bash "${SCRIPT_DIR}/run_rf_10k_nopc_1thread.sh" "$@"
