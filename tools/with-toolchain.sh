#!/usr/bin/env bash
# Build/test helper for this machine.
#
# R's default C++ compiler here is anaconda clang 14, which cannot parse the
# macOS SDK libc++ headers ("no template named 'is_integral'").  Apple clang 17
# at /usr/bin/g++ does work.  This wrapper puts it first on PATH and provides a
# project-local Makevars override, then runs whatever command you pass.
#
#   tools/with-toolchain.sh Rscript -e 'devtools::test()'
#
# A conda fallback is documented in README.md if this ever stops working.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export PATH="/usr/bin:${PATH}"

if [[ ! -f "$REPO_ROOT/tools/Makevars.local" ]]; then
  cat > "$REPO_ROOT/tools/Makevars.local" <<'EOF'
CC=/usr/bin/g++ -arch arm64
CXX=/usr/bin/g++ -arch arm64 -std=gnu++17
SHLIB_CXXLD=/usr/bin/g++ -arch arm64
CXX17=/usr/bin/g++ -arch arm64
EOF
fi
export R_MAKEVARS_USER="${R_MAKEVARS_USER:-$REPO_ROOT/tools/Makevars.local}"

exec "$@"
