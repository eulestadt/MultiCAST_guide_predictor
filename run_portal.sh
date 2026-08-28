#!/usr/bin/env bash
# Local web portal for structure-guided MultiCAST guide design.
set -euo pipefail
cd "$(dirname "$0")"
exec streamlit run app.py --server.port "${PORT:-8501}" "$@"
