#!/usr/bin/env bash
set -eo pipefail

if ! command -v mamba >/dev/null 2>&1; then
  if ! command -v module >/dev/null 2>&1; then
    for modules_init in /etc/profile.d/modules.sh /usr/share/lmod/lmod/init/bash; do
      if [[ -r "${modules_init}" ]]; then
        # shellcheck source=/dev/null
        source "${modules_init}"
        break
      fi
    done
  fi
  module load Miniforge3
fi

set -u
exec mamba run -n cice_setup python "$@"
