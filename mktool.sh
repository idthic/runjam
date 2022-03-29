#!/usr/bin/env bash

function cmd:update-commit-hash {
  [[ -f config.cpp ]] || return 0

  local hash=$(git show -s --format=%h 2>/dev/null)
  if [[ $hash ]]; then
    hash=+$hash
    git status -s | grep -qE '^ ?M' && hash=$hash*
  fi

  grep -qF \""$hash"\" config.cpp && return
  sed "s/package_hash = \".*\"/package_hash = \"$hash\"/" config.cpp > config.part &&
    mv config.part config.cpp
}

function cmd:jam2-check-resodata {
  local fileListAll=listAll.txt
  if [[ ! -s $fileListAll ]]; then
    echo "$fileListAll: not found" >&2
    return 1
  fi

  local resodata=$1
  local IFS=$'\n' pdg line
  for pdg in $(awk 'function abs(x) { return x < 0? -x : x; } !/^#/ {for(i=9;i<=NF;i++)print abs($i);}' "$resodata" | sort -nu); do
    line=($(awk -v pdg="$pdg" '$2 !~ /^[0-9]+$/ && $1 == pdg' "$fileListAll"))
    if ((${#line[@]} != 1)); then
      if ((${#line[@]} > 1)); then
        echo "pdg=$pdg: multiple lines"
        printf '> %s\n' "${line[@]}"
      else
        echo "pdg=$pdg: not found"
      fi
    fi
  done
}
#cmd:jam2-check-resodata data/ResonanceJam2.dat

if declare -f "cmd:$1" &>/dev/null; then
  "cmd:$@"
else
  echo "usage: $0 subcommand..."
  return 2
fi
