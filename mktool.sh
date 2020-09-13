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

if declare -f "cmd:$1" &>/dev/null; then
  "cmd:$@"
fi
