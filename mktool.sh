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
  local fileListAll=out/listAll.txt
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

function cmd:jam2-generate-resodata-full {
  local fileListAll=out/listAll.txt
  if [[ ! -s $fileListAll ]]; then
    echo "$fileListAll: not found" >&2
    return 1
  fi

  {
    echo "#MASS_MEV    DEG    DEGEFF   MU  BF ANTI          KEY            NAME           PDGCODES..."
    echo "#--------    ---    ------   --  -- ----          ---            ----           -----------"
    awk '
      $2 !~ /^[0-9]+$/ && $0 != "" {

        # exclude elementary particles and other exotic objects
        if ($1 <= 110) next;

        # exclude K0_L and K0_S (duplicate with K0, K0bar)
        if ($1 == 130 || $1 == 310) next;

        # exclude massless and 100 GeV+ particles
        i=3; while (i < 10 && $i !~ /^[0-9]+$/) i++;
        spin = $i;
        BF = spin % 2 == 0 ? 1 : 2; # boson: 1, fermion: 2
        mass = $(i+3);
        if (!(0 < mass && mass < 100)) next;

        # exclude gravitino and piv/rhov
        if ($1 ~ /^10000..$/ || $1 ~ /^4900...$/) next;

        # exclude nuclei
        if (1e9 <= $1 && $1 < 2e9) next;

        if (i > 4) print "more than two particle names?" > "/dev/stderr";
        for (j = 2; j < i; j++) {
          pdg = $1;
          if (j > 2) pdg = -pdg;

          name = $j
          name = gensub(/\y0$|([*)]|\y)[-+]$/, "\\1", "g", name);
          sub(/bar0$/, "bar", name);
          if (j > 2 && name !~ /bar$/) name = name "bar";

          key = name;
          key = gensub(/'\''_([0-9]+)/, "\\1p", "g", key);
          gsub(/'\''/, "p", key);
          gsub(/_/, "", key);
          gsub(/\*/, "star", key);
          key = gensub(/\(([0-9]+)\)$/, "_\\1", "g", key);
          key = gensub(/\(([0-9]+)\)bar$/, "bar_\\1", "g", key);
          gsub(/[^a-zA-Z0-9_]/, "_", key);
          gsub(/__+/, "_", key);
          gsub(/_$/, "", key);
          key = toupper(key);

          printf("%.2f %d %g %g %d %d %s %s %d\n", mass*1e3, spin, 0.0, 0.0, BF, (j > 2 ? 1 : 0), key, name, pdg);
        }
      }
    ' "$fileListAll" | sort -n -k1,5
  } | column -t -R 1,2,3,4,5,6 > data/ResonanceJam2Full.dat
}

if declare -f "cmd:$1" &>/dev/null; then
  "cmd:$@"
else
  echo "usage: $0 subcommand..."
  return 2
fi
