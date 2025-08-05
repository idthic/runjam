#!/usr/bin/env bash

function cmd:update-commit-hash {
  local file=${1-config.cpp}
  [[ -f $file ]] || return 0

  local hash=$(git show -s --format=%h 2>/dev/null)
  if [[ $hash ]]; then
    hash=+$hash
    git status -s | grep -qE '^ ?M' && hash=$hash*
  fi

  grep -qF \""$hash"\" "$file" && return
  sed "s/package_hash = \".*\"/package_hash = \"$hash\"/" "$file" > "$file.part" &&
    mv "$file.part" "$file"
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
  local fileListAll=${1-out/listAll.txt}
  if [[ ! -s $fileListAll ]]; then
    echo "$fileListAll: not found" >&2
    return 1
  fi

  local -a awk2_args=()
  if [[ -s feeddown.txt ]]; then
    awk2_args=(mode=feeddown feeddown.txt)
  fi
  awk2_args+=(mode=input -)

  {
    echo "#MASS_MEV    DEG    DEGEFF   MU  BF ANTI          KEY            NAME           PDGCODES..."
    echo "#--------    ---    ------   --  -- ----          ---            ----           -----------"
    awk '
      $2 !~ /^[0-9]+$/ && $0 != "" {

        # exclude no number
        if ($1 !~ /[0-9]/) next;

        # exclude elementary particles and other exotic objects
        if ($1 <= 110) next;

        # exclude K0_L and K0_S (duplicate with K0, K0bar)
        if ($1 == 130 || $1 == 310) next;

        # exclude massless and 100 GeV+ particles
        i=3; while (i < 10 && $i !~ /^[0-9]+$/) i++;
        spin = $i;
        BF = spin % 2 == 1 ? 1 : 2; # boson: 1, fermion: 2
        mass = $(i+3);
        if (!(0 < mass && mass < 100)) next;

        # exclude gravitino and piv/rhov
        if ($1 ~ /^10000..$/ || $1 ~ /^4900...$/) next;

        # exclude diquarks
        if ($1 ~ /^-?[1-9][0-9]0[0-9]$/) next;

        # exclude nuclei
        if (1e9 <= $1 && $1 < 2e9) next;

        # exclude qv1qv1
        if ($1 ~ /^-?4901103$/) next;

        # exclude Deltav1
        if ($1 ~ /^-?490[1-8]114$/) next;

        if (i > 4) print "more than two particle names?" > "/dev/stderr";
        for (j = 2; j < i; j++) {
          pdg = $1;
          if (j > 2) pdg = -pdg;
          anti = BF == 2 && j > 2 ? 1 : 0; # Note: anti-baryon

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

          printf("%.2f %d %g %g %d %d %s %s %d\n", mass*1e3, spin, 0.0, 0.0, BF, anti, key, name, pdg);
        }
      }
    ' "$fileListAll" | sort -n -k1 -k2 -k5 | awk '
      mode == "feeddown" {
        name = $1;
        npi = $6; # charged pion numbers
        sub(/\(.*\)/, "", npi);
        npi = npi + 0.0;
        g_degeff[name] = npi;
        next;
      }

      function flush_line() {
        if (pnam == "") return;

        # Get degeff
        degeff = 0.0;
        if (g_degeff[pkey] != "") degeff = g_degeff[pkey];
        degeff = sprintf("%.3f", degeff);

        print mass, deg, degeff, 0, bf, anti, pkey, pnam, pdg;
        pnam = "";
      }

      {
        key = $1 "," $2 "," $5 "," $6;
        if (key != okey) {
          flush_line();
          okey = key;

          mass = $1;
          bf = $5;
          anti = $6;
          pkey = $7;
          pnam = $8;

          deg = $2;
          pdg = $9;
        } else {
          deg += $2;
          pdg = pdg " " $9
        }
      }
      END { flush_line(); }
    ' "${awk2_args[@]}"
  } | column -t -R 1,2,3,4,5,6,9,10,11,12,13,14 | sed 's/[[:space:]]\{1,\}$//' > data/ResonanceJam2Full.dat
}

if declare -f "cmd:$1" &>/dev/null; then
  "cmd:$@"
else
  echo "usage: $0 subcommand..."
  return 2
fi
