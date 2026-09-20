#!/bin/bash
# Add a tabulated EAM potential to critic2's catalogue (dat/eam), so it can be
# found automatically instead of having to be named on every EDIT RELAX EAM or
# TRICK ENERGY EAM line.
#
# critic2 ships potentials for Ag Al Au Co Cr Cu Fe Mg Mo Ni Pb Pd Pt Ta Ti V W
# Zr (see dat/eam/README). Everything else -- the alkali and alkaline-earth
# metals, Sc Mn Zn Y Nb Tc Ru Rh Cd Hf Re Os Ir Hg, Ga In Sn Tl Bi, and the
# lanthanides and actinides -- has to come from elsewhere, because no source
# with a licence compatible with critic2's GPL-3 provides them. This script is
# how you add one.
#
# The best source is the NIST Interatomic Potentials Repository:
#
#     https://www.ctcms.nist.gov/potentials/
#
# Browse to the element you want, pick an entry, and copy the download link of
# its "eam.alloy" or "eam.fs" file. NIST states no licence for the files it
# hosts, only that users are "encouraged to download and use interatomic
# potentials, with proper acknowledgement" -- so use them, cite the paper the
# entry names, and do NOT redistribute them with critic2 or anything else.
# OpenKIM (https://openkim.org) also hosts EAM models, but most of its
# collection is CDDL-1.0; check the LICENSE inside an item before relying on it.
#
# Only the setfl (.eam.alloy) and eam/fs (.eam.fs) formats are supported.
# Formats that start like setfl but carry extra tables (.adp, .cdeam) and the
# single-element DYNAMO funcfl layout (.eam) are rejected by critic2.
#
# Usage:
#   tools/add-eam-potential.sh URL_OR_FILE [...]
#
# Each argument is either a URL to download or a local file to copy. The
# potential is placed in the catalogue directory, its element list is read from
# its own header, and a line is appended to the index. Examples:
#
#   tools/add-eam-potential.sh https://www.ctcms.nist.gov/potentials/Download/.../Nb.eam.fs
#   tools/add-eam-potential.sh ~/Downloads/Hf.eam.alloy ~/Downloads/Zn.eam.alloy
#
# Environment overrides:
#   EAMDIR   catalogue directory (default: the dat/eam of this source tree, or
#            $CRITIC_HOME/eam if that is set and the source tree is not found)

set -e

usage() {
   sed -n '2,/^set -e/p' "$0" | sed 's/^# \{0,1\}//; /^set -e/d'
   exit "${1:-0}"
}

[ $# -eq 0 ] && usage 1
case "$1" in -h|--help|-help) usage 0 ;; esac

## locate the catalogue directory
if [ -n "$EAMDIR" ]; then
   dir="$EAMDIR"
else
   here=$(cd "$(dirname "$0")/.." && pwd)
   if [ -d "$here/dat/eam" ]; then
      dir="$here/dat/eam"
   elif [ -n "$CRITIC_HOME" ] && [ -d "$CRITIC_HOME/eam" ]; then
      dir="$CRITIC_HOME/eam"
   else
      echo "error: cannot find the catalogue directory; set EAMDIR" >&2
      exit 1
   fi
fi
index="$dir/index"
[ -d "$dir" ] || { echo "error: $dir is not a directory" >&2; exit 1; }

## fetch or copy, then validate and register
for src in "$@"; do
   name=$(basename "$src")
   dest="$dir/$name"

   case "$name" in
      *.eam.alloy|*.eam.fs) ;;
      *) echo "error: $name is not a .eam.alloy or .eam.fs file (only setfl and"\
              "eam/fs are supported)" >&2; exit 1 ;;
   esac
   if [ -e "$dest" ]; then
      echo "error: $dest already exists; remove it first" >&2; exit 1
   fi

   case "$src" in
      http://*|https://*|ftp://*)
         echo "fetching $name"
         if command -v curl > /dev/null; then
            curl -fsSL -o "$dest" "$src"
         elif command -v wget > /dev/null; then
            wget -q -O "$dest" "$src"
         else
            echo "error: neither curl nor wget is available" >&2; exit 1
         fi
         ;;
      *)
         [ -f "$src" ] || { echo "error: no such file: $src" >&2; exit 1; }
         cp "$src" "$dest"
         ;;
   esac

   ## the element line is line 4: a count followed by that many symbols
   nel=$(sed -n '4p' "$dest" | awk '{print $1}')
   els=$(sed -n '4p' "$dest" | awk '{$1=""; print}' | xargs)
   nsym=$(echo "$els" | wc -w)
   if ! [ "$nel" -eq "$nel" ] 2> /dev/null || [ "$nel" -lt 1 ] || [ "$nsym" -ne "$nel" ]; then
      rm -f "$dest"
      echo "error: $name does not look like a setfl/eam-fs file: line 4 should be"\
           "a count followed by that many element symbols" >&2
      exit 1
   fi

   ## have critic2 read it, so a bad file is caught here and not at run time
   if command -v critic2 > /dev/null || [ -x "$here/build/src/critic2" ]; then
      c2=$(command -v critic2 || echo "$here/build/src/critic2")
      tmp=$(mktemp -d)
      first=$(echo "$els" | awk '{print $1}')
      printf 'crystal\n cell 3.6 3.6 3.6 90 90 90 angstrom\n spg 1\n neq 0 0 0 %s\nendcrystal\ntrick energy eam %s\n' \
         "$first" "$dest" > "$tmp/t.cri"
      if ! "$c2" "$tmp/t.cri" 2>&1 | grep -q "energy at input geometry"; then
         rm -rf "$tmp"; rm -f "$dest"
         echo "error: critic2 could not read $name; not added" >&2
         exit 1
      fi
      rm -rf "$tmp"
   else
      echo "  (critic2 not found, skipping the read-back check)"
   fi

   ## register it, before the general-purpose entries so it wins for its elements
   if grep -q "^[[:space:]]*$name[[:space:]]" "$index" 2> /dev/null; then
      echo "  $name is already in the index, not adding a second line"
   else
      printf '%-28s %s\n' "$name" "$els" >> "$index"
      echo "  added: $name  ($els)"
      echo "  NOTE: appended at the end of $index. The first entry covering a"
      echo "        system wins, so move the line up if this potential should"
      echo "        take precedence for its elements."
   fi
done
