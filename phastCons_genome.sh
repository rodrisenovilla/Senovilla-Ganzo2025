#!/usr/bin/env bash
set -euo pipefail

# Expects one MAF per chromosome in mm10/<chrom>.maf and galGal6/<chrom>.maf,
# named to match the entries in mm10.chrom.sizes / galGal6.chrom.sizes.
# Chromosomes listed in the .chrom.sizes files that have no matching .maf are
# skipped (with a warning) rather than crashing the whole run.
# Every step is skip-if-output-exists, so a rerun after a failure resumes
# instead of redoing already-finished chromosomes.

LOG="phastCons_genome.$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOG") 2>&1

step() {
  echo "=== [$(date '+%Y-%m-%d %H:%M:%S')] START: $* ==="
}
done_step() {
  echo "=== [$(date '+%Y-%m-%d %H:%M:%S')] DONE ==="
}
skip_step() {
  echo "=== [$(date '+%Y-%m-%d %H:%M:%S')] SKIP (already done): $* ==="
}

MM10_SPECIES=rn5,dipOrd1,hetGla2,cavPor3,speTri2,oryCun2,ochPri2,hg19,panTro4,gorGor3,ponAbe2,nomLeu2,rheMac3,papHam1,calJac3,saiBol1,tarSyr1,micMur1,otoGar3,tupBel1,susScr3,vicPac1,turTru2,oviAri1,bosTau7,felCat5,canFam3,ailMel1,equCab2,myoLuc2,pteVam1,eriEur1,sorAra1,loxAfr3,proCap1,echTel1,triMan1,dasNov3,choHof1
GALGAL6_SPECIES=cotJap2,melGal5,tytAlb1,bucRhi1,anaPla1,apaVit1,calAnn1,cucCan1,chaVoc2,fulGla1,tauEry1,opiHoa1,phoRub1,colLiv1,lepDis1,merNub1,pelCri1,phaCar1,phaLep1,pteGut1,nipNip1,egrGar1,pygAde1,aptFor1,carCri1,mesUni1,eurHel1,balPav1,chlUnd1,falChe1,falPer1,aquChr2,halAlb1,halLeu1,corBra1,corCor1,acaChl1,ficAlb2,serCan1,zonAlb1,geoFor1,taeGut2,pseHum1,gavSte1,capCar1,melUnd1,amaVit1,araMac1,colStr1,picPub1,strCam1,tinGut2

# msa_view takes exactly one <infile> positional arg; combining several
# per-chromosome .ss files requires --aggregate <comma-separated-species>,
# which also re-orders/re-subsets sequences per file as needed (chromosomes
# can list their sequences in different orders). Grab the species list from
# any one .ss file already on disk rather than hardcoding it.
ss_names() {
  local ss_file="$1"
  grep '^NAMES = ' "$ss_file" | sed 's/^NAMES = //'
}

# Only keep chrom.sizes entries that have a .maf on disk AND that .maf
# actually contains at least one alignment block (a header-only .maf, with
# no "a " lines, has no alignment data and is treated the same as missing;
# maf_parse errors out on those otherwise). Preserves the size-descending
# order of the .chrom.sizes file for the ones that pass.
filter_chroms_with_maf() {
  local sizes_file="$1" maf_dir="$2"
  awk 'NR==FNR{have[$1]=1; next} ($1 in have){print $1}' \
    <(grep -l '^a ' "$maf_dir"/*.maf 2>/dev/null | xargs -n1 basename | sed 's/\.maf$//') \
    "$sizes_file"
}

MM10_ALL=$(cut -f1 mm10.chrom.sizes)
GALGAL6_ALL=$(cut -f1 galGal6.chrom.sizes)
MM10_CHROMS=$(filter_chroms_with_maf mm10.chrom.sizes mm10)
GALGAL6_CHROMS=$(filter_chroms_with_maf galGal6.chrom.sizes galGal6)

MM10_MISSING=$(comm -23 <(echo "$MM10_ALL" | sort) <(echo "$MM10_CHROMS" | sort))
GALGAL6_MISSING=$(comm -23 <(echo "$GALGAL6_ALL" | sort) <(echo "$GALGAL6_CHROMS" | sort))
if [[ -n "$MM10_MISSING" ]]; then
  echo "WARNING: mm10.chrom.sizes lists $(echo "$MM10_MISSING" | wc -l) chrom(s) with no mm10/<chrom>.maf; skipping: $(echo "$MM10_MISSING" | tr '\n' ' ')"
fi
if [[ -n "$GALGAL6_MISSING" ]]; then
  echo "WARNING: galGal6.chrom.sizes lists $(echo "$GALGAL6_MISSING" | wc -l) chrom(s) with no galGal6/<chrom>.maf; skipping: $(echo "$GALGAL6_MISSING" | tr '\n' ' ')"
fi

# 1. Prune trees once at species level (not per chromosome)
step "tree_doctor mm10 (prune Eutherians)"
if [[ -s mm10.60way.noEutherians.nh ]]; then
  skip_step "tree_doctor mm10"
else
  tree_doctor -p "$MM10_SPECIES" mm10.60way.nh > mm10.60way.noEutherians.nh
fi
done_step

step "tree_doctor galGal6 (prune birds)"
if [[ -s galGal6.77way.noBirds.nh ]]; then
  skip_step "tree_doctor galGal6"
else
  tree_doctor -p "$GALGAL6_SPECIES" galGal6.77way.nh > galGal6.77way.noBirds.nh
fi
done_step

# 2. Per-chromosome species exclusion + MAF->SS conversion
step "maf_parse + msa_view, mm10 (per chromosome)"
for chrom in $MM10_CHROMS; do
  if [[ -s mm10/${chrom}.noEutherians.ss ]]; then
    echo "--- [$(date '+%Y-%m-%d %H:%M:%S')] mm10/${chrom} SKIP (already done) ---"
    continue
  fi
  echo "--- [$(date '+%Y-%m-%d %H:%M:%S')] mm10/${chrom} ---"
  maf_parse mm10/${chrom}.maf --seqs "$MM10_SPECIES" --exclude -o MAF > mm10/${chrom}.noEutherians.maf
  msa_view --in-format MAF --out-format SS mm10/${chrom}.noEutherians.maf > mm10/${chrom}.noEutherians.ss
done
done_step

step "maf_parse + msa_view, galGal6 (per chromosome)"
for chrom in $GALGAL6_CHROMS; do
  if [[ -s galGal6/${chrom}.noBirds.ss ]]; then
    echo "--- [$(date '+%Y-%m-%d %H:%M:%S')] galGal6/${chrom} SKIP (already done) ---"
    continue
  fi
  echo "--- [$(date '+%Y-%m-%d %H:%M:%S')] galGal6/${chrom} ---"
  maf_parse galGal6/${chrom}.maf --seqs "$GALGAL6_SPECIES" --exclude -o MAF > galGal6/${chrom}.noBirds.maf
  msa_view --in-format MAF --out-format SS galGal6/${chrom}.noBirds.maf > galGal6/${chrom}.noBirds.ss
done
done_step

# 3. Concatenate all chromosomes and fit ONE genome-wide neutral model per species
step "msa_view concat mm10"
if [[ -s mm10.genome.noEutherians.ss ]]; then
  skip_step "msa_view concat mm10"
else
  MM10_SS_NAMES=$(ss_names "$(ls mm10/*.noEutherians.ss | head -1)")
  msa_view --aggregate "$MM10_SS_NAMES" --unordered-ss --out-format SS mm10/*.noEutherians.ss > mm10.genome.noEutherians.ss
fi
done_step

step "msa_view concat galGal6"
if [[ -s galGal6.genome.noBirds.ss ]]; then
  skip_step "msa_view concat galGal6"
else
  GALGAL6_SS_NAMES=$(ss_names "$(ls galGal6/*.noBirds.ss | head -1)")
  msa_view --aggregate "$GALGAL6_SS_NAMES" --unordered-ss --out-format SS galGal6/*.noBirds.ss > galGal6.genome.noBirds.ss
fi
done_step

step "phyloFit mm10 genome-wide neutral model"
if [[ -s mm10.neutral.mod ]]; then
  skip_step "phyloFit mm10"
else
  phyloFit --tree mm10.60way.noEutherians.nh --msa-format SS --out-root mm10.neutral mm10.genome.noEutherians.ss
fi
done_step

step "phyloFit galGal6 genome-wide neutral model"
if [[ -s galGal6.neutral.mod ]]; then
  skip_step "phyloFit galGal6"
else
  phyloFit --tree galGal6.77way.noBirds.nh --msa-format SS --out-root galGal6.neutral galGal6.genome.noBirds.ss
fi
done_step

# 4. Run phastCons per chromosome, reusing the shared genome-wide model
step "phastCons mm10 (per chromosome)"
for chrom in $MM10_CHROMS; do
  if [[ -s mm10/${chrom}.noEutherians.wig && -s mm10/${chrom}.noEutherians.elements.bed ]]; then
    echo "--- [$(date '+%Y-%m-%d %H:%M:%S')] mm10/${chrom} SKIP (already done) ---"
    continue
  fi
  echo "--- [$(date '+%Y-%m-%d %H:%M:%S')] mm10/${chrom} ---"
  phastCons --target-coverage 0.25 --expected-length 12 --rho 0.3 \
    --most-conserved mm10/${chrom}.noEutherians.elements.bed \
    mm10/${chrom}.noEutherians.ss mm10.neutral.mod > mm10/${chrom}.noEutherians.wig
done
done_step

step "phastCons galGal6 (per chromosome)"
for chrom in $GALGAL6_CHROMS; do
  if [[ -s galGal6/${chrom}.noBirds.wig && -s galGal6/${chrom}.noBirds.elements.bed ]]; then
    echo "--- [$(date '+%Y-%m-%d %H:%M:%S')] galGal6/${chrom} SKIP (already done) ---"
    continue
  fi
  echo "--- [$(date '+%Y-%m-%d %H:%M:%S')] galGal6/${chrom} ---"
  phastCons --target-coverage 0.25 --expected-length 12 --rho 0.3 \
    --most-conserved galGal6/${chrom}.noBirds.elements.bed \
    galGal6/${chrom}.noBirds.ss galGal6.neutral.mod > galGal6/${chrom}.noBirds.wig
done
done_step

# 5. Merge per-chromosome outputs into genome-wide files
step "merge mm10 elements.bed + wig"
cat mm10/chr*.noEutherians.elements.bed > mm10.genome.noEutherians.elements.bed
cat mm10/chr*.noEutherians.wig > mm10.genome.noEutherians.wig
done_step

step "merge galGal6 elements.bed + wig"
cat galGal6/chr*.noBirds.elements.bed > galGal6.genome.noBirds.elements.bed
cat galGal6/chr*.noBirds.wig > galGal6.genome.noBirds.wig
done_step

step "wigToBigWig mm10 genome-wide"
if [[ -s mm10.genome.noEutherians.bw ]]; then
  skip_step "wigToBigWig mm10"
else
  wigToBigWig mm10.genome.noEutherians.wig mm10.chrom.sizes mm10.genome.noEutherians.bw
fi
done_step

step "wigToBigWig galGal6 genome-wide"
if [[ -s galGal6.genome.noBirds.bw ]]; then
  skip_step "wigToBigWig galGal6"
else
  wigToBigWig galGal6.genome.noBirds.wig galGal6.chrom.sizes galGal6.genome.noBirds.bw
fi
done_step

echo "=== [$(date '+%Y-%m-%d %H:%M:%S')] ALL STEPS COMPLETED ==="
