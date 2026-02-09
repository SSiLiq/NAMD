#!/bin/bash

# ==============================================================================
# MASTER SCRIPT: NAMD SETUP (VERSÃO FIEL AO EM.CONF + RESTRAINTS)
# ==============================================================================

clear
echo "========================================================"
echo "      NAMD SYSTEM"
echo "========================================================"

# --- 1. ENTRADA DE DADOS ---
read -p "Archive name - PDB (peptide): " PEPT_PDB_NAME
read -p "Archive name - PSF (peptide): " PEPT_PSF_NAME

PEPT_PDB="peptide/$PEPT_PDB_NAME"
PEPT_PSF="peptide/$PEPT_PSF_NAME"

# --- 2. CONFIGURAÇÃO GERAL ---
TOPPAR_DIR="charmm-gui/toppar"
MEMB_PDB="membrane/membrane.pdb"
OUT_NAME="system_final"
PROJ_NAME="system_$(date +%Y%m%d_%H%M)"
TEMP=298

# ==============================================================================
# STEP 1: BUILD SYSTEM & RESTRAINTS & BOX CALC (VMD)
# ==============================================================================
echo ">>> [1/2] Building system, restraints and calculating box..."

cat << EOF > build_system.tcl
package require psfgen
package require solvate
package require autoionize

resetpsf
topology $TOPPAR_DIR/top_all36_prot.rtf
topology $TOPPAR_DIR/top_all36_lipid.rtf
topology $TOPPAR_DIR/top_all36_carb.rtf
topology $TOPPAR_DIR/top_all36_na.rtf
topology $TOPPAR_DIR/top_all36_cgenff.rtf

pdbalias residue HIS HSD
pdbalias atom ILE CD1 CD
pdbalias atom HOH O OH2
pdbalias residue HOH TIP3

segment MEMB { pdb $MEMB_PDB }
coordpdb $MEMB_PDB MEMB
segment PROT { pdb $PEPT_PDB }
coordpdb $PEPT_PDB PROT

guesscoord
writepsf temp_dry.psf
writepdb temp_dry.pdb

solvate temp_dry.psf temp_dry.pdb -minmax {{-35 -35 -60} {35 35 60}} -o temp_solv
autoionize -psf temp_solv.psf -pdb temp_solv.pdb -sc 0.15 -cation SOD -anion CLA -o temp_ion

mol new temp_ion.psf; mol addfile temp_ion.pdb
set all [atomselect top all]

# Centralização
set center_all [measure center \$all]
\$all moveby [vecscale -1.0 \$center_all]

# Limpeza de águas na membrana
set bad_water [atomselect top "water and z > -14 and z < 14"]
if { [\$bad_water num] > 0 } {
    set bad_res [lsort -unique [\$bad_water get residue]]
    set keep [atomselect top "not residue \$bad_res"]
    \$keep writepsf "$OUT_NAME.psf"
    \$keep writepdb "$OUT_NAME.pdb"
} else {
    \$all writepsf "$OUT_NAME.psf"
    \$all writepdb "$OUT_NAME.pdb"
}

# --- CÁLCULO DE CAIXA (MinMax) ---
mol delete top; mol new "$OUT_NAME.psf"; mol addfile "$OUT_NAME.pdb"
set final_sel [atomselect top all]
set minmax [measure minmax \$final_sel]
set final_center [measure center \$final_sel]
set vec [vecsub [lindex \$minmax 1] [lindex \$minmax 0]]

set fp [open "box_params.txt" w]
puts \$fp "BOX_X=[lindex \$vec 0]"
puts \$fp "BOX_Y=[lindex \$vec 1]"
puts \$fp "BOX_Z=[lindex \$vec 2]"
puts \$fp "CEN_X=[lindex \$final_center 0]"
puts \$fp "CEN_Y=[lindex \$final_center 1]"
puts \$fp "CEN_Z=[lindex \$final_center 2]"
close \$fp

# --- CRIAÇÃO DO ARQUIVO DE RESTRIÇÃO ---
\$final_sel set beta 0.0
set bb [atomselect top "protein and backbone"]
if {[\$bb num] > 0} { \$bb set beta 1.0 }
\$final_sel writepdb restraints.pdb

exit
EOF

vmd -dispdev text -e build_system.tcl > vmd_build.log 2>&1
source box_params.txt

# ==============================================================================
# STEP 2: GENERATING em_final.conf (CÓPIA FIEL DO SEU ORIGINAL)
# ==============================================================================
echo ">>> [2/2] Generating em_final.conf from your original template..."

cat << EOF > em.conf
# Configuration file for NAMD simulations

set restart_job 0
set do_pbc 1
set do_pme 1
set do_lgv_pressure 1
set do_freeze 0
set do_restraints 0
set do_imd 0
set do_minimize 1
set do_md 0
set do_abf 0
set do_metadyn 0
set do_ramd 0
set temperature $TEMP
set pressure 1.01325
set basename $OUT_NAME
set inputname $OUT_NAME
set outputname em_output
set em_steps 15000
set md_steps 5000000

structure          \$basename.psf
coordinates        \$basename.pdb
outputName         \$outputname
binaryoutput       yes

if {\$restart_job} {
    binCoordinates     \$inputname.coor
    binVelocities      \$inputname.vel
    extendedSystem     \$inputname.xsc
} else {
    temperature        \$temperature
}

firsttimestep      0

paraTypeCharmm     on
parameters         ../$TOPPAR_DIR/par_all36m_prot.prm
parameters         ../$TOPPAR_DIR/par_all36_lipid.prm
parameters         ../$TOPPAR_DIR/par_all36_carb.prm
parameters         ../$TOPPAR_DIR/par_all36_na.prm
parameters         ../$TOPPAR_DIR/par_all36_cgenff.prm
parameters         ../$TOPPAR_DIR/toppar_water_ions.str


if {\$do_pbc} {
cellBasisVector1   $BOX_X 0 0
cellBasisVector2   0 $BOX_Y 0
cellBasisVector3   0 0 $BOX_Z
cellOrigin         $CEN_X $CEN_Y $CEN_Z
wrapWater          on
wrapAll            on
}

exclude            scaled1-4
1-4scaling         1.0
cutoff             12.0
switching          on
switchdist         10.0
pairlistdist       14.0
timestep           2.0
rigidBonds         all
nonbondedFreq      1
fullElectFrequency 2  
stepspercycle      20
pairlistsPerCycle  2 

if {\$do_pme} {
PME                yes
PMEGridSpacing     1.0
}

langevin           on
langevinDamping    1
langevinTemp       \$temperature
langevinHydrogen   off

if {\$do_lgv_pressure} {
useGroupPressure      yes 
useFlexibleCell       yes 
useConstantArea       yes 

langevinPiston        on
langevinPistonTarget  \$pressure
langevinPistonPeriod  200.0
langevinPistonDecay   100.0
langevinPistonTemp    \$temperature
}

restartfreq        10000
dcdfreq            10000
xstFreq            10000
outputEnergies     1000
outputPressure     1000

# Minimization
if {\$do_minimize} {
minimize            \$em_steps
reinitvels          \$temperature
}

# MD
if {\$do_md} {
run           \$md_steps
}
EOF

# Organização final
mkdir -p "$PROJ_NAME"
mv "$OUT_NAME.psf" "$OUT_NAME.pdb" "restraints.pdb" "em_final.conf" "$PROJ_NAME/"
rm temp_* build_system.tcl box_params.txt vmd_build.log

echo "========================================================"
echo "  📦 Pasta: $PROJ_NAME"
echo "  ⚙️  Restraints.pdb gerado (backbone beta=1.0)"
echo "========================================================"
