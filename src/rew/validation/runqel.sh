#!/bin/bash

##### SETTINGS #####

ENERGY_MIN=500
ENERGY_MAX=1500
ENERGY_STEP=500

MA_MIN=800
MA_MAX=1200
MA_STEP=200

OUTDIR="qelma/"
mkdir -p $OUTDIR

##### SIMULATIONS ######

# for energy in $(seq $ENERGY_MIN $ENERGY_STEP $ENERGY_MAX)
# do
#     for ma in $(seq $MA_MIN $MA_STEP $MA_MAX)
#     do
#         ../../../bin/nuwro -i "parqel.txt" \
#             -p "beam_particle = 14" \
#             -p "beam_energy = $energy" \
#             -p "qel_cc_axial_mass = $ma" \
#             -p "qel_nc_axial_mass = $ma" \
#             -p "dyn_qel_cc = 1" \
#             -p "dyn_qel_nc = 0" \
#             -o "${OUTDIR}/numu${energy}_cc_ma${ma}.root" &

#         ../../../bin/nuwro -i "parqel.txt" \
#             -p "beam_particle = 14" \
#             -p "beam_energy = $energy" \
#             -p "qel_cc_axial_mass = $ma" \
#             -p "qel_nc_axial_mass = $ma" \
#             -p "dyn_qel_cc = 0" \
#             -p "dyn_qel_nc = 1" \
#             -o "${OUTDIR}/numu${energy}_nc_ma${ma}.root" &

#         ../../../bin/nuwro -i "parqel.txt" \
#             -p "beam_particle = -14" \
#             -p "beam_energy = $energy" \
#             -p "qel_cc_axial_mass = $ma" \
#             -p "qel_nc_axial_mass = $ma" \
#             -p "dyn_qel_cc = 1" \
#             -p "dyn_qel_nc = 0" \
#             -o "${OUTDIR}/numubar${energy}_cc_ma${ma}.root" &

#         ../../../bin/nuwro -i "parqel.txt" \
#             -p "beam_particle = -14" \
#             -p "beam_energy = $energy" \
#             -p "qel_cc_axial_mass = $ma" \
#             -p "qel_nc_axial_mass = $ma" \
#             -p "dyn_qel_cc = 0" \
#             -p "dyn_qel_nc = 1" \
#             -o "${OUTDIR}/numubar${energy}_nc_ma${ma}.root"
#     done
# done

##### REWEIGHTING #####

# for energy in $(seq $ENERGY_MIN $ENERGY_STEP $ENERGY_MAX)
# do
#     for ma in $(seq $MA_MIN $MA_STEP $MA_MAX)
#     do
#         ../../../bin/reweight_to \
#             "${OUTDIR}/numu${energy}_cc_ma800.root" \
#             -o "${OUTDIR}/numu${energy}_cc_rewto_ma${ma}.root" \
#             -p qel_cc_axial_mass $ma \
#             -p qel_nc_axial_mass $ma &

#         ../../../bin/reweight_to \
#             "${OUTDIR}/numu${energy}_nc_ma800.root" \
#             -o "${OUTDIR}/numu${energy}_nc_rewto_ma${ma}.root" \
#             -p qel_cc_axial_mass $ma \
#             -p qel_nc_axial_mass $ma &

#         ../../../bin/reweight_to \
#             "${OUTDIR}/numubar${energy}_cc_ma800.root" \
#             -o "${OUTDIR}/numubar${energy}_cc_rewto_ma${ma}.root" \
#             -p qel_cc_axial_mass $ma \
#             -p qel_nc_axial_mass $ma &

#         ../../../bin/reweight_to \
#             "${OUTDIR}/numubar${energy}_nc_ma800.root" \
#             -o "${OUTDIR}/numubar${energy}_nc_rewto_ma${ma}.root" \
#             -p qel_cc_axial_mass $ma \
#             -p qel_nc_axial_mass $ma
#     done
# done

###### RUN MACRO ######

# for energy in $(seq $ENERGY_MIN $ENERGY_STEP $ENERGY_MAX)
# do
#     for ma in $(seq $MA_MIN $MA_STEP $MA_MAX)
#     do
#         ../../../bin/myroot -b -q 'qel.c("'${OUTDIR}'/numu'${energy}'_cc_ma'${ma}'")' &
#         ../../../bin/myroot -b -q 'qel.c("'${OUTDIR}'/numu'${energy}'_nc_ma'${ma}'")' &
#         ../../../bin/myroot -b -q 'qel.c("'${OUTDIR}'/numu'${energy}'_cc_rewto_ma'${ma}'")' &
#         ../../../bin/myroot -b -q 'qel.c("'${OUTDIR}'/numu'${energy}'_nc_rewto_ma'${ma}'")'
#         ../../../bin/myroot -b -q 'qel.c("'${OUTDIR}'/numubar'${energy}'_cc_ma'${ma}'")' &
#         ../../../bin/myroot -b -q 'qel.c("'${OUTDIR}'/numubar'${energy}'_nc_ma'${ma}'")' &
#         ../../../bin/myroot -b -q 'qel.c("'${OUTDIR}'/numubar'${energy}'_cc_rewto_ma'${ma}'")' &
#         ../../../bin/myroot -b -q 'qel.c("'${OUTDIR}'/numubar'${energy}'_nc_rewto_ma'${ma}'")'
#     done
# done

##### COMBO TEST #####

# ../../../bin/nuwro -i "parqel.txt" \
#     -p "beam_particle = 14" \
#     -p "beam_energy = 1000" \
#     -p "qel_cc_axial_mass = 1000" \
#     -p "qel_nc_axial_mass = 1000" \
#     -p "qel_s_axial_mass = 1000" \
#     -p "delta_s = 0" \
#     -p "dyn_qel_cc = 1" \
#     -p "dyn_qel_nc = 1" \
#     -o "${OUTDIR}/numu1000.root"

# ../../../bin/nuwro -i "parqel.txt" \
#     -p "beam_particle = 14" \
#     -p "beam_energy = 1000" \
#     -p "qel_cc_axial_mass = 1200" \
#     -p "qel_nc_axial_mass = 1200" \
#     -p "qel_s_axial_mass = 1200" \
#     -p "delta_s = -0.2" \
#     -p "dyn_qel_cc = 1" \
#     -p "dyn_qel_nc = 1" \
#     -o "${OUTDIR}/numu1200.root"

# ../../../bin/reweight_to \
#     "${OUTDIR}/numu1200.root" \
#     -o "${OUTDIR}/numu1000r.root" \
#     -p qel_cc_axial_mass 1000 \
#     -p qel_nc_axial_mass 1000 \
#     -p qel_s_axial_mass 1000 \
#     -p delta_s 0

# ../../../bin/reweight_to \
#     "${OUTDIR}/numu1000.root" \
#     -o "${OUTDIR}/numu1200r.root" \
#     -p qel_cc_axial_mass 1200 \
#     -p qel_nc_axial_mass 1200 \
#     -p qel_s_axial_mass 1200 \
#     -p delta_s -0.2

# ../../../bin/myroot -b -q 'qel.c("'${OUTDIR}'/numu1000")' &
# ../../../bin/myroot -b -q 'qel.c("'${OUTDIR}'/numu1000r")' &
# ../../../bin/myroot -b -q 'qel.c("'${OUTDIR}'/numu1200")' &
# ../../../bin/myroot -b -q 'qel.c("'${OUTDIR}'/numu1200r")' &

##### STRANGE TEST #####

../../../bin/nuwro -i "parqel.txt" \
    -p "dyn_qel_cc = 0" \
    -p "dyn_qel_nc = 1" \
    -p "delta_s = 0" \
    -o "${OUTDIR}/numu_s0.root" &

../../../bin/nuwro -i "parqel.txt" \
    -p "dyn_qel_cc = 0" \
    -p "dyn_qel_nc = 1" \
    -p "delta_s = -0.2" \
    -o "${OUTDIR}/numu_s-0.2.root" &

../../../bin/nuwro -i "parqel.txt" \
    -p "dyn_qel_cc = 0" \
    -p "dyn_qel_nc = 1" \
    -p "delta_s = 0.2" \
    -o "${OUTDIR}/numu_s+0.2.root"

../../../bin/reweight_to \
    "${OUTDIR}/numu_s0.root" \
    -o "${OUTDIR}/numu_s-0.2r.root" \
    -p delta_s -0.2

../../../bin/reweight_to \
    "${OUTDIR}/numu_s0.root" \
    -o "${OUTDIR}/numu_s+0.2r.root" \
    -p delta_s 0.2

../../../bin/myroot -b -q 'qel.c("'${OUTDIR}'/numu_s-0.2")' &
../../../bin/myroot -b -q 'qel.c("'${OUTDIR}'/numu_s+0.2")' &
../../../bin/myroot -b -q 'qel.c("'${OUTDIR}'/numu_s-0.2r")' &
../../../bin/myroot -b -q 'qel.c("'${OUTDIR}'/numu_s+0.2r")' &
