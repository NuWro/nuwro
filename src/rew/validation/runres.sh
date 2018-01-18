#!/bin/bash

##### SETTINGS #####

ENERGY_MIN=1000
ENERGY_MAX=5000
ENERGY_STEP=1000

MA_MIN=900
MA_MAX=1100
MA_STEP=100

C5A_MIN=1.1
C5A_MAX=1.3
C5A_STEP=0.1

OUTDIR="res/"
mkdir -p $OUTDIR

##### SIMULATIONS ######

for energy in $(seq $ENERGY_MIN $ENERGY_STEP $ENERGY_MAX)
do
    for ma in $(seq $MA_MIN $MA_STEP $MA_MAX)
        for c5a in $(seq $C5A_MIN $C5A_STEP $C5A_MAX)
        do
            ../../../bin/nuwro -i "parres.txt" \
                -p "beam_particle = 14" \
                -p "beam_energy = $energy" \
                -p "pion_axial_mass = $ma" \
                -p "pion_C5A = $c5a" \
                -p "dyn_qel_cc = 1" \
                -p "dyn_qel_nc = 0" \
                -o "${OUTDIR}/numu${energy}_cc_ma${ma}_c5a${c5a}.root" &

            ../../../bin/nuwro -i "parres.txt" \
                -p "beam_particle = 14" \
                -p "beam_energy = $energy" \
                -p "pion_axial_mass = $ma" \
                -p "pion_C5A = $c5a" \
                -p "dyn_qel_cc = 0" \
                -p "dyn_qel_nc = 1" \
                -o "${OUTDIR}/numu${energy}_nc_ma${ma}_c5a${c5a}.root" &

            ../../../bin/nuwro -i "parres.txt" \
                -p "beam_particle = -14" \
                -p "beam_energy = $energy" \
                -p "pion_axial_mass = $ma" \
                -p "pion_C5A = $c5a" \
                -p "dyn_qel_cc = 1" \
                -p "dyn_qel_nc = 0" \
                -o "${OUTDIR}/numubar${energy}_cc_ma${ma}_c5a${c5a}.root" &

            ../../../bin/nuwro -i "parres.txt" \
                -p "beam_particle = 14" \
                -p "beam_energy = $energy" \
                -p "pion_axial_mass = $ma" \
                -p "pion_C5A = $c5a" \
                -p "dyn_qel_cc = 0" \
                -p "dyn_qel_nc = 1" \
                -o "${OUTDIR}/numubar${energy}_nc_ma${ma}_c5a${c5a}.root" &
        done
    done
done

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