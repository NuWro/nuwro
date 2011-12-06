#!/bin/bash
for targ in targFe_Ben targAr_Ben; do # targC_Ben targO_Ben targO_GSF targCa_GSF targAr_GSF targFe_Ben; do
#  for momDistrib in md_C12_Ben, md_O16_CdA, md_O16_Ben, md_O16_GCo, md_Ca40_CdA, md_Ca40_GCo, md_Ca48_GCo; do
    for neutrino_flavour in flav_e flav_mu flav_ae flav_amu; do # flav_mless,  flav_e,  flav_mu,  flav_tau, flav_amless, flav_ae, flav_amu, flav_atau; do
      for form_factors in dipoleFF BBBA05FF BBA03FF JLabFF; do
        for switchtoEM in 0; do
          for switchPauliBlocking in 1; do
            for neutrino_energy in 800; do
              for number_of_repetitions in 7200000; do
  ./Main -tg $targ -nf $neutrino_flavour -ff $form_factors -em $switchtoEM -pb $switchPauliBlocking -ne $neutrino_energy -nr $number_of_repetitions # -md $momDistrib 
              done
            done
          done
        done
      done
    done
#  done
done