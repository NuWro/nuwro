
add_custom_command(OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/params_all.h
  COMMAND echo '\#define PARAMS_ALL\(\)\\' > ${CMAKE_CURRENT_SOURCE_DIR}/params_all.h
  COMMAND sed -f ${CMAKE_CURRENT_SOURCE_DIR}/params.sed ${CMAKE_CURRENT_SOURCE_DIR}/params.xml >> ${CMAKE_CURRENT_SOURCE_DIR}/params_all.h
  COMMAND echo '' >> ${CMAKE_CURRENT_SOURCE_DIR}/params_all.h
  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/params.xml ${CMAKE_CURRENT_SOURCE_DIR}/params.sed
)

add_custom_target(params_header DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/params_all.h)

set(DIS charge.cc LeptonMass.cc parameters.cc grv94_bodek.cc dis_cr_sec.cc  dis_nc.cc dis_cc_neutron.cc delta.cc dis2res.cc 
	 dis_cc_proton.cc fragmentation.cc fragmentation_nc.cc fragmentation_cc.cc singlepion.cc 
   disevent.cc resevent2.cc singlepionhadr.cc alfa.cc res_kinematics.cc res_xsec.cc
)

set(COMMON_SOURCES  
  cohevent2.cc cohdynamics2.cc  qelevent1.cc mecdynamics.cc mecevent.cc mecevent_tem.cc mecevent_Nieves.cc mecevent_SuSA.cc particle.cc
  mecevent_common.cc mecdynamics2.cc mecevent2.cc qel_sigma.cc kinsolver.cc input_data.cc data_container.cc kinematics.cc pdg.cc target_mixer.cc
  nucleus.cc sfevent.cc ff.cc rpa_2013.cc nucleus_data.cc isotopes.cc elements.cc beam.cc nd280stats.cc beamHist.cc coh.cc fsi.cc pitab.cc scatter.cc
  kaskada7.cc hipevent.cc Interaction.cc e_el_event.cc e_el_sigma.cc e_spp_event.cc  espp/crossection.cc espp/DM.cc rew/PythiaQuiet.f
  pauli.cc generatormt.cc rew/rewparams.cc
)

set(EVENT_SOURCES event1.cc generatormt.cc particle.cc pdg.cc dirs.cc  event1Dict.cc )

set(FSI_SOURCES scatter.cc kaskada7.cc Interaction.cc input_data.cc data_container.cc nucleus.cc  nucleus_data.cc isotopes.cc elements.cc fsi.cc pitab.cc )  

set(REW_SOURCES  rew/Reweighters.cc rew/rewQEL.cc rew/rewRES.cc rew/rewNorm.cc rew/PythiaQuiet.f)

foreach(FILE ${DIS})
  set(DIS_SOURCES ${DIS_SOURCES} dis/${FILE})
endforeach()

set(ESPP_SURCES e_spp_event.cc)

AUX_SOURCE_DIRECTORY(espp ESPP_SOURCES)

AUX_SOURCE_DIRECTORY(sf SF_SOURCES)

include_directories($ENV{ROOTSYS}/include .)
link_directories(${ROOT_LIB_DIR};$ENV{PYTHIA6})
link_libraries(${EXTERNAL_LIBS})

add_custom_command(  
 OUTPUT ${CMAKE_SOURCE_DIR}/src/event1Dict.cc 
 WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/src/
 COMMAND rootcint
 ARGS   -f event1Dict.cc -c event1.h event1LinkDef.h 
 DEPENDS ${CMAKE_SOURCE_DIR}/src/event1.h ${CMAKE_SOURCE_DIR}/src/event1LinkDef.h ${CMAKE_SOURCE_DIR}/src/vec.h
)

add_library(event SHARED ${EVENT_SOURCES})
add_dependencies(event params_header)

add_library(event1 STATIC ${EVENT_SOURCES})
add_dependencies(event1 params_header)

add_library(common1 STATIC ${COMMON_SOURCES})
add_dependencies(common1 params_header)

add_library(fsi1 STATIC ${FSI_SOURCES} pdg.cc)
add_dependencies(fsi1 params_header)

add_library(sf1 STATIC ${SF_SOURCES} pdg.cc)
add_dependencies(sf1 params_header)

add_library(dis1 STATIC ${DIS_SOURCES})
add_dependencies(dis1 params_header)

add_library(espp1 STATIC ${ESPP_SOURCES})
add_dependencies(espp1 params_header)

add_library(rew1 STATIC ${REW_SOURCES})
add_dependencies(rew1 params_header)
target_link_libraries(rew1 ${EXTERNAL_LIBS} common1 dis1)



add_executable(nuwro main.cc nuwro.cc)
target_link_libraries(nuwro ${EXTERNAL_LIBS} event1 common1 dis1 espp1 sf1)

add_executable(kaskada scatter.cc kaskada7.cc Interaction.cc input_data.cc data_container.cc 
                nucleus.cc kaskada.cc fsi.cc pitab.cc nucleus_data.cc isotopes.cc elements.cc) 
target_link_libraries(kaskada ${EXTERNAL_LIBS} event1)

add_executable(myroot myroot.cc event1Dict.cc)
add_dependencies(myroot params_header)
target_link_libraries(myroot ${EXTERNAL_LIBS} event1 )

add_executable(nuwro2neut nuwro2neut.cc)
target_link_libraries(nuwro2neut ${EXTERNAL_LIBS} event1)

add_executable(nuwro2nuance nuwro2nuance.cc )
target_link_libraries(nuwro2nuance ${EXTERNAL_LIBS} event1)

add_executable(nuwro2rootracker nuwro2rootracker.cc )
target_link_libraries(nuwro2rootracker ${EXTERNAL_LIBS} event1)

add_executable(glue glue.cc  )
target_link_libraries(glue ${EXTERNAL_LIBS} event1)

add_executable(dumpParams dumpParams.cc dirs.cc)

add_executable(reweight_to rew/reweight_to.cc)
target_link_libraries(reweight_to ${EXTERNAL_LIBS} event1 common1 fsi1 dis1 sf1 rew1 )

add_executable(reweight_along rew/reweight_along.cc)
target_link_libraries(reweight_along ${EXTERNAL_LIBS} event1 common1 fsi1 rew1 sf1 dis1)

add_executable(whist whist.cc)
target_link_libraries(whist ${EXTERNAL_LIBS} event1)

add_executable(fsi  fsi_main.cc calculations.cc simulations.cc plots.cc mplots.cc vivisection.cc)
target_link_libraries(fsi ${EXTERNAL_LIBS} event1 fsi1 common1)

add_executable(ladek_topologies ladek_topologies.cc fsi.cc pitab.cc)
target_link_libraries(ladek_topologies ${EXTERNAL_LIBS} event1)

add_executable(ganalysis  ganalysis.cc nuwro.cc)
target_link_libraries(ganalysis ${EXTERNAL_LIBS} event1 fsi1 common1 dis1 rew1 sf1 espp1)

add_executable(mb_nce_run   mb_nce_run.cc mb_nce.cc mb_nce_fit.cc scatter.cc)
target_link_libraries(mb_nce_run ${EXTERNAL_LIBS} event1)

add_executable(niwg niwg.cc niwg_ccqe.cc niwg_tech.cc niwg_ccpi.cc niwg.cc calculations.cc )
target_link_libraries(niwg ${EXTERNAL_LIBS} event1 common1 )

add_executable(test_nucleus  generatormt.cc nucleus.cc test_nucleus.cc  nucleus_data.cc isotopes.cc elements.cc)
target_link_libraries(test_nucleus ${EXTERNAL_LIBS} event1)

add_executable(test_beam	test_beam.cc)
target_link_libraries(test_beam ${EXTERNAL_LIBS} event1)

add_executable(test_beam_rf test_beam_rf.cc )
target_link_libraries(test_beam_rf ${EXTERNAL_LIBS} event1 )

add_executable(test_makehist test_makehist.cc nd280stats.cc)
target_link_libraries(test_makehist ${EXTERNAL_LIBS} event1 )



install(TARGETS nuwro kaskada myroot nuwro2neut nuwro2nuance nuwro2rootracker  
        dumpParams glue reweight_to reweight_along 
        whist niwg fsi ganalysis mb_nce_run ladek_topologies 
        test_nucleus test_beam test_beam_rf test_makehist
        DESTINATION bin)
install(TARGETS event DESTINATION lib)
install(FILES event1Dict_rdict.pcm DESTINATION bin OPTIONAL) # OPTIONAL necessary for ROOT5
install(FILES event1Dict_rdict.pcm DESTINATION lib OPTIONAL) # OPTIONAL necessary for ROOT5
install(DIRECTORY ../data/  DESTINATION data) # data files are needed for nuwro run 
install(FILES params.h vec.h dirs.h params_all.h particle.h generatormt.h 
              jednostki.h  pdg.h vect.h vecrand.h event1.h
        DESTINATION include) # may be useful for ROOT6
