#---------------------------------------------------
#QTINCDIR =$(wildcard /usr/include/qt4/ ${QTDIR}/include/)
#QTINCLUDEDIRS = -I. -I/usr/include/qt4/ -I/usr/include/qt4/QtCore -I/usr/include/qt4/QtGui -I${QTDIR}/include/QtCore -I${QTDIR}/include/QtGui
#QTINCLUDEDIRS = -I$(QTINCDIR)/ -I$(QTINCDIR)/QtCore -I$(QTINCDIR)/QtGui 
#QTLIBS =   -L/usr/lib -lQtGui -lQtCore -lpthread 

#DEBUG         = 1
#DEBUGON = -g
#CXXFLAGS      = `${ROOTSYS}/bin/root-config --cflags` -fPIC -O2 -I src 
CXXFLAGS      = `${ROOTSYS}/bin/root-config --cflags` -fPIC -O2 $(DEBUGON) -I src -Wl,--no-as-needed $(QTINCLUDEDIRS)
#LDFLAGS       = `${ROOTSYS}/bin/root-config --libs` -lPythia6 -lEG -lEGPythia6 -lCore  -lCint -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lGeom -lpthread -lm -ldl -rdynamic -lHist $(QTLIBS)
LDFLAGS       = `${ROOTSYS}/bin/root-config --libs` -lPythia6  -lEG -lEGPythia6 -lGeom -lMinuit $(QTLIBS)
LD	      = g++
CXX	      = g++
CC 	      = g++


TRGTS =         $(addprefix $(BIN)/,nuwro kaskada myroot glue event1.so nuwro2neut nuwro2nuance \
                dumpParams test_beam_rf test_makehist test_nucleus test_beam \
                fsi niwg ladek_topologies test mb_nce_run\
                )

DIS=    charge.o LeptonMass.o parameters.o grv94_bodek.o dis_cr_sec.o  dis_nc.o dis_cc_neutron.o delta.o dis2res.o \
		dis_cc_proton.o fragmentation.o fragmentation_nc.o fragmentation_cc.o singlepion.o \
	    disevent.o resevent2.o singlepionhadr.o alfa.o

SF_OBJS = $(patsubst %.cc,%.o,$(wildcard src/sf/*.cc))
GUI_OBJS = $(patsubst %.cc,%.o,$(wildcard src/gui/*.cc))
GUI_OBJS += $(patsubst src/gui/C%.cc,src/gui/moc_C%.o,$(wildcard src/gui/C*.cc))

DIS_OBJS= $(addprefix src/dis/,$(DIS))
BIN=bin
#BIN=.

EVENT_OBJS =  $(addprefix src/, event1.o event1dict.o pdg.o particle.o generatormt.o dirs.o)


all:            $(TRGTS)

#$(BIN)/boone1:         src/generatormt.o src/boone1.o
#		$(LINK.cc) $^ -o $@


$(BIN)/nuwro:   $(addprefix src/, event1.o event1dict.o generatormt.o particle.o pauli.o cohevent2.o cohdynamics2.o qelevent1.o mecdynamics.o mecevent.o\
        qel_sigma.o kinsolver.o kinematics.o pdg.o target_mixer.o nucleus.o  sfevent.o ff.o dirs.o rpa_lib.o nucleus_data.o isotopes.o elements.o \
        nuwro.o beam.o nd280stats.o beamHist.o coh.o fsi.o pitab.o scatter.o kaskada7.o Interaction.o main.o) \
        $(SF_OBJS) $(DIS_OBJS)
		$(LINK.cc) $^ -o $@

#src/gui/moc_%.cc:	src/gui/%.h 
#			moc $< -o $@

#QtNuwro:  $(addprefix src/, event1.o event1dict.o generatormt.o particle.o pauli.o cohevent2.o cohdynamics2.o qelevent1.o dirs.o \
#        qel_sigma.o kinsolver.o kinematics.o pdg.o scatter.o kaskada7.o  Interaction.o target_mixer.o nucleus.o  sfevent.o ff.o nuwro.o\
#        beam.o nd280stats.o beamHist.o coh.o fsi.o pitab.o) $(SF_OBJS) $(DIS_OBJS)  $(GUI_OBJS)\
#        gui/CMainWindow.o gui/moc_CMainWindow.o gui/QtNuwro.o gui/CGeneralTab.o gui/moc_CGeneralTab.o gui/CBeamTab.o gui/moc_CBeamTab.o \
#        gui/CTargetTab.o gui/moc_CTargetTab.o gui/CPhysicsTab.o gui/moc_CPhysicsTab.o gui/CElement.o gui/moc_CElement.o \
#        gui/CQel.o gui/moc_CQel.o gui/CRes.cc gui/moc_CRes.o gui/CCoh.o gui/moc_CCoh.o gui/CKaskada.o gui/moc_CKaskada.o 
#		$(LINK.cc) $^ -o $@


$(BIN)/kaskada:   $(addprefix src/, scatter.o generatormt.o particle.o event1.o event1dict.o kaskada7.o Interaction.o dirs.o\
				  pdg.o nucleus.o kaskada.o fsi.o pitab.o  nucleus_data.o isotopes.o elements.o)
		$(LINK.cc) $^ -o $@

$(BIN)/myroot:  $(EVENT_OBJS) src/myroot.o
		$(LINK.cc) $^ -o $@

$(BIN)/event1.so: $(EVENT_OBJS)
	        $(LD) -shared  $(LDFLAGS) $^ $(LIBS) -o $@

#$(BIN)/event1.a: $(EVENT_OBJS)
#	         ar r $@ $^  

$(BIN)/glue:    $(EVENT_OBJS) src/glue.o
		$(LINK.cc) $^ -o $@

$(BIN)/nuwro2neut:  $(EVENT_OBJS) src/nuwro2neut.o 
		$(LINK.cc) $^ -o $@

$(BIN)/nuwro2nuance: $(EVENT_OBJS) src/nuwro2nuance.o
		 $(LINK.cc) $^ -o $@

$(BIN)/fsi:   src/scatter.o src/generatormt.o src/particle.o src/event1.o src/event1dict.o src/kaskada7.o src/Interaction.o src/pdg.o src/dirs.o  src/nucleus.o  src/nucleus_data.o src/isotopes.o src/elements.o\
       src/fsi.o src/pitab.o src/calculations.o src/simulations.o src/vivisection.o src/plots.o  src/mplots.o  src/dirs.o src/fsi_main.o 
		$(LINK.cc) $^ -o $@

$(BIN)/mb_nce_run:   src/mb_nce.o src/mb_nce_run.o src/event1.o src/event1dict.o src/mb_nce_fit.o src/pdg.o src/scatter.o src/generatormt.o src/dirs.o src/particle.o
		$(LINK.cc) $^ -o $@

$(BIN)/niwg:   src/scatter.o src/generatormt.o src/particle.o src/event1.o src/event1dict.o src/kaskada7.o src/Interaction.o src/pdg.o src/dirs.o  src/nucleus.o  src/nucleus_data.o src/isotopes.o src/elements.o\
        src/fsi.o src/pitab.o src/calculations.o src/niwg_ccqe.o src/niwg_tech.o src/niwg_ccpi.o src/niwg.o  
		$(LINK.cc) $^ -o $@

$(BIN)/ladek_topologies: src/event1.o src/event1dict.o src/pdg.o src/particle.o  src/generatormt.o src/ladek_topologies.o src/dirs.o \
          src/fsi.o src/pitab.o
		$(LINK.cc) $^ -o $@

$(BIN)/test: src/event1.o src/event1dict.o src/pdg.o src/particle.o  src/generatormt.o src/test.o src/dirs.o
		$(LINK.cc) $^ -o $@

#$(BIN)/plots:           src/event1.o src/event1dict.o src/pdg.o src/particle.o src/generatormt.o src/dirs.o

$(BIN)/dumpParams:      src/dumpParams.o src/dirs.o
		$(LINK.cc) $^ -o $@

$(BIN)/test_nucleus:   src/generatormt.o src/nucleus.o src/test_nucleus.o src/pdg.o src/dirs.o  src/nucleus_data.o src/isotopes.o src/elements.o
		$(LINK.cc) $^ -o $@

$(BIN)/test_beam:	src/generatormt.o src/pdg.o src/test_beam.o 
		$(LINK.cc) $^ -o $@

$(BIN)/test_beam_rf:      src/test_beam_rf.o src/particle.o  src/generatormt.o 
		$(LINK.cc) $^ -o $@

$(BIN)/test_makehist:    src/test_makehist.o src/nd280stats.o
		$(LINK.cc) $^ -o $@

$(BIN)/test_balancer:       src/test_balancer.cc  src/generatormt.o
		$(LINK.cc) $^ -o $@

clean:;         @rm -f          *.o *.d src/event1dict.* core src/dis/*.o src/dis/*.d src/sf/*.o src/sf/*.d src/*.o src/*.d\
		src/gui/*.o src/gui/*.d src/gui/moc_*


distclean:;     @rm -f $(TRGTS) *.o *.d src/event1dict.* core src/dis/*.o src/dis/*.d src/sf/*.o src/sf/*.d src/*.o src/*.d\
		src/gui/*.o src/gui/*.d src/gui/moc_*  *.root *.root.txt


src/event1dict.h src/event1dict.cc:  src/params_all.h src/params.h src/event1.h src/event1LinkDef.h src/event1.o
		@echo "Generating dictionary ..."
		cd src;$(ROOTSYS)/bin/rootcint -f event1dict.cc -c event1.h event1LinkDef.h;cd ..
		

src/params_all.h:  src/params.xml src/params.h src/params.sed Makefile
		@echo "Building params_all.h"
		@echo "#define PARAMS_ALL()\\">src/params_all.h
#		@sed '/<!--.*-->/d;/<!--/,/-->/d;s/.*<param .*name="\([^"]*\)".*ctype="\([^"]*\)".*default="\([^"]*\)".*/PARAM(\2,\1,\3)\\/;tx;d;:x s/\(PARAM(\(vec\|line\|string\),[^,]*,\)\(.*\))\\/\1"\3")\\/'  src/params.xml >> src/params_all.h 
		@sed -f src/params.sed src/params.xml >> src/params_all.h 
		@echo "" >> src/params_all.h



%.d: %.cc
	@echo Making dependencies for $<
	@$(SHELL) -ec '$(CC) -MM -MT "$@ $<" $(CXXFLAGS) $< \
	| sed s/.cc\:/.o:/ > $@;\
	[ -s $@ ] || rm -f $@'


ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),distclean)
 -include $(patsubst %.cc,%.d,$(wildcard  src/dis/*.cc src/sf/*.cc src/*.cc src/gui/*.cc)) src/event1dict.d
endif
endif

# DO NOT DELETE
