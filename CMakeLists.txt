cmake_minimum_required (VERSION 2.6 FATAL_ERROR)

project (NUWRO)

enable_language(Fortran)

set (VERBOSE TRUE)

set (CMAKE_SKIP_BUILD_RPATH TRUE)

IF( NOT CMAKE_BUILD_TYPE )
   SET( CMAKE_BUILD_TYPE Release)
ENDIF()

execute_process(COMMAND git describe --tags 
          OUTPUT_VARIABLE NUWRO_VERSION 
          OUTPUT_STRIP_TRAILING_WHITESPACE)


set(CMAKE_INSTALL_PREFIX ${CMAKE_SOURCE_DIR})

if (NOT DEFINED ENV{ROOTSYS})
  message (FATAL_ERROR "$ROOTSYS is not defined, please set up root first.")
endif()

function(GenROOTDictionary)
endfunction()

execute_process (COMMAND root-config --cflags   OUTPUT_VARIABLE ROOT_CXX_FLAGS OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process (COMMAND root-config --libdir   OUTPUT_VARIABLE ROOT_LIBDIR    OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process (COMMAND root-config --version  OUTPUT_VARIABLE ROOT_VERSION   OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process (COMMAND root-config --features OUTPUT_VARIABLE ROOT_FEATURES  OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process (COMMAND root-config --libs     OUTPUT_VARIABLE ROOT_LIBS0     OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process (COMMAND root-config --ldflags  OUTPUT_VARIABLE ROOT_LD_FLAGS  OUTPUT_STRIP_TRAILING_WHITESPACE)


# remove libaries and directories  from ROOT_LIBS0 to obtain ROOT_LD_FLAGS
string(REGEX REPLACE " -(l|L)[^ ]* "  " "   CMAKE_LINK_FLAGS " ${ROOT_LIBS0} ") 

# obtain root lib directory
string(REGEX MATCH "-L([^ ]*)" DUMMY ${ROOT_LIBS0} )
set(ROOT_LIB_DIR ${CMAKE_MATCH_1}) 

# remove -l from ROOT_LIBS0 to obtain ROOT_LIBS
set(ROOT_LIBS0 " ${ROOT_LIBS0} -lPythia6  -lEG -lEGPythia6 -lGeom -lMinuit -lgfortran " )
string(REPLACE " -l" " "  ROOT_LIBS1 ${ROOT_LIBS0})
# remove non libraries from ROOT_LIBS
string(REGEX REPLACE " -[^ ]* "    " "  ROOT_LIBS1 " ${ROOT_LIBS1} ") 
# strip whitespace from ROOT_LIBS
string(STRIP "${ROOT_LIBS1}" ROOT_LIBS) 
# semicolon separated list of ROOT libraries
string(REPLACE " " ";"  ROOT_LIBS ${ROOT_LIBS})

message(STATUS  ROOT_LIBS ${ROOT_LIBS})

message ( STATUS "[ROOT]: root-config --version: "  ${ROOT_VERSION}  )
message ( STATUS "[ROOT]: root-config --cflags: "   ${ROOT_CXX_FLAGS})
message ( STATUS "[ROOT]: root-config --libdir: "   ${ROOT_LIBDIR} )
message ( STATUS "[ROOT]: root-config --features: " ${ROOT_FEATURES} )
message ( STATUS "[ROOT]: root-config --libs: "     ${ROOT_LIBS0} )
message ( STATUS "[ROOT]: root libraries: "         ${ROOT_LIBS1} )
message ( STATUS "[ROOT]: root lib dir: "           ${ROOT_LIB_DIR} )


############################### External Libs ##################################
set(EXTERNAL_LIBS ${ROOT_LIBS})
################################## COMPILER ####################################

message(STATUS "System:" ${CMAKE_SYSTEM_NAME})

set(CXX_WARNINGS "-Wall -Wno-unused-variable -Wno-sign-compare -Wno-unused-function    -Wno-unknown-warning-option -Wno-unused-but-set-variable -Wmissing-braces ")

if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
  set(CMAKE_CXX_FLAGS " ${CMAKE_CXX_FLAGS} ${ROOT_CXX_FLAGS} ${CXX_WARNINGS} ${NUWRO_CFLAGS} -std=c++0x -fPIC -O2 -Wl,--no-as-needed -DVERSION=\\\"${NUWRO_VERSION}\\\"  " )
endif()

if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  set(CMAKE_CXX_FLAGS " ${CMAKE_CXX_FLAGS} ${ROOT_CXX_FLAGS} ${CXX_WARNINGS} -Wmissing-braces ${NUWRO_CFLAGS}  -fPIC -O2  -DVERSION=\\\"${NUWRO_VERSION}\\\"  " )
endif()

message (STATUS "C++ Compiler      : " ${CXX_COMPILER_NAME})
message (STATUS "    Release flags : " ${CMAKE_CXX_FLAGS_RELEASE})
message (STATUS "    Debug flags   : " ${CMAKE_CXX_FLAGS_DEBUG})
message (STATUS "    Link Flags    : " ${CMAKE_LINK_FLAGS})


################################################################################
#                            Specify Target Subdirs
################################################################################

add_subdirectory(src)

################################################################################
#                            Generate setup.sh
################################################################################

file(WRITE  ${CMAKE_INSTALL_PREFIX}/setup.sh  "source $ENV{ROOTSYS}/bin/thisroot.sh\n")
file(APPEND ${CMAKE_INSTALL_PREFIX}/setup.sh "export NUWRO=${PROJECT_SOURCE_DIR}\n")
file(APPEND ${CMAKE_INSTALL_PREFIX}/setup.sh "export NUWRO_INC=${CMAKE_INSTALL_PREFIX}/include\n")
file(APPEND ${CMAKE_INSTALL_PREFIX}/setup.sh "export PATH=${CMAKE_INSTALL_PREFIX}/bin:\${PATH}\n")
if(DEFINED ENV{PYTHIA6})
  file(APPEND ${CMAKE_INSTALL_PREFIX}/setup.sh "export PYTHIA6=$ENV{PYTHIA6}\n")
  file(APPEND ${CMAKE_INSTALL_PREFIX}/setup.sh "export LD_LIBRARY_PATH=\${PYTHIA6}:\${LD_LIBRARY_PATH}\n")
endif()
file(APPEND ${CMAKE_INSTALL_PREFIX}/setup.sh "export LD_LIBRARY_PATH=${CMAKE_INSTALL_PREFIX}/lib:\${LD_LIBRARY_PATH}")

