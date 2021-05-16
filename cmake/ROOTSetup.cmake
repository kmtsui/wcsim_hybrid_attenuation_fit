if ( NOT DEFINED ENV{ROOTSYS} )
  cmessage (FATAL_ERROR "$ROOTSYS is not defined, please set up ROOT first.")
else()
  cmessage(STATUS "Using ROOT installed at $ENV{ROOTSYS}")
  set(CMAKE_ROOTSYS $ENV{ROOTSYS})
endif()

# Get cflags from ROOT
execute_process (COMMAND root-config
  --cflags OUTPUT_VARIABLE ROOT_CXX_FLAGS_RAW OUTPUT_STRIP_TRAILING_WHITESPACE)
string(REPLACE " " ";" ROOT_CXX_FLAGS "${ROOT_CXX_FLAGS_RAW}")
# Get libdir from ROOT
execute_process (COMMAND root-config
  --libdir OUTPUT_VARIABLE ROOT_LIBDIR OUTPUT_STRIP_TRAILING_WHITESPACE)
# Get version from ROOT
execute_process (COMMAND root-config
  --version OUTPUT_VARIABLE ROOT_VERSION OUTPUT_STRIP_TRAILING_WHITESPACE)
# Get features from ROOT
execute_process (COMMAND root-config
  --features OUTPUT_VARIABLE ROOT_FEATURES OUTPUT_STRIP_TRAILING_WHITESPACE)

string(REGEX MATCH "^6.*" ROOT_SIX ${ROOT_VERSION})

LIST(APPEND ROOT_LIBS
  Core
  RIO
  Hist
  -Wl,--whole-archive HistPainter -Wl,--no-whole-archive
  Graf
  Graf3d
  Gpad
  Tree
  Matrix
  Physics
  MathCore)

#Check what ROOT thinks the standard is, set that project-wide
# and then remove it from ROOT_CXX_FLAGS
list (FIND ROOT_CXX_FLAGS "-std=c++11" CPP11_INDEX)
list (FIND ROOT_CXX_FLAGS "-std=c++14" CPP14_INDEX)
list (FIND ROOT_CXX_FLAGS "-std=c++17" CPP17_INDEX)
if (${CPP11_INDEX} GREATER -1)
  SET(CMAKE_CXX_STANDARD 11)
elseif (${CPP14_INDEX} GREATER -1)
  SET(CMAKE_CXX_STANDARD 14)
elseif (${CPP17_INDEX} GREATER -1)
  SET(CMAKE_CXX_STANDARD 17)
else()
  SET(CMAKE_CXX_STANDARD 11)
endif()
list(REMOVE_ITEM ROOT_CXX_FLAGS "-std=c++11")
list(REMOVE_ITEM ROOT_CXX_FLAGS "-std=c++14")
list(REMOVE_ITEM ROOT_CXX_FLAGS "-std=c++17")

cmessage ( STATUS "[ROOT]: root-config --version: ${ROOT_VERSION} ")
cmessage ( STATUS "[ROOT]: root-config --cflags : ${ROOT_CXX_FLAGS} ")
cmessage ( STATUS "[ROOT]: root-config --libdir : ${ROOT_LIBDIR} ")

link_directories(${ROOT_LIBDIR})

function(LinkToROOT Target)
  get_target_property(ALREADY_LINKED ${Target} LINKED_TO_ROOT)

  if(NOT "${ALREADY_LINKED} " STREQUAL "1 " )
    target_compile_options(${Target} PUBLIC ${ROOT_CXX_FLAGS})    

    target_link_libraries(${Target} ${ROOT_LIBS})
    set_target_properties(${Target} PROPERTIES LINKED_TO_ROOT 1)
  endif()
endfunction()

#Helper functions for building dictionaries
function(GenROOTDictionary OutputDictName Header LinkDef)

  get_directory_property(incdirs INCLUDE_DIRECTORIES)
  string(REPLACE ";" ";-I" LISTDIRINCLUDES "-I${incdirs};${CMAKE_CURRENT_SOURCE_DIR}")
  string(REPLACE " " ";" LISTCPPFLAGS "${ROOT_CXX_FLAGS}")

  message(STATUS "LISTCPPFLAGS: ${LISTCPPFLAGS}")
  message(STATUS "LISTINCLUDES: ${LISTDIRINCLUDES}")
  #Learn how to generate the Dict.cc and Dict.hxx
  cmessage(STATUS "Outputname: ${CMAKE_CURRENT_BINARY_DIR}/${OutputDictName}.cc")
  add_custom_command(
    OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/${OutputDictName}.cc" "${CMAKE_CURRENT_BINARY_DIR}/${OutputDictName}.h"
    COMMAND rootcint
    ARGS -f ${OutputDictName}.cc -c
    -p ${LISTDIRINCLUDES} ${LISTCPPFLAGS} ${Header} ${LinkDef}
    DEPENDS ${Header};${LinkDef})
endfunction()
