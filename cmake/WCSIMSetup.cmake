if ( NOT DEFINED ENV{WCSIMDIR} )
  cmessage (FATAL_ERROR "$WCSIMDIR is not defined, please set up WCSIM first.")
else()
  cmessage(STATUS "Using WCSIM installed at $ENV{WCSIMDIR}")
  set(CMAKE_WCSIMDIR $ENV{WCSIMDIR})
endif()

SET(WCSIM_CXX_FLAGS "-I${CMAKE_WCSIMDIR}/src;-I${CMAKE_WCSIMDIR}/include")
SET(WCSIM_LIBDIR "${CMAKE_WCSIMDIR}")

LIST(APPEND WCSIM_LIBS
  WCSimRoot
)

cmessage ( STATUS "[WCSIM]: WCSIM_CXX_FLAGS : ${WCSIM_CXX_FLAGS} ")
cmessage ( STATUS "[WCSIM]: WCSIM_LIBDIR : ${WCSIM_LIBDIR} ")
cmessage ( STATUS "[WCSIM]: WCSIM_LIBS : ${WCSIM_LIBS} ")

link_directories(${WCSIM_LIBDIR})

function(LinkToWCSIM Target)
  get_target_property(ALREADY_LINKED ${Target} LINKED_TO_WCSIM)

  if(NOT "${ALREADY_LINKED} " STREQUAL "1 " )
    target_compile_options(${Target} PUBLIC ${WCSIM_CXX_FLAGS})    

    target_link_libraries(${Target} ${WCSIM_LIBS})
    set_target_properties(${Target} PROPERTIES LINKED_TO_WCSIM 1)
  endif()
endfunction()