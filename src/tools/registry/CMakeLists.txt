if (DEFINED ENV{iAMAS_TOOL_DIR})
message(STATUS "*** Using iAMAS tools/parse from $ENV{iAMAS_TOOL_DIR}/registry ***")
foreach(_core IN LISTS iAMAS_CORES)
  set(PARSE_PATH $ENV{iAMAS_TOOL_DIR}/registry/parse)
  file(COPY ${PARSE_PATH} DESTINATION ${CMAKE_BINARY_DIR}/bin FILE_PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ)
  file(RENAME ${CMAKE_BINARY_DIR}/bin/parse ${CMAKE_BINARY_DIR}/bin/parse_${_core})
  add_custom_target(parse_${_core} ALL DEPENDS ${CMAKE_BINARY_DIR}/bin/parse_${_core})
endforeach()

else()
message(STATUS "*** Building iAMAS tools/parse from source ***")
#Parsing library core-independent code
add_library(parselib dictionary.c fortprintf.c utility.c)
target_link_libraries(parselib PUBLIC ${PROJECT_NAME}::external::ezxml)
target_link_libraries(parselib PUBLIC ${PROJECT_NAME}::external::esmf)

# Generate parser for each core
#
# Note: One parser is required per-core because the gen_inc.c depends on
# a pre-processor define iAMAS_NAMELIST_SUFFIX which is core specific
foreach(_core IN LISTS iAMAS_CORES)
    add_executable(parse_${_core} parse.c gen_inc.c)
    target_link_libraries(parse_${_core} PUBLIC parselib)
    target_compile_definitions(parse_${_core} PRIVATE iAMAS_NAMELIST_SUFFIX=${_core}
                                                           iAMAS_GIT_VERSION=${iAMAS_GIT_VERSION}
                                                           iAMAS_EXE_NAME=${_core}_model)
endforeach()
endif()
