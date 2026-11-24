if (DEFINED ENV{iAMAS_TOOL_DIR})
message(STATUS "*** Using iAMAS namelist&streams_gen from $ENV{iAMAS_TOOL_DIR}/input_gen ***")

#add_executable(iamas_namelist_gen IMPORTED)
#set_target_properties(iamas_namelist_gen PROPERTIES IMPORTED_LOCATION "$ENV{iAMAS_TOOL_DIR}/input_gen/namelist_gen")
#add_executable(${PROJECT_NAME}::iamas_namelist_gen ALIAS iamas_namelist_gen)
#
#add_executable(iamas_streams_gen IMPORTED)
#set_target_properties(iamas_streams_gen PROPERTIES IMPORTED_LOCATION "$ENV{iAMAS_TOOL_DIR}/input_gen/streams_gen")
#add_executable(${PROJECT_NAME}::iamas_streams_gen ALIAS iamas_streams_gen)

set(NAMELIST_GEN_PATH $ENV{iAMAS_TOOL_DIR}/input_gen/namelist_gen)
file(COPY ${NAMELIST_GEN_PATH} DESTINATION ${CMAKE_BINARY_DIR}/bin FILE_PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ)
file(RENAME ${CMAKE_BINARY_DIR}/bin/namelist_gen ${CMAKE_BINARY_DIR}/bin/namelist_gen)
add_custom_target(namelist_gen ALL DEPENDS ${CMAKE_BINARY_DIR}/bin/namelist_gen)

set(STREAMS_GEN_PATH $ENV{iAMAS_TOOL_DIR}/input_gen/streams_gen)
file(COPY ${STREAMS_GEN_PATH} DESTINATION ${CMAKE_BINARY_DIR}/bin FILE_PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ)
file(RENAME ${CMAKE_BINARY_DIR}/bin/streams_gen ${CMAKE_BINARY_DIR}/bin/streams_gen)
add_custom_target(streams_gen ALL DEPENDS ${CMAKE_BINARY_DIR}/bin/streams_gen)
else()
message(STATUS "*** Building iAMAS namelist&streams_gen from source ***")
add_executable(namelist_gen namelist_gen.c test_functions.c)
target_link_libraries(namelist_gen PUBLIC ${PROJECT_NAME}::external::ezxml)

add_executable(streams_gen streams_gen.c test_functions.c)
target_link_libraries(streams_gen PUBLIC ${PROJECT_NAME}::external::ezxml)
endif()
