set(TEST_INCLUDE_DIRECTORIES ${CMAKE_BINARY_DIR}/include ${CMAKE_SOURCE_DIR}/include ${CMAKE_CURRENT_SOURCE_DIR}/include)

include(CTest)
include(Catch)

file(GLOB_RECURSE test_sources RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "test_*.cpp")
foreach(testfile ${test_sources})
 get_filename_component(program_name ${testfile} NAME_WE)

 message(STATUS "Add test '${program_name}'")
 add_executable(${program_name} ${testfile})
 target_link_libraries(${program_name} PRIVATE ${DEPENDENCIES} ${TEST_DEPENDENCIES})
 target_include_directories(${program_name} PRIVATE ${TEST_INCLUDE_DIRECTORIES})
 catch_discover_tests(${program_name})
endforeach()
