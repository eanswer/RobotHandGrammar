set(SOURCE_DIR ${CMAKE_CURRENT_LIST_DIR})

file(GLOB_RECURSE SOURCE ${CMAKE_CURRENT_LIST_DIR}/*.h ${CMAKE_CURRENT_LIST_DIR}/*.hpp ${CMAKE_CURRENT_LIST_DIR}/*.cpp ${CMAKE_CURRENT_LIST_DIR}/*.natvis)

add_executable(design ${SOURCE})

target_link_libraries(design PRIVATE OpenGL::GL)
target_link_libraries(design PRIVATE glad)
target_link_libraries(design PRIVATE glfw)
target_link_libraries(design PRIVATE imgui)
target_link_libraries(design PRIVATE tetgen)
target_link_libraries(design PRIVATE pugixml)
if (UNIX)
target_link_libraries(design PRIVATE stdc++fs)
endif (UNIX)

if (MSVC)
target_link_directories(design PRIVATE ${EXTERNAL_LIBRARY})
endif (MSVC)

target_include_directories(design PRIVATE ${EXTERNAL_HEADER})
target_include_directories(design PRIVATE ${SOURCE_DIR})

if (WIN32)
	add_definitions(/bigobj)
endif()
