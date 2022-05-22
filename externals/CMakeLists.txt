set(EXTERNAL_HEADER
  "${CMAKE_CURRENT_LIST_DIR}"
  "${CMAKE_CURRENT_LIST_DIR}/eigen"
  "${CMAKE_CURRENT_LIST_DIR}/glad/include"
  "${CMAKE_CURRENT_LIST_DIR}/glfw/include"
  "${CMAKE_CURRENT_LIST_DIR}/tetgen"
  "${CMAKE_CURRENT_LIST_DIR}/pugixml/src"
)

if (MSVC)
list(APPEND EXTERNAL_HEADER "${CMAKE_CURRENT_LIST_DIR}/freeglut/include")
set(EXTERNAL_LIBRARY
  "${CMAKE_CURRENT_LIST_DIR}/freeglut/lib/x64"
)
endif (MSVC)

# Adds library target imgui
add_subdirectory(libigl-imgui)
# Adds library target glfw
add_subdirectory(glfw)
# Adds library target glad
add_subdirectory(glad)
# Adds Library target tetgen
add_subdirectory(tetgen)
# pugixml
add_subdirectory(pugixml)

