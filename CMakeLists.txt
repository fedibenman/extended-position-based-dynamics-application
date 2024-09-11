cmake_minimum_required(VERSION 3.5.0)
project(graphics VERSION 0.1.0 LANGUAGES C CXX)
find_package(OpenGL REQUIRED)

set(OpenGL_GL_PREFERENCE GLVND CACHE STRING "")
cmake_policy(SET CMP0077 NEW)

add_executable(graphics src/main.cpp) 

target_link_libraries(graphics glfw OpenGL::GL glm::glm)

find_package(glm CONFIG REQUIRED)

target_include_directories(graphics PRIVATE include)



set_property(TARGET graphics PROPERTY CXX_STANDARD 17)

target_compile_options(graphics PRIVATE
    $<$<CONFIG:Debug>:-g>
    $<$<CONFIG:RelWithDebInfo>:-g -fsanitize=address>
    $<$<CONFIG:Release>:-DNDEBUG>
)

include(CTest)
enable_testing()

