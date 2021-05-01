include_directories(src)
include_directories(include)
file(GLOB HEADERS include/*.h)
file(GLOB SOURCES src/*.cpp)

add_executable(${PROJECTNAME} ${SOURCES} main.cpp)

target_link_libraries(${PROJECTNAME}
                      ${PCL_LIBRARIES}
                      ${CGAL_LIBRARIES}
                      ${OpenCV_LIBS}
                      ${Boost_LIBRARIES}
                      ${SuiteSparse_CHOLMOD_LIBRARY_RELEASE})
