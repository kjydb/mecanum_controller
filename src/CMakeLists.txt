project(mecanum_controller LANGUAGES CXX)

set(MECANUM_SOURCE_FILES
    mecanum_controller.cpp
    main.cpp)

add_executable(mecanum ${MECANUM_SOURCE_FILES})
target_include_directories(mecanum PUBLIC ${CMAKE_SOURCE_DIR}/include)
