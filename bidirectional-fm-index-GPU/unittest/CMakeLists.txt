include_directories(${PROJECT_SOURCE_DIR}/src)

set(CMAKE_CXX_STANDARD 14) # google test requires c++14
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(unittest bitvectest.cpp bwtreprtest.cpp encodetextTest.cpp suffixArrayTest.cpp)
add_executable(consistencyUnitTest consistencyTest.cpp)

target_link_libraries(unittest AuxiliaryClasses)
target_link_libraries(consistencyUnitTest AuxiliaryClasses)

find_package(GTest)

if (NOT GTEST_FOUND)
    message(STATUS "GTest not found, fetching it using FetchContent.")

    
    include(FetchContent)
    FetchContent_Declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG v1.15.2
    )
    # For Windows: Prevent overriding the parent project's compiler/linker settings
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
    set(BUILD_GMOCK OFF CACHE BOOL "Disable Google Mock" FORCE)
    set(BUILD_GTEST ON CACHE BOOL "Enable Google Test" FORCE)

    FetchContent_MakeAvailable(googletest)


else()
    message(STATUS "GTest locally found.")
    # Print where the GTest libraries and headers are found
    message(STATUS "GTest Include Directory: ${GTEST_INCLUDE_DIR}")
    message(STATUS "GTest Library: ${GTEST_LIBRARY}")
    message(STATUS "GTest Main Library: ${GTEST_MAIN_LIBRARY}")

    include_directories(${GTEST_INCLUDE_DIRS})
    link_directories(${GTEST_LIBRARY_DIRS})
   
    target_link_libraries(unittest ${GTEST_LIBRARIES})
    target_link_libraries(consistencyUnitTest ${GTEST_LIBRARIES})

endif()

target_link_libraries(unittest gtest_main)
target_link_libraries(consistencyUnitTest gtest_main)

#optional linking for RLC
if(RUN_LENGTH_COMPRESSION)
    target_link_libraries(consistencyUnitTest  pthread)
endif()

if (ZLIB_FOUND)
        target_link_libraries(unittest ${ZLIB_LIBRARY})
        target_link_libraries(consistencyUnitTest ${ZLIB_LIBRARY})
endif (ZLIB_FOUND)

# Link SDSL if RUN_LENGTH_COMPRESSION is enabled and SDSL is found
if(RUN_LENGTH_COMPRESSION AND SDSL_FOUND)
    target_link_libraries(unittest ${SDSL_LIBRARY})

endif()

