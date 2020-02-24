
INCLUDE_DIRECTORIES( ${PROJECT_SOURCE_DIR}/include )
FILE(GLOB SC_FILES "./PROST/*.c" "./PROST/*.h" "*.H" "*.C")

# generate library
add_library (eosH2ONaCl ${SC_FILES})