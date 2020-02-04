
aux_source_directory(. DIR_LIB_SRCS)
include_directories(./ ../PROST ../freesteam)

# generate library
add_library (eosH2ONaCl ${DIR_LIB_SRCS})
set_property(TARGET eosH2ONaCl PROPERTY IMPORTED_LOCATION /Users/zguo/MyData/Research/3-CodeProject/PROST-4.1/libsteam4.a)
