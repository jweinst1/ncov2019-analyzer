# examples cmake file

add_executable(DNATypes DNATypes.cpp)
add_test(DNATypeTests DNATypes test)

execute_process(COMMAND DNATypes show)