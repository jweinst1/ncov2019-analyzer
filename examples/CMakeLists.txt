# examples cmake file

add_executable(DNATypes DNATypes.cpp)
add_test(DNATypeTests DNATypes test)

# Way to run the example 
add_custom_command(TARGET DNATypes POST_BUILD COMMAND DNATypes show )

add_executable(DNATrie DNATrie.cpp)