# examples cmake file

add_executable(DNATypes DNATypes.cpp)
add_test(DNATypeTests DNATypes test)

# Way to run the example 
add_custom_command(TARGET DNATypes POST_BUILD COMMAND DNATypes show )

add_definitions(-DTRIE_TEST_DEBUG)
add_executable(DNATrie DNATrie.cpp)
add_test(DNATrieTests DNATrie test)
# visual the test.
add_custom_command(TARGET DNATrie POST_BUILD COMMAND DNATrie test )