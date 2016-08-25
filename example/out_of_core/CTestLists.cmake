#
# Check Example out_of_core
#

set(TESTLIST
    out_of_core
    )

foreach(test ${TESTLIST})
    add_test(example_ooc_${test} ./${test})
endforeach()
