###
#
# @copyright (c) 2009-2014 The University of Tennessee and The University 
#                          of Tennessee Research Foundation. 
#                          All rights reserved.
# @copyright (c) 2012-2014 Inria. All rights reserved.
# @copyright (c) 2012-2014 IPB. All rights reserved. 
#
###
#
#  @file CMakeLists.txt
#
#  @project MORSE
#  MORSE is a software package provided by:
#     Inria Bordeaux - Sud-Ouest,
#     Univ. of Tennessee,
#     King Abdullah Univesity of Science and Technology
#     Univ. of California Berkeley,
#     Univ. of Colorado Denver. 
# 
#  @version 1.1.0
#  @author Cedric Castagnede
#  @author Emmanuel Agullo
#  @author Mathieu Faverge
#  @date 13-07-2012
#   
###

cmake_minimum_required(VERSION 2.8)

# Looking for doxygen
# -------------------
find_package(Doxygen)

# Generate documentation
# ----------------------
if(DOXYGEN_EXECUTABLE)
    # Create file magma_morse.dox
    # ---------------------------
    if(MORSE_SCHED_STARPU)
        if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/../../../scratch)
            configure_file("conf-dev-starpu.dox.in"
                           "magma_morse.dox"
                           @ONLY)
        else()
            configure_file("conf-usr-starpu.dox.in"
                           "magma_morse.dox"
                           @ONLY)
        endif()

    elseif(MORSE_SCHED_QUARK)
        if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/../../../scratch)
            configure_file("conf-dev-quark.dox.in"
                           "magma_morse.dox"
                            @ONLY)
        else()
            configure_file("conf-usr-quark.dox.in"
                           "magma_morse.dox"
                           @ONLY)
        endif()

    else()
        if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/../../../scratch)
            configure_file("conf-dev-starpu.dox.in"
                           "magma_morse.dox"
                           @ONLY)
        else()
            configure_file("conf-usr-starpu.dox.in"
                           "magma_morse.dox"
                           @ONLY)
        endif()

    endif()

    # Message during configure
    #-------------------------
    message(STATUS "Looking for doxygen - found")

    # User's documentation
    # --------------------
    add_custom_command(OUTPUT  ${CMAKE_CURRENT_BINARY_DIR}/out
                       COMMAND ${DOXYGEN_EXECUTABLE}
                       ARGS    ${CMAKE_CURRENT_BINARY_DIR}/magma_morse.dox
                       DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/magma_morse.dox
                      )
    add_custom_target(doxygen-out ALL
                      DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/out)

    # Compilation and installation for latex version
    # ----------------------------------------------
    find_package(LATEX)
    if(PDFLATEX_COMPILER AND MAKEINDEX_COMPILER)
        # Message
        # -------
        message(STATUS "Looking for pdflatex - found")
        message(STATUS "Looking for makeindex - found")

        # Compile pdf
        # -----------
        add_custom_command(OUTPUT            refman.pdf
                           COMMAND           ${CMAKE_MAKE_PROGRAM}
                           ARGS              refman.pdf
                           WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/out/latex
                           )
        add_custom_target(doxygen-user_guide-magma_morse ALL
                          DEPENDS refman.pdf ${CMAKE_CURRENT_BINARY_DIR}/out) 

        # Installation of produced pdf
        # ----------------------------
        if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/../../../scratch)
            install(FILES       ${CMAKE_CURRENT_BINARY_DIR}/out/latex/refman.pdf
                    RENAME      magma-morse-devguide-${MAGMA_VERSION_NUMBER}.pdf
                    DESTINATION share/pdf)
        else()
            install(FILES       ${CMAKE_CURRENT_BINARY_DIR}/out/latex/refman.pdf
                    RENAME      magma-morse-usersguide-${MAGMA_VERSION_NUMBER}.pdf
                    DESTINATION share/pdf)
        endif()

    endif()

    # Installation for html version
    # -----------------------------
    install(DIRECTORY   ${CMAKE_CURRENT_BINARY_DIR}/out/html
            DESTINATION share/doxygen)

    # Installation for man version
    # ----------------------------
    install(DIRECTORY   ${CMAKE_CURRENT_BINARY_DIR}/out/man
            DESTINATION share/)

else()
    message(STATUS "Looking for doxygen - not found")

endif()

###
### END CMakeLists.txt
###
