# Compile MDI
add_subdirectory(mdi)

# Add an __init__.py to the MDI directory, so that it can be used as a package
file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/mdi/__init__.py "")

# Copy the driver file into the compile directory
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/MDI_PyMeta.py ${CMAKE_CURRENT_BINARY_DIR}/MDI_PyMeta.py COPYONLY)
