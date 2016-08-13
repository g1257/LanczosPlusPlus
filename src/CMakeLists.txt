cmake_minimum_required(VERSION 3.4)
project (LanczosPlusPlus)

cmake_policy(SET CMP0056 NEW)
cmake_policy(SET CMP0065 NEW)
cmake_policy(SET CMP0042 NEW)

message(${PsimagLite_SOURCE_DIR})

include_directories(${PsimagLite_SOURCE_DIR}
  ${PsimagLite_SOURCE_DIR}/../
  Engine
  )

add_executable (lanczos lanczos.cpp)
add_executable (thermal thermal.cpp)
add_executable (lorentzian lorentzian.cpp)

target_link_libraries (lanczos
  PsimagLite
  ${ADDITIONAL_LIBRARIES}
  )

target_link_libraries (thermal
  PsimagLite
  ${ADDITIONAL_LIBRARIES}
  )

target_link_libraries (lorentzian
  PsimagLite
  ${ADDITIONAL_LIBRARIES}
  )
