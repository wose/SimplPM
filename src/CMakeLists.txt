set(simplpm_SRC
  basics.cpp
  common.cpp
  deep.cpp
  get_el.cpp
  sdp4.cpp
  sdp8.cpp
  sgp.cpp
  sgp4.cpp
  sgp8.cpp
  )

set(simplpm_HDR
  norad.h
  norad_in.h
  )

add_library(simplpm ${simplpm_SRC})
target_link_libraries(simplpm m)
