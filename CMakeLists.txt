SET (src_dir ${CMAKE_CURRENT_SOURCE_DIR})

SET (coalescer_sources
  ${src_dir}/CoalescerMN.cc
  ${src_dir}/CoalescerMN.h
  PARENT_SCOPE
)

SET (coalescer_data
  ${src_dir}/CoalescerMN.xsd
  PARENT_SCOPE
)
