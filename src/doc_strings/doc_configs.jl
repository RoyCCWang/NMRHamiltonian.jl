# SPDX-License-Identifier: GPL-3.0-only
# Copyright © 2024 Roy Chih Chung Wang <roy.c.c.wang@proton.me>




DOCSTRING_cs_sys_mixture(T, show_type = true) = """
$(show_type ? "`::Vector{Vector{Vector{$(T)}}}`" : "")


The set of chemical shift values for non-singlet spin systems for multiple molecule entries. See `PhysicalParamsType`.
"""

DOCSTRING_cs_singlets_mixture(T, show_type = true) = """
$(show_type ? "`::Vector{Vector{$(T)}}`" : "")


Singlet chemical shift values for multiple molecule entries. See `PhysicalParamsType`.
"""

DOCSTRING_Phys(T, show_type = true) = """
$(show_type ? "`::Vector{PhysicalParamsType{$(T)}}`" : "")


Physical parameters for multiple molecule entries.
"""



############## functions.

