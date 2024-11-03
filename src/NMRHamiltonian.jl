# SPDX-License-Identifier: GPL-3.0-only
# Copyright © 2024 Roy Chih Chung Wang <roy.c.c.wang@proton.me>

module NMRHamiltonian

using LinearAlgebra
import Kronecker, Graphs
import JSON3

import SingleLinkagePartitions as SL

#using Serialization
using Statistics

# constant values.
function twopi(::Type{T}) where T <: AbstractFloat
    return T(2)*T(π)
end


include("./doc_strings/doc_types.jl")
include("./doc_strings/doc_configs.jl")


include("./types.jl")
include("./utils.jl")

# data conversion, disk input/output.
include("./IO/conversion_helpers.jl")
include("./IO/IO_mag_eq.jl")
include("./IO/IO_csJ.jl")
include("./IO/cs_parse.jl")
include("./IO/IO_phys.jl")


include("./IO/maximal_ME.jl") # generate Phys from an existing Phys, such that ME is maximized.

# spin Hamiltonian.
include("./SH/SH.jl")
include("./SH/Hamiltonian.jl")
include("./SH/operators.jl")
include("./SH/singlets.jl")


# overall simulator.
include("./engine/reduce_nuclei.jl")
include("./engine/partition.jl")
include("./engine/front_end.jl") # front end for scripts in the SH folder.

include("convinence.jl")

export 
assemble_physical_parameters,
PhysicalParamsType,
extract_cs,

get_partitions,
get_rcs,
get_group_cd,
get_component_cd,
num_groups,
get_Δc_matrix,
get_sorted_Δc_matrix,

simulate,
SHConfig,
SHType,

get_chem_shifts,
write_chem_shifts!,
extract_ME_nuclei,
get_num_groups

end
