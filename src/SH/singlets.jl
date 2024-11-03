# SPDX-License-Identifier: GPL-3.0-only
# Copyright © 2024 Roy Chih Chung Wang <roy.c.c.wang@proton.me>

# ### not pursuing the unification of singlets and spin systems because solvent might need a separate model.
# a separate singlet model might be easier to extend for the solvent in the future.
function setupsingletspinsystem!(
    αs::Vector{Vector{T}},
    Ωs::Vector{Vector{T}},
    Δc::Vector{Vector{Vector{T}}},
    parts::Vector{Vector{Vector{Int}}},
    Δc_bar::Vector{Vector{Vector{T}}},
    cs_Δc::Vector{Vector{T}},
    ME_nucs::Vector{Vector{Vector{Int}}}, # here and above all mutates.
    
    N_spins_sys::Vector{Int}, # this and the above mutates.
    cs_singlets::Vector{T},
    H_inds_singlets::Vector{Vector{Int}},
    H_IDs::Vector{Int},

    #η::Vector{Vector{Vector{T}}},
    αs_singlets::Vector{T},
    Ωs_singlets::Vector{T},
    N_spins_singlet::Vector{Int},
    ) where T <: AbstractFloat

    @assert length(αs_singlets) == length(Ωs_singlets) == length(N_spins_singlet)

    # :αs_singlets, :Ωs_singlets, :fs, :SW, :ν_0ppm

    if isempty(αs_singlets)
        # no singlets in current molecule.
        return nothing
    end

    push!(N_spins_sys, N_spins_singlet...)

    for i in eachindex(αs_singlets, cs_singlets, H_inds_singlets)
        
        ### constant params.
        
        singlet_α = Vector{T}(undef, 1) # one resonance component in a singlet.
        singlet_α[begin] = αs_singlets[i]
        push!(αs, singlet_α)

        singlet_Ω = Vector{T}(undef, 1) # one resonance component in a singlet.
        singlet_Ω[begin] = Ωs_singlets[i]
        push!(Ωs, singlet_Ω)

        # one spin system (first index), one resonance group (second index)
        singlet_Δc = Vector{Vector{T}}(undef, 1)
        singlet_Δc[begin] = -ones(T, 1)

        singlet_part_inds = Vector{Vector{Int}}(undef, 1)
        singlet_part_inds[begin] = ones(Int, 1)

        push!(Δc, singlet_Δc)
        push!(Δc_bar, singlet_Δc) # since only one resonance component, set Δc_bar to Δc.
        push!(parts, singlet_part_inds)

        push!(cs_Δc, [cs_singlets[i];])

        # tmp_inds = Vector{Vector{Int}}(undef, 1)
        # tmp_inds[1] = copy(H_inds_singlets[i])
        tmp = Vector{Vector{Int}}(undef, 1)
        #tmp[1] = #collect( [H_IDs[h_id];] for h_id in H_inds_singlets[i])
        tmp[1] = H_IDs[H_inds_singlets[i]]
        push!(ME_nucs, tmp)
    end

    return nothing
end