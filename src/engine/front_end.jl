# SPDX-License-Identifier: GPL-3.0-only
# Copyright © 2024 Roy Chih Chung Wang <roy.c.c.wang@proton.me>

"""
```
simulate(
    Phys::Vector{PhysicalParamsType{T}},
    molecule_entries::Vector{String},
    fs::T,
    SW::T,
    ν_0ppm::T,
    config::SHConfig;
) where T <: AbstractFloat
```

Simulate the resonance components' intensities and frequencies, compute order of coherences for each component, and group the components into resonance groups.

### Inputs

- `Phys` -- A `Vector` of physical chemistry parameters, e.g. the output of assemble_physical_parameters().
- `molecule_entries` -- A `Vector` of compound entries.
- `fs` -- the sampling frequency in Hz for use in the simulation.
- `SW` -- the spectral window in ppm for use in the simulation.
- `ν_0ppm` -- the 0 ppm peak frequency in Hz in the spectrum, for use in the simulation.
- `config` -- a configuration file of type `SHConfig`.

### Outputs

- `As::Vector{SHType{Float64}}` -- the simulated resonance intensities and frequencies, sub-system order of coherences, and resonance groups.

- `MSPs` -- the spin Hamiltonian matrices and coherence-related quantities. For diagnostic purposes.

"""
function simulate(
    Phys::Vector{PhysicalParamsType{T}},
    molecule_entries::Vector{String},
    fs::T,
    SW::T,
    ν_0ppm::T,
    configs::Vector{SHConfig{T}},
    ) where T <: AbstractFloat

    # set up.
    N_molecules = length(molecule_entries)
    @assert length(Phys) == N_molecules

    ppm2hzfunc = pp->(ν_0ppm + pp*fs/SW)

    # assemble.
    As = Vector{SHType{T}}(undef, N_molecules)
    #Rs = Vector{Vector{PartitionSearchRecord{T}}}(undef, N_molecules)
    MSPs = Vector{MoleculeSpinSystem{T}}(undef, N_molecules)

    # check
    if length(unique(molecule_entries)) != length(molecule_entries)
        println("Error: Please supply a list of unique compound molecule_entries. Exit without simulation.")
        return As, MSPs
    end
    
    # simulate.
    for n in eachindex(Phys)

        # create quantities based on non-singlets.
        αs, Ωs, parts, Δc, Δc_bar, N_spins_sys, cs_Δc, ME_nucs,
        αs_singlets, Ωs_singlets, N_spins_singlet, 
        MSPs[n] = setupmoleculeSH(Phys[n], ppm2hzfunc, configs[n])

        # incorporate info from singlets.
        setupsingletspinsystem!(
            αs,
            Ωs,
            Δc,
            parts,
            Δc_bar,
            cs_Δc,
            ME_nucs, # here and above mutates if there are singlets.

            N_spins_sys, # this and the above mutates.
            Phys[n].cs_singlets,
            Phys[n].H_inds_singlets,
            Phys[n].H_IDs,
            
            # here, N_spins_singlet gets appended to N_spins_sys, as with αs_singlets to αs, etc.
            αs_singlets,
            Ωs_singlets,
            N_spins_singlet,
        )

        # package up.
        As[n] = SHType(
            αs,
            Ωs,
            Δc,
            parts,
            Δc_bar,
            cs_Δc,
            ME_nucs,
            N_spins_sys,
            #αs_singlets, Ωs_singlets,
            fs, SW, ν_0ppm
        )
        
        #Rs[n] = search_results
    end

    return As, MSPs
end

# use the same config for all compounds.
function simulate(
    Phys::Vector{PhysicalParamsType{T}},
    molecule_entries::Vector{String},
    fs::T,
    SW::T,
    ν_0ppm::T,
    config::SHConfig;
    ) where T <: AbstractFloat

    return simulate(
        Phys, molecule_entries, fs, SW, ν_0ppm,
        collect( config for _ = 1:length(molecule_entries))
    )
end

function setupmoleculeSH(
    phy::PhysicalParamsType{T},
    ppm2hzfunc,
    config::SHConfig,
    ) where T <: AbstractFloat

    cs_sys = phy.cs_sys

    ### spin Hamiltonian simulation.
    N_spins_singlet = length.(phy.H_inds_singlets)

    N_spins_sys = collect( length(cs_sys[m]) for m in eachindex(cs_sys) )
    intermediates_sys = prepcouplingalgorithm(T, N_spins_sys)
    
    # αs_inp, Ωs_inp, coherence_mat_sys, eig_vals_sys, Q_sys,
    # coherence_state_pairs_sys, H_sys, states_sys, ms_sys, M_sys = computeSH(
    MSP = computeSH( # MSP stands for molecule spin system.
        cs_sys,
        phy.J_vals_sys,
        phy.J_inds_sys_local,
        intermediates_sys,
        ppm2hzfunc,
        phy.cs_singlets,
        N_spins_singlet,
        config,
    )
    
    αs, Ωs, parts, Δc, Δc_bar, cs_Δc, ME_nucs, _ = partitionresonances(
        MSP.spin_systems,
        N_spins_sys,
        config,
        phy.ME,
        cs_sys,
        phy.H_IDs,
        phy.H_inds_sys,
    )

    return αs, Ωs, parts, Δc, Δc_bar, N_spins_sys, cs_Δc, ME_nucs,
    MSP.singlet_intensities, MSP.singlet_frequencies, N_spins_singlet, MSP
end


####################### for convinence. Based on presets from real data.

function getpresetspectrometer(::Type{T}, tag) where T <: AbstractFloat

    fs = T(14005.602240896402)
    SW = T(20.0041938620844)
    ν_0ppm = T(10656.011933076665)

    if tag == "700"
        # machine values taken from the BMRB 700 MHz 20 mM glucose experiment.
        fs = T(14005.602240896402)
        SW = T(20.0041938620844)
        ν_0ppm = T(10656.011933076665)

    elseif tag == "600"
        ## machine values from a 600 MHz experiment: bmse000915, methionine.
        fs = T(9615.38461538462)
        SW = T(16.022093454391)
        ν_0ppm = T(6685.791496181313)

    elseif tag == "900"
        ## machine values from a 900 MHz experiment.
        fs = T(14423.0769230769)
        SW = T(16.0300195009073)
        ν_0ppm = T(10160.027322585376)

    elseif tag == "500"
        fs = T(6493.50649350649)
        SW = T(12.9911090156122)
        ν_0ppm = T(4035.6644246816795)

    elseif tag == "400"

        ### 400 MHz, bmse000297, ethanol.
        fs = T(4807.69230769231)
        SW = T(12.0152693165838)
        ν_0ppm = T(2884.905244600881)
    end

    return fs, SW, ν_0ppm
end

function loadandsimulate(
    fs::T,
    SW::T,
    ν_0ppm::T,
    molecule_entries::Vector{String},
    H_params_path::String,
    molecule_mapping_file_path::String,
    config::SHConfig{T};
    unique_cs_digits::Int = 6,
    )::Tuple{
        Vector{PhysicalParamsType{T}},
        Vector{SHType{T}},
        Vector{MoleculeSpinSystem{T}},
    } where T <: AbstractFloat

    Phys = assemble_physical_parameters(
        T,
        molecule_entries,
        H_params_path,
        molecule_mapping_file_path;
        unique_cs_digits = unique_cs_digits,
    )

    As, MSPs = simulate(
        Phys,
        molecule_entries,
        fs,
        SW,
        ν_0ppm,
        config,
    )

    return Phys, As, MSPs
end

# Utilities
function get_num_groups(As::Vector{NMRHamiltonian.SHType{T}}) where T <: AbstractFloat
    
    N_groups = 0
    for A in As
        for Δc_bar in A.Δc_bar
            N_groups += length(Δc_bar)
        end
    end
    return N_groups
end