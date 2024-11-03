# SPDX-License-Identifier: GPL-3.0-only
# Copyright © 2024 Roy Chih Chung Wang <roy.c.c.wang@proton.me>

function createorderingfromeqinds(eq_inds::Vector{Vector{Int}}, N::Int)::Tuple{Vector{Int},Int}

    j = 0 # degrees of freedom counter.
    out = zeros(Int, N)

    if isempty(eq_inds)
        return collect(1:N), N
    end

    # check if `eq_inds` only contains unique values from 1:N.
    if all( all( (eq_inds[k][l] in 1:N) for l in eachindex(eq_inds[k])) for k in eachindex(eq_inds) )
        for k in eachindex(eq_inds)
            j += 1

            out[eq_inds[k]] .= j
        end

        # fill the rest.
        for n = 1:N
            if out[n] == 0
                # out[n] is unassigned.
                j += 1
                out[n] = j
            end
        end

        return out, j
    end

    println("Invalid ME, using default.")
    return collect(1:N), N
end

function condensenuclei(
    x::Vector{T},
    ordering::Vector{Int},
    N::Int
    ) where T

    @assert length(x) == length(ordering)
    @assert norm(collect(1:N) - sort(unique(ordering))) < 1e-14

    y = zeros(T, N)
    for i in eachindex(x)

        k = ordering[i]
        y[k] += x[i] # add the coherence contributions among the equivalent resonances.
    end

    return y
end

function condensenuclei(
    x::Vector{T},
    ordering::Vector{Int}
    ) where T
    
    @assert length(x) == length(ordering)

    N = length(unique(ordering))

    return condensenuclei(x, ordering, N)
end

function reduceΔc(
    Δc_m::Vector{Vector{T}},
    ME_m::Vector{Vector{Int}},
    N_spins::Integer,
    cs::Vector{T},
    ) where T

    !isempty(ME_m) || error("ME_m cannot be empty.")

    # prepare.
    ordering, DOF = createorderingfromeqinds(ME_m, N_spins)

    cs_reduced = Vector{T}(undef, DOF)
    #nucs_reduced = Vector{Vector{Int}}(undef, DOF)
    for i in eachindex(ordering)
        k = ordering[i]
        cs_reduced[k] = cs[i]
        #nucs_reduced = nucs[i]
    end
    
    # condense.
    dc = Vector{Vector{T}}(undef, length(Δc_m))
    for l in eachindex(Δc_m) # over existing groups in Δc.
        dc[l] = condensenuclei(Δc_m[l], ordering, DOF)
    end

    return dc, cs_reduced
end


function getcoherencedropviolations(
    Δc::Vector{Vector{Vector{T}}};
    atol::T = convert(T, 1e-1),
    ) where T

    no_violations = true
    out = Vector{Tuple{Int,Int}}(undef, 0)
    for i in eachindex(Δc)

        for k in eachindex(Δc[i])
            if !isapprox(sum(Δc[i][k]), -one(T); atol = atol)
                push!(out, (i,k))
                no_violations = false
            end
        end
    end

    return out, no_violations
end

function getcoherencemagnitudeviolations(
    Δc::Vector{Vector{Vector{T}}};
    atol::T = convert(T, 1e-1),
    ) where T

    lb = -one(T) -atol
    ub = -lb

    no_violations = true

    out = Vector{Tuple{Int,Int}}(undef, 0)
    for i in eachindex(Δc)

        for k in eachindex(Δc[i])
            
            x = Δc[i][k]

            if !all(lb .< x .< ub)
                no_violations = false
                push!(out, (i,k))
            end
        end
    end

    return out, no_violations
end

function getNsys(A::SHType)::Int
    return length(A.Δc_bar)
end

# for a mixture
struct CoherenceDiagnostics

    inds_drop_Δc_m::Vector{Vector{Vector{Tuple{Int,Int}}}}
    inds_drop_Δc_bar::Vector{Vector{Vector{Tuple{Int,Int}}}}
    
    valid_drop_Δc_m::Bool
    valid_drop_Δc_bar::Bool

    inds_magnitude_Δc_m::Vector{Vector{Vector{Tuple{Int,Int}}}}
    inds_magnitude_Δc_bar::Vector{Vector{Vector{Tuple{Int,Int}}}}
    
    valid_magnitude_Δc_m::Bool
    valid_magnitude_Δc_bar::Bool
end

function extractdropstatus(A::CoherenceDiagnostics)
    return A.valid_drop_Δc_m, A.valid_drop_Δc_bar
end

function extractmagnitudestatus(A::CoherenceDiagnostics)
    return A.valid_magnitude_Δc_m, A.valid_magnitude_Δc_bar
end

function getΔcdiagnostics(
    As::Vector{SHType{T}};
    atol::T = convert(T, 1e-1),
    )::CoherenceDiagnostics where T

    N = length(As)
    
    inds_drop_Δc_m = Vector{Vector{Vector{Tuple{Int,Int}}}}(undef, N)
    inds_drop_Δc_bar = Vector{Vector{Vector{Tuple{Int,Int}}}}(undef, N)
    
    valid_drop_Δc_m = true
    valid_drop_Δc_bar = true

    inds_magnitude_Δc_m = Vector{Vector{Vector{Tuple{Int,Int}}}}(undef, N)
    inds_magnitude_Δc_bar = Vector{Vector{Vector{Tuple{Int,Int}}}}(undef, N)
    
    valid_magnitude_Δc_m = true
    valid_magnitude_Δc_bar = true

    for n in eachindex(As)
        N_sys = getNsys(As[n])

        inds_drop_Δc_m[n] = Vector{Vector{Tuple{Int,Int}}}(undef, N_sys)
        inds_drop_Δc_bar[n] = Vector{Vector{Tuple{Int,Int}}}(undef, N_sys)

        inds_magnitude_Δc_m[n] = Vector{Vector{Tuple{Int,Int}}}(undef, N_sys)
        inds_magnitude_Δc_bar[n] = Vector{Vector{Tuple{Int,Int}}}(undef, N_sys)

        for i in eachindex(As[n].Δc)

            inds_drop_Δc_m[n][i], status = getcoherencedropviolations(As[n].Δc)
            valid_drop_Δc_m = valid_drop_Δc_m && status

            inds_drop_Δc_bar[n][i], status = getcoherencedropviolations(As[n].Δc_bar)
            valid_drop_Δc_bar = valid_drop_Δc_bar && status

            inds_magnitude_Δc_m[n][i], status = getcoherencemagnitudeviolations(
                As[n].Δc;
                atol = atol,
            )
            valid_magnitude_Δc_m = valid_magnitude_Δc_m && status

            inds_magnitude_Δc_bar[n][i], status = getcoherencemagnitudeviolations(
                As[n].Δc_bar;
                atol = atol,
            )
            valid_magnitude_Δc_bar = valid_magnitude_Δc_bar && status
        end
    end

    return CoherenceDiagnostics(
        inds_drop_Δc_m,
        inds_drop_Δc_bar,
        valid_drop_Δc_m,
        valid_drop_Δc_bar,
        inds_magnitude_Δc_m,
        inds_magnitude_Δc_bar,
        valid_magnitude_Δc_m,
        valid_magnitude_Δc_bar,
    )
end


function checkcoherences(
    As::Vector{SHType{T}};
    atol::T = convert(T, 1e-1),
    ) where T
    
    coherence_diagnostics = getΔcdiagnostics(As; atol = atol)
    c_drop, bar_drop = extractdropstatus(coherence_diagnostics)
    c_mag, bar_mag = extractmagnitudestatus(coherence_diagnostics)
    
    #out::Bool = c_mag && bar_mag && c_drop && bar_drop
    out::Bool = c_drop && bar_drop # the sum(Δc) = -1 is more important than the magnitude < 1.
    
    return out, coherence_diagnostics
end

function checkcoherences(
    MSPs::Vector{MoleculeSpinSystem{T}};
    coherence_sum_zero_tol::T = convert(T, 1e-14),
    ) where T

    for n in eachindex(MSPs)
        for i in eachindex(MSPs[n].spin_systems)
            sp = MSPs[n].spin_systems[i]

            LHS = norm(sum.(sp.partial_quantum_numbers)-sp.quantum_numbers)
            if LHS > coherence_sum_zero_tol
                return false
            end
        end
    end

    return true
end