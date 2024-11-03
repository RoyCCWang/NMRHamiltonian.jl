# SPDX-License-Identifier: GPL-3.0-only
# Copyright © 2024 Roy Chih Chung Wang <roy.c.c.wang@proton.me>


function combinevectors(x::Vector{Vector{T}}) where T

    if isempty(x)
        return Vector{T}(undef, 0)
    end

    N = sum(length(x[i]) for i in eachindex(x))

    y = Vector{T}(undef,N)

    st_ind = 0
    fin_ind = 0
    for i in eachindex(x)
        st_ind = fin_ind + 1
        fin_ind = st_ind + length(x[i]) - 1

        y[st_ind:fin_ind] = x[i]
    end

    return y
end

function convertcompactdomain(x::T, a::T, b::T, c::T, d::T)::T where T <: AbstractFloat

    return (x-a)*(d-c)/(b-a)+c
end

function convertcompactdomain(x::Vector{T}, a::T, b::T, c::T, d::T) where T <: AbstractFloat

    return collect( convertcompactdomain(x[i], a, b, c, d) for i in eachindex(x) )
end

