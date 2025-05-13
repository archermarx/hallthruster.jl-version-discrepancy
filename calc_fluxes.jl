using HallThruster: HallThruster as het


abstract type Fluid end
abstract type AbstractIsothermalFluid <: Fluid end
abstract type AbstractContinuityFluid <: Fluid end

struct IsothermalFluid3 <: AbstractIsothermalFluid
    # Conservative variables
    density::Vector{Float64}
    momentum::Vector{Float64}

    # Edge states
    dens_L::Vector{Float64}
    dens_R::Vector{Float64}
    mom_L::Vector{Float64}
    mom_R::Vector{Float64}    

    # Fluxes
    flux_dens::Vector{Float64}
    flux_mom::Vector{Float64}

    # Wave speed
    wave_speed::Array{Float64, 0}

    # Data
    species::het.Species
    sound_speed::Float64

    function IsothermalFluid3(species, temp, num_cells)
        R = species.element.R
        γ = species.element.γ
        a = sqrt(γ * R * temp)

        return new(
            # Conservative variables
            zeros(num_cells+2),
            zeros(num_cells+2),

            # Edge states
            zeros(num_cells+1),
            zeros(num_cells+1),
            zeros(num_cells+1),
            zeros(num_cells+1),

            # Fluxes
            zeros(num_cells+1),
            zeros(num_cells+1),

            # Wave speed
            fill(a),

            # Data
            species,
            a
        )
    end
end

struct ContinuityFluid2 <: AbstractContinuityFluid
    # Conservative varibales
    density::Vector{Float64}

    # Edge states
    dens_L::Vector{Float64}
    dens_R::Vector{Float64}

    # Fluxes
    flux_dens::Vector{Float64}

    # Wave speed
    wave_speed::Array{Float64, 0}

    # Data
    species::het.Species
    sound_speed::Float64
    velocity::Float64

    function ContinuityFluid2(species, vel, temp, num_cells)
        R = species.element.R
        γ = species.element.γ
        a = sqrt(γ * R * temp)

        return new(
            # Conservative variables
            zeros(num_cells+2),

            # Edge states
            zeros(num_cells+1),
            zeros(num_cells+1),

            # Fluxes
            zeros(num_cells+1),

            # Wave speed
            fill(max(abs(vel+a), abs(vel-a))),

            # Data
            species,
            a,
            vel,
        )
    end
end

@inline function reconstruct(uⱼ₋₁, uⱼ, uⱼ₊₁, limiter)
    r = (uⱼ₊₁ - uⱼ) / (uⱼ - uⱼ₋₁)
    Δu = 0.25 * limiter(r) * (uⱼ₊₁ - uⱼ₋₁)
    return uⱼ - Δu, uⱼ + Δu
end

function compute_edge_states!(fluid::AbstractContinuityFluid, limiter, do_reconstruct)
    (;density, dens_L, dens_R) = fluid
    N = length(fluid.density)

    if do_reconstruct
        @inbounds for i in 2:N-1
            iL, iR = het.left_edge(i), het.right_edge(i)

            # Reconstruct density
            u₋ = density[i - 1]
            uᵢ = density[i]
            u₊ = density[i + 1]
            dens_R[iL], dens_L[iR] = reconstruct(u₋, uᵢ, u₊, limiter,)
        end
    else
        @inbounds for i in 2:N-1
            iL, iR = het.left_edge(i), het.right_edge(i)
            dens_L[iR] = density[i]
            dens_R[iL] = density[i]
        end
    end

    fluid.dens_L[1] = fluid.density[1]
    fluid.dens_R[end] = fluid.density[end]
end

function compute_edge_states!(fluid::AbstractIsothermalFluid, limiter, do_reconstruct)
    (;density, momentum, dens_L, dens_R, mom_L, mom_R) = fluid
    N = length(fluid.density)

    if do_reconstruct
        @inbounds for i in 2:N-1
            iL, iR = het.left_edge(i), het.right_edge(i)

            # Reconstruct density
            u₋ = density[i - 1]
            uᵢ = density[i]
            u₊ = density[i + 1]
            dens_R[iL], dens_L[iR] = reconstruct(u₋, uᵢ, u₊, limiter,)

            # Reconstruct velocity as primitive variable instead of momentum density
            u₋ = momentum[i-1] / u₋
            uᵢ = momentum[i] / uᵢ
            u₊ = momentum[i+1] / u₊
            uR, uL = reconstruct(u₋, uᵢ, u₊, limiter)
            mom_L[iR] = uL * dens_L[iR]
            mom_R[iL] = uR * dens_R[iL]
        end
    else
        @inbounds for i in 2:N-1
            iL, iR = het.left_edge(i), het.right_edge(i)
            dens_L[iR] = density[i]
            dens_R[iL] = density[i]
            mom_L[iR] = momentum[i]
            mom_R[iL] = momentum[i]
        end
    end

    fluid.dens_L[1] = fluid.density[1]
    fluid.dens_R[end] = fluid.density[end]
    fluid.mom_L[1] = fluid.momentum[1]
    fluid.mom_R[end] = fluid.momentum[end]
end

function compute_fluxes!(
    F,
    continuity::Vector{T},
    isothermal::Vector{U},
    grid,
    limiter,
    do_reconstruct;
    apply_boundary_conditions = false
    ) where {T<:AbstractContinuityFluid, U <: AbstractIsothermalFluid}

    # Compute edge states (not including boundary conditions)
    for (i, fluid) in enumerate(continuity)
        compute_edge_states!(fluid, limiter, do_reconstruct)
    end

    for (i, fluid) in enumerate(isothermal)
        compute_edge_states!(fluid, limiter, do_reconstruct)
    end

    # Compute fluxes
    for (i, fluid) in enumerate(continuity)
        compute_fluxes!(fluid, grid)
        F[i, :] .= fluid.flux_dens
    end

    start = length(continuity)
    for fluid in isothermal
        compute_fluxes!(fluid, grid)
        F[start+1, :] .= fluid.flux_dens
        F[start+2, :] .= fluid.flux_mom
        start += 2
    end
end

function compute_fluxes!(fluid::AbstractContinuityFluid, grid)
    (;flux_dens, dens_L, dens_R, wave_speed, velocity) = fluid
    smax = wave_speed[]
    @inbounds for i in eachindex(grid.edges)
        ρ_L, ρ_R = dens_L[i], dens_R[i]
        flux_dens[i] = 0.5 * (velocity * (ρ_L + ρ_R) - smax * (ρ_R - ρ_L))
    end
end

function compute_fluxes!(fluid::AbstractIsothermalFluid, grid)
    (;flux_dens, flux_mom, dens_L, dens_R, mom_L, mom_R, wave_speed) = fluid
    a = fluid.sound_speed
    RT = a^2 / fluid.species.element.γ

    max_wave_speed = 0.0

    @inbounds for i in eachindex(grid.edges)
        ρ_L, ρ_R = dens_L[i], dens_R[i]
        ρu_L, ρu_R = mom_L[i], mom_R[i]

        u_L, u_R = ρu_L / ρ_L, ρu_R / ρ_R

        smax = max(abs(u_L - a), abs(u_L + a), abs(u_R - a), abs(u_R + a))
        max_wave_speed = max(smax, max_wave_speed)

        flux_mom_L = ρ_L * (u_L^2 + RT) 
        flux_mom_R = ρ_R * (u_R^2 + RT)

        flux_dens[i] = 0.5 * ((ρu_L + ρu_R) - smax * (ρ_R - ρ_L))
        flux_mom[i] = 0.5 * ((flux_mom_L + flux_mom_R) - smax * (ρu_R - ρu_L))
    end

    wave_speed[] = max_wave_speed
end



