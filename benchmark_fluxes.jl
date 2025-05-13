
using HallThruster: HallThruster as het

include("profile.jl")
include("calc_fluxes.jl")

#sim = run(50, 3, 1e-5)

function allocate_arrays(sim)
    scheme = sim.config.scheme
    F = zeros(size(sim.params.cache.F))
    UL = copy(F)
    UR = copy(F)
    U = copy(sim.params.cache.U)
    _, ncells = size(U)

    (;index, fluids, grid, ncharge, is_velocity_index) = sim.params
    (; λ_global, dt_u) = sim.params.cache
    cache = (; λ_global, dt_u)
    params = (;index, fluids, grid, cache, ncharge, is_velocity_index)

    continuity_fluids = [
        ContinuityFluid2(fluids[1].species, fluids[1].u, fluids[1].T, ncells-2)
    ]

    isothermal_fluids = [
        IsothermalFluid3(f.species, f.T, ncells-2) for f in fluids[2:end]
    ]

    continuity_fluids[1].density .= U[1, :]
    for (i, f) in enumerate(isothermal_fluids)
        f.density .= U[2*i, :]
        f.momentum .= U[2*i+1, :]
    end

    fluid_containers = (;
        continuity = continuity_fluids,
        isothermal = isothermal_fluids,
    )

    return (; F, UL, UR, U, params, scheme), fluid_containers
end

fluid_containers, F, F_new = let
    sim = run(100, 3, 1e-4)
    (; F, UL, UR, U, params, scheme), fluid_containers = allocate_arrays(sim)

    reconstruct = true
    scheme = het.HyperbolicScheme(
        scheme.flux_function, scheme.limiter, reconstruct
    )
    grid = sim.params.grid
    #@btime het.compute_fluxes!($F, $UL, $UR, $U, $params, $scheme)
    @time het.compute_fluxes!(F, UL, UR, U, params, scheme)

    (;isothermal, continuity) = fluid_containers
    F_new = zeros(size(F))
    #@btime compute_fluxes!($F_new, $continuity, $isothermal, $grid, $scheme.limiter, $reconstruct)
    @time compute_fluxes!(F_new, continuity, isothermal, grid, scheme.limiter, reconstruct)
    
    #@profview for i in 1:100000
    #    het.compute_fluxes!(F, UL, UR, U, params, scheme)
    #end

    fluid_containers, F, F_new
end;



