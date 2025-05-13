using HallThruster: HallThruster as het
using BenchmarkTools

function run(num_cells::Int = 100, max_charge::Int = 1, duration=1e-3)
    config = het.Config(
        ncharge = max_charge,
        thruster = het.SPT_100,
        domain = (0.0, 0.08),
        discharge_voltage = 300.,
        anode_mass_flow_rate = 5e-6,
        anom_model = het.ScaledGaussianBohm(
            anom_scale = 0.0625,
            barrier_scale = 0.98,
            width = 0.25,
            center = 1.0,
        )
    )

    sim_params = het.SimParams(
        num_save = round(Int, 1000 * duration / 1e-3),
        adaptive = true,
        grid = het.EvenGrid(num_cells),
        dt = 1e-9,
        duration = duration,
        verbose=false
    )

    sol = het.run_simulation(config, sim_params)

    return sol
end

#run(200)

# for ncells in [50, 100, 150, 200]
#     for ncharge in 1:3
#         stats = @timed run(ncells, ncharge)
#         println("$(ncells) cells, $(ncharge) charges: $(stats.time) seconds")
#     end
# end

# @profview run(100, 3, 1e-3)

#@benchmark run(200)