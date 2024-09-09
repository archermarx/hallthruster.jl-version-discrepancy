using HallThruster
using CairoMakie

# Utility function for interactive work : just type "rl()" to reload sim
rl() = include("test.jl")

# Other utility functions
mean(x; start=1) = sum(x[i] for i in start:length(x)) / (length(x) - start + 1)

result = @timed HallThruster.run_simulation("variables.json") 
sim = result.value
time = result.time

start_iter = 500

discharge_current = HallThruster.discharge_current(sim)
discharge_current_mean = mean(discharge_current, start = start_iter)

ion_current = HallThruster.ion_current(sim)
ion_current_mean = mean(ion_current, start = start_iter)

thrust = HallThruster.thrust(sim)
thrust_mean = mean(thrust, start = start_iter) * 1000

using Printf

println("===========================")
println("          Summary")
println("===========================")
@printf("Runtime:           %.2f s\n", time)
@printf("Discharge current: %.2f A\n", discharge_current_mean)
@printf("Ion current:       %.2f A\n", ion_current_mean)
@printf("Thrust:            %.1f mN\n", thrust_mean)

