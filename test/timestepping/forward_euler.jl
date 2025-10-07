using Terrarium

@testset "Forward Euler" begin
    Δt = 10.0
    euler = ForwardEuler(; Δt)
    @test !is_adaptive(euler)
    @test default_dt(euler) == Δt
    # Forward Euler is simple so we can just directly test it here
    progvar = 1.0
    tendency = 0.1
    @test Terrarium.step(euler, progvar, tendency, Δt) ≈ progvar + Δt*tendency
end