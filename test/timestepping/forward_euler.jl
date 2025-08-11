using Terra

@testset "Forward Euler" begin
    dt = 10.0
    euler = ForwardEuler(; dt)
    @test !is_adaptive(euler)
    @test get_dt(euler) == dt
    # Forward Euler is simple so we can just directly test it here
    progvar = 1.0
    tendency = 0.1
    @test Terra.step(euler, progvar, tendency, dt) ≈ progvar + dt*tendency
end