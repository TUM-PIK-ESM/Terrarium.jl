"""
    $TYPEDSIGNATURES

Return a `KernelFunctionOperation` that computes the skin temperature residual over all grid cells.
"""
diagnose_skin_temperature_residual(state, model::LandModel) = diagnose_skin_temperature_residual(state, model, model.surface_energy_balance, model.snow)
function diagnose_skin_temperature_residual(
        state,
        model::AbstractModel,
        seb::SurfaceEnergyBalance,
        snow::Optional{AbstractSnow} = nothing
    )
    return kernel_operation2D(state, model.grid, seb.skin_temperature, model.constants, snow) do i, j, grid, fields, skinT::ImplicitSkinTemperature, args...
        Ts_implicit = compute_skin_temperature(i, j, grid, fields, skinT, args...)
        Ts_prev = state.skin_temperature[i, j, 1]
        return Ts_prev - Ts_implicit
    end
end
