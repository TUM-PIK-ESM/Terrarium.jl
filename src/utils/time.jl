"""
    $SIGNATURES

Convert `Δt`s to the given target `Period` or numeric type. Numeric types are always
assumed to be in seconds.
"""
convert_dt(Δt::NF) where {NF <: Number} = convert_dt(NF, Δt)
convert_dt(::Type{NF}, Δt::Number) where {NF <: Number} = convert(NF, Δt) # assume already in seconds
convert_dt(::Type{NF}, Δt::Period) where {NF <: Number} = Second(Δt).value
convert_dt(::Type{P}, Δt::NF) where {P <: Period, NF <: Number} = convert(P, Millisecond(round(Int, Δt * NF(1.0e3))))
convert_dt(::Type{P}, Δt::Period) where {P <: Period} = convert(P, Δt)

"""
    $SIGNATURES

Convert the given `time`stamp to a time delta relative to `reftime`.
"""
timedelta(time::Union{Number, Period}) = timedelta(zero(time), time)
timedelta(reftime, time) = timedelta(Second, reftime, time)
timedelta(::Type{TT}, time::Union{Number, Period}) where {TT} = timedelta(TT, zero(time), time)
timedelta(::Type{TT}, reftime::TimeType, time::TimeType) where {TT} = convert_dt(TT, time - reftime)
timedelta(::Type{TT}, reftime::Number, time::Number) where {TT} = convert_dt(TT, time - reftime)

"""
    $SIGNATURES

Convert the time delta `Δt` to a suitable `timestamp` of type `TT` relative to `reftime`. If `TT` is a `TimeType`
(i.e. `Date` or `DateTime`), then this method returns `reftime + Δt`. If `TT` is a `Period` or numeric type, `reftime`
is ignored and the original `Δt` converted to the units implied by `TT` is returned.
"""
timestamp(::Type{TT}, reftime::TimeType, Δt::Period) where {TT <: TimeType} = convert(TT, reftime + Δt)
timestamp(::Type{TT}, reftime::TimeType, Δt::Number) where {TT <: TimeType} = convert(TT, reftime + convert_dt(Millisecond, Δt))
timestamp(::Type{TT}, reftime::TimeType, Δt) where {TT <: Union{Period, Number}} = convert_dt(TT, Δt)
timestamp(::Type{TT}, reftime::TimeType, Δt::TT) where {TT <: Union{Period, Number}} = Δt # if Δt type already matches, return as is
timestamp(::Type{TT}, reftime::TT, Δt::TT) where {TT <: Union{Period, Number}} = reftime + Δt
