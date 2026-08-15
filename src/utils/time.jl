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

Convert the time delta `Δt` to a suitable `timestamp` of type `TT` relative to `reftime`. If `TT` is a `TimeType`
(i.e. `Date` or `DateTime`), then this method returns `reftime + Δt`.
"""
timestamp(::Type{TT}, reftime::TimeType, Δt::Period) where {TT <: TimeType} = convert(TT, reftime + Δt)
timestamp(::Type{TT}, reftime::TimeType, Δt::Number) where {TT <: TimeType} = convert(TT, reftime + convert_dt(Millisecond, Δt))
# If the requested type TT is a delta (number or period), convert the given Δt
timestamp(::Type{TT}, reftime::TimeType, Δt) where {TT <: Union{Period, Number}} = convert_dt(TT, Δt)
timestamp(::Type{TT}, reftime::TimeType, Δt::TT) where {TT <: Union{Period, Number}} = Δt # if Δt type already matches, return as is
timestamp(::Type{TT}, reftime::Number, Δt::Number) where {TT <: Union{Period, Number}} = convert_dt(TT, reftime + Δt)
# Fallback: convert the given timestamp to a delta
timestamp(::Type{TT}, reftime::TimeType, time::TimeType) where {TT <: Union{Period, Number}} = convert_dt(TT, time - reftime)
