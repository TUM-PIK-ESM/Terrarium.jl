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
timedelta(::Type{TT}, reftime::TimeType, time::TimeType) where {TT} = convert_dt(TT, time - reftime)
# Fallback cases: the caller provided us with a Number or Period `time` which must already be a delta
# If the reference time is also a delta, we treat it as an offset; otherwise, ignore it
timedelta(::Type{TT}, Δt::Union{Number, Period}) where {TT} = timedelta(TT, zero(time), time)
timedelta(::Type{TT}, reftime::Union{Number, Period}, Δt::Union{Number, Period}) where {TT} = convert_dt(TT, Δt) - convert_dt(TT, reftime)
timedelta(::Type{TT}, reftime::TimeType, Δt::Union{Number, Period}) where {TT} = timestamp(TT, reftime, Δt)

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
timestamp(::Type{TT}, reftime::TimeType, time::TimeType) where {TT <: Union{Period, Number}} = timedelta(TT, reftime, time)
