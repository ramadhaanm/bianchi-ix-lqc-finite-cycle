"""Memoize immutable Vern9 coefficients for Double64 without changing their bits.

OrdinaryDiffEqVerner 2.4.2's generic constructors rebuild BigFloat-derived constants
on every dense interpolation. Event guards call interpolation very frequently.
This optional execution optimization preserves the original numerical algorithm.
"""
module VernerDoubleCache
using OrdinaryDiffEqVerner, DoubleFloats
const V = OrdinaryDiffEqVerner
const interpolation = V.Vern9InterpolationCoefficients(Double64)
const tableau = V.Vern9Tableau(Double64, Double64)
const extra = V.Vern9ExtraStages(Double64, Double64)
all(isbitstype(typeof(x)) for x in (interpolation,tableau,extra)) || error("Coefficient tables must be immutable bits types")

@eval V Vern9InterpolationCoefficients(::Type{$Double64}) = $interpolation
@eval V Vern9Tableau(::Type{$Double64}, ::Type{$Double64}) = $tableau
@eval V Vern9ExtraStages(::Type{$Double64}, ::Type{$Double64}) = $extra

"""Compare cached values against the original generic constructor, byte for byte."""
function verify()
    pairs = (
        (interpolation, invoke(V.Vern9InterpolationCoefficients,Tuple{Type},Double64)),
        (tableau, invoke(V.Vern9Tableau,Tuple{Type,Type},Double64,Double64)),
        (extra, invoke(V.Vern9ExtraStages,Tuple{Type,Type},Double64,Double64)))
    all(reinterpret(UInt8,[a]) == reinterpret(UInt8,[b]) for (a,b) in pairs)
end
verify() || error("Cached coefficients differ from the original constructors")
end
