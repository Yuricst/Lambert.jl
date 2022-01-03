"""
Canonical parameters class
"""

abstract type CanonicalParameters end

struct SunEarthCanonical <: CanonicalParameters
	lstar::Float64
	vstar::Float64
	tstar::Float64
	et0::Float64
end