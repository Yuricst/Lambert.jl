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


# struct PlaceHoldCanonical <: CanonicalParameters
# 	lstar=1.0
# 	vstar=1.0
# 	tstar=1.0
# 	et0=0.0
# end