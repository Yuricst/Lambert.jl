"""
Spice helpers
"""


"""
	get_spice_locate_function(
		body_id::Int,
		et0::Float64, 
		etf::Float64,
		canonical_params=nothing,
		frame::String="ECLIPJ2000"
	)

Create function with signature:
`locate_spice(epoch) -> r::Vector{Float64}, v::Vector{Float64}`
"""
function get_spice_locate_function(
	body_id::Int,
	et0::Float64, 
	canonical_params=nothing,
	frame::String="ECLIPJ2000"
)
	# if canonical parameters is not provided, create Sun-Earth based ones
	if isnothing(canonical_params)
		# dynamics constants
		MU_SUN = 1.3271244004193938e11
		mu =  1.0
		lstar = 1.495978707e8  # 1AU in km
		vstar = sqrt(MU_SUN / lstar)
		tstar = lstar / vstar
		canonical_params = SunEarthCanonical(
			lstar, vstar, tstar, et0
		)
	end

	# construct functions to locate positions
	function locate_spice(epoch::Float64)
		sv = spkssb(body_id, et0 + epoch*canonical_params.tstar, frame)
		return sv[1:3]/lstar, sv[4:6]/canonical_params.vstar
	end
	return locate_spice, canonical_params
end


"""
	get_spice_locate_function(
		body_id::Int,
		et0_str::String="2022-01-01T00:00:00.00", 
		canonical_params=nothing,
		frame::String="ECLIPJ2000"
	)

Alias with string-based inputs for earliest epoch
"""
function get_spice_locate_function(
	body_id::Int,
	et0_str::String="2022-01-01T00:00:00.00", 
	canonical_params=nothing,
	frame::String="ECLIPJ2000"
)
	# convert epochs
	et0 = utc2et(et0_str)
	return get_spice_locate_function(
		body_id, et0, canonical_params, frame
	)
end
