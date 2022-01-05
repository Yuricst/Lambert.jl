"""
Spice helpers
"""


"""
Create function with signature `locate_spice(epoch) -> r::Vector{Float64}, v::Vector{Float64}`
"""
function get_spice_locate_function()

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

	# scaled time-windows
	t0 = 0.0
	tf = (etf - et0)/canonical_params.tstar

	# construct functions to locate positions
	function locate_spice(epoch::Float64)
		sv = spkssb(bodies[1], et0 + epoch*canonical_params.tstar, "ECLIPJ2000")
		return sv[1:3]/lstar, sv[4:6]/canonical_params.vstar
	end

	return locate_spice
end