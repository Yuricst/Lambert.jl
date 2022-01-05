"""
Utility expressions for numerical considerations
"""

function sqrt_vector(vec)
    return [sqrt(el) for el in vec]
end


function sin_vec(vec)
    return [sin(el) for el in vec]
end


function acos_safe(val)
    return acos(max(-1.0, min(1.0, val)))
end


function asin_safe(val)
    return asin(max(-1.0, min(1.0, val)))
end
