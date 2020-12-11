module MotifRecogntion

include("CollisionMatrix.jl")


"""
    error_dist(x1, x2)

Returns the number of differences between two sorted patterns 'x1' and 'x2'.
"""
function error_dist(x1, x2)
    return sum(x1 .!= x2)
end

end # module
