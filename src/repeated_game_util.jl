#=
File contains a variety of utilities that make solving repeated games
easier
=#

"""
Places `npts` equally spaced points along the
2 dimensional circle and returns the points
with x coordinates in first column and y
coordinates in second column

i.e. if you wanted point i, it would be pts[i, :]
"""
function unitcircle(npts::Int)
    # Want our points placed on [0, 2π]
    incr = 2π / npts
    degrees = 0.0:incr:2π

    # Points on circle
    pts = Array(Float64, npts, 2)
    for i=1:npts
        x = degrees[i]
        pts[i, 1] = cos(x)
        pts[i, 2] = sin(x)
    end
    return pts
end

