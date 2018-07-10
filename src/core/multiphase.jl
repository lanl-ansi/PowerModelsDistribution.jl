function wrapto180(degrees)
    return degrees - 360*floor((degrees + 180)/360)
end

function wraptopi(radians)
    return radians - 2*pi*floor((radians+pi)/(2*pi))
end
