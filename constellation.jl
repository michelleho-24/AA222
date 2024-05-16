using SatelliteToolbox
using SateliteDynamics
using Plots 



function surface_area(h,re,cone_angle)
    # h is altitude above surface
    # re is radius of planet
    # cone_angle is angle associated with satellite (assigned)
    
    # Calculate angle at the center of planet
    beta = rad2deg(atan(h*tan(deg2rad(cone_angle))/re))
    Ae = 4*pi*re^2 # surface area of planet

    x = re*(1-cos(deg2rad(beta)))
    A = 2*pi*re*x # area covered by satellite 
    percent_A = A/Ae
    return A, percent_A
end 

function ground_track(lat::Vector,lon::Vector)
    # lat is nx1 array of coordinates in degrees
    # lon is n x 1 array of coordinates in degrees 
    lon = lon .+ 360*(lon .< 0)
    lon = lon .- 360*(lon .> 180)
    plot(lon, lat, marker = :circle, markersize = 3, color = :red, linestyle = :solid)

    plot!(xlims = (-180,180), ylims = (-90,90), 
        xticks = [-180, -120, -60, 0, 60, 120, 180], 
        yticks = [-90, -60, -30, 0, 30, 60, 90],
        aspect_ratio = :equal)

    scatter!([lon[1]],[lat[1]], color = :green, marker = :circle, markersize = 5, label = "Start")
    scatter!([lon[end]], [lat[end]], color = :orange, marker = :circle, markersize = 5, label = "End")
    im = load("Blue_Marble_2002.png")
    plot!(im, xlims = (-180, 180), ylims = (-90, 90), legend = false)

    xlabel!("Longitude [Deg]")
    ylabel!("Latitude [Deg]")
    savefig("ground_track.png")
end 