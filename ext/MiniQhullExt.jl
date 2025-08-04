module MiniQhullExt

#   This file is part of FlowGeometry.jl
#   It is licensed under the AGPL license
#   FlowGeometry Copyright (C) 2020 Michael Reed
#       _           _                         _
#      | |         | |                       | |
#   ___| |__   __ _| | ___ __ __ ___   ____ _| | __ _
#  / __| '_ \ / _` | |/ / '__/ _` \ \ / / _` | |/ _` |
# | (__| | | | (_| |   <| | | (_| |\ V / (_| | | (_| |
#  \___|_| |_|\__,_|_|\_\_|  \__,_| \_/ \__,_|_|\__,_|
#
#   https://github.com/chakravala
#   https://crucialflow.com

using Grassmann, Cartan, FlowGeometry
isdefined(FlowGeometry, :Requires) ? (import FlowGeometry: MiniQhull) : (using MiniQhull)

function FlowGeometry.spheresurf(r=1)
    p = PointCloud(sphere(r))
    t = sphere(sphere(∂(MiniQhull.delaunay(p)),r),r)
    p(t)
end
function FlowGeometry.cubesphere(r=1,c=2)
    p = PointCloud(sphere(r))
    t = sphere(sphere(∂(MiniQhull.delaunay(p)),r),r)
    push!(fullpoints(p),cube(c*r)...)
    np = length(p)
    out = [topology(t);topology(∂(MiniQhull.delaunay(p,np-7:np)))]
    p(SimplexTopology(0,out,np))
end

end # module
