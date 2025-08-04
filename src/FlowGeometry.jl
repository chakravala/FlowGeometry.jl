module FlowGeometry

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
#            ________)
#           (, /         /)
#             /___,     //   ____   _
#          ) /         (/_  (_) (_(/
#         (_/
#
#        _____)
#      /
#     /   ___     _  ______    _ _/_  __
#    /     / )  _(/_(_) // (__(/_(__ / (_  (_/_
#   (____ /                               .-/
#                                        (_/

using Grassmann, Cartan, LinearAlgebra

import Cartan: PointCloud, points, edges, initpoints, initedges
import Base: @pure

export chord, chordedges, edgeslist, edgeslist!, addbound, airfoilbox, airfoiledges, rectcirc

include("profiles.jl")
include("airfoils.jl")

initpoints(N::A) where A<:Profile{p} where p = initpoints(interval(N))
initedges(N::A) where A<:Profile{P} where P = initedges(PointCloud(initpoints(N)))

function edgeslist(p::PointCloud)
    Values{2,Int}.(1:length(p),[2:length(p);1])
end
function edgeslist!(p,r)
    l,n = length(p),length(r)
    push!(p,r...)
    Values{2,Int}.(l+1:l+n,l.+[2:n;1])
end

# need PointCloud etc operations for refinement (same ID) or union etc (new ID)

addbound(e,r=rectangle(xn,xm,yn,ym)) = [e;edgeslist!(points(e),r)]
airfoilbox(n) = addbound(airfoiledges(n),rectangle(xn,xm,yn,ym))
airfoiledges(n) = edgeslist(PointCloud(points(n)))

function rectcirc(n,xn,xm,yn,ym,c::Chain{V}=Chain{Submanifold(ℝ^3),1}(1.0,0.0,0.0)) where V
    x = float.([yn,xm,ym,xn])
    r = rectangle(xn,xm,yn,ym)
    rc = r.-c
    at = [atan(k[3]/k[2]) for k ∈ rc] .+ [π,2π,0,π]
    out = typeof(r)(undef,0)
    for i ∈ 1:4
        push!(out,r[i])
        if isodd(i)
            v = x[i]-c[3]
            for k ∈ 1:n-2
                w = v/tan(at[i]+k*((at[(i%4)+1]-at[i]+2π)%(2π))/(n-1))
                push!(out,Chain{V,1}(1.0,w,x[i]))
            end
        else
            v = x[i]-c[2]
            for k ∈ 1:n-2
                w = v*tan(at[i]+k*((at[(i%4)+1]-at[i]+2π)%(2π))/(n-1))
                push!(out,Chain{V,1}(1.0,x[i],w))
            end

        end
    end
    return out
end

# rectangular structured triangles

export rectangletriangle, rectangletriangles, FittedPoint

rectangletriangle(i,m) = rectangletriangle(((i-1)%(2(m-1)))+1,((i-1)÷(2(m-1)))+1,m)
function rectangletriangle(i,j,m)
    k=(j-1)*m+(i÷2)+1;n=m+k
    isodd(i) ? Values(k,n,k+1) : Values(k,n-1,n)
end
#rectangletriangle(i,j,m) = (k=(j-1)*m+(i÷2)+1;n=m+k; iseven(i) ? Values(k,n,n+1) : Values(k,k-1,n))

function rectangletriangles(m=51,JL=51)
    SimplexTopology([rectangletriangle(i,JL) for i ∈ 1:2*(m-1)*(JL-1)],m*JL)
end
function rectanglebounds(n=51,JL=51)
    b = [1:JL:JL*n; JL*(n-1)+2:JL*n; JL*(n-1):-JL:JL; JL-1:-1:2]; push!(b,1)
    SimplexTopology(Values{2,Int}.([(b[i],b[i+1]) for i ∈ 1:length(b)-1]),n*JL)
end

FittedPoint(k,JL=51) = Chain{Submanifold(ℝ^3),1}(1.0,(k-1)÷JL,(k-1)%JL)

# Rakich stretch mesh

export initrakich, arc, arcslope, FittedPoint

function RakichNewton(D=50,JL=51,Δy=6e-3)
    κ,j = 1,1/(JL-1)
    for n ∈ 1:10
        eκ,ejκ = ℯ^κ, ℯ^(j*κ)
        κ -= (((ejκ-1)*D - (eκ-1)*Δy)*(eκ-1)) / ((ejκ*(eκ*(j-1)-j) + eκ)*D)
    end
    return κ
end

Rakich(κ,j,y0=0,D=50,JL=51) = y0 + (D-y0)*(exp(κ*(j-1)/(JL-1))-1)/(exp(κ)-1)
RakichLine(y=0,D=50,JL=51,Δy=6e-3,κ=RakichNewton(D-y,JL,Δy)) = Rakich.(κ,1:JL,y,D,JL)

function RakichPlate(::Profile{n},D=50,JL=51) where n
    x = interval(n)
    r = RakichLine(-x[1],D,((JL-n)÷2)+1,Float64(x.step))
    [-reverse(r);x[2:end-1];r.+((x[1]>0 ? 2 : -1)*x[1]+x[end])]
end

function RakichPoint(k,x,y,s,κ,JL=legnth(y))
    xk,yk = ((k-1)÷JL)+1,((k-1)%JL)+1
    Chain{Submanifold(ℝ^3),1}(1.0,x[xk],s[xk]≠0 ? Rakich(κ[xk],yk,s[xk],y[end],JL) : y[yk])
end

function rakichpoints(P::CircularArc{T,m}=CircularArc{6,21}(),D=50,n=51,JL=51) where {T,m}
    t,c,x0 = T/100,1,0
    x = RakichPlate(P,D,n)
    y = RakichLine(0,D,JL,t/10)
    s = profile.(Ref(P),x,c,x0)
    κ = [k≠0 ? RakichNewton(D-k,JL,t/10) : 0.0 for k ∈ s]
    PointCloud([RakichPoint(k,x,y,s,κ,JL) for k ∈ 1:n*JL])
end

function initrakich(P=CircularArc{6,61}(),D=50,n=101,JL=51)
    p = rakichpoints(P,D,n,JL)
    p(rectangletriangles(n,JL)),p(rectanglebounds(n,JL))
end

# points

export icosahedron, cube, box, sphere

square(x) = square(-x,x)
square(xn,xm) = rectangle(xn,xm,xn,xm)
rectangle(xn,xm,yn,ym) = Chain{Submanifold(ℝ^3),1}.(Values{3,Float64}.(
    [(1.0,xn, yn), (1.0,xm, yn),
     (1.0,xm, ym), (1.0,xn, ym)]))

cube(x) = cube(-x,x)
cube(xn,xm) = box(xn,xm,xn,xm,xn,xm)
box(xn,xm,yn,ym,zn,zm) = Chain{Submanifold(ℝ^4),1}.(Values{4,Float64}.(
    [(1.0,xn, yn, zn), (1.0,xm, yn, zn),
     (1.0,xm, ym, zn), (1.0,xn, ym, zn),
     (1.0,xn, yn, zm), (1.0,xm, yn, zm),
     (1.0,xm, ym, zm), (1.0,xn, ym, zm)]))

icosahedron(a=1,b=a*Irrational{:φ}()) =
    Chain{Submanifold(ℝ^4),1}.(Values.(
     [(1,0,a,b),(1,b,0,a),(1,a,b,0),
      (1,0,a,-b),(1,-b,0,a),(1,a,-b,0),
      (1,0,-a,b),(1,b,0,-a),(1,-a,b,0),
      (1,0,-a,-b),(1,-b,0,-a),(1,-a,-b,0)]))

@generated function circlemid(x,r)
    v = Λ(Submanifold(ℝ^4)).v1
    :(r*Grassmann.unit(x/2-$v) + $v)
end

sphere(r=1) = icosahedron(r/sqrt(1+Irrational{:φ}()^2))
function sphere(fac::SimplexBundle,r=1,p=fullpoints(fac))
    M = Manifold(fac)
    out = Values{3,Int}[]
    for f ∈ topology(fac)
        p0,p1,p2 = p[f]
        p3 = circlemid(p0+p1,r)
        p4 = circlemid(p1+p2,r)
        p5 = circlemid(p2+p0,r)
        v0,v1,v2 = f
        v3,v4,v5 = length(p).+(1,2,3)
        push!(fullpoints(p),p3,p4,p5)
        push!(out,Values(v0, v3, v5))
        push!(out,Values(v3, v1, v4))
        push!(out,Values(v4, v2, v5))
        push!(out,Values(v3, v4, v5))
    end
    np = length(p)
    return fac(SimplexTopology(0,out,Base.OneTo(np),np))
end

function wing(N::Airfoil,λ=0.7,σ=0.5)
    U,L = imag.(upper(N)),imag.(lower(N))
    np = length(U)
    out1 = zeros(np,2np-1)
    out2 = zeros(np,2np-1)
    out3 = zeros(np,2np-1)
    int = interval(N)
    rint = reverse(int)
    pf = profile(N.c)
    for i ∈ 1:np-1
        out1[:,i] .= (σ*int[i]).+(λ+(1-λ)*rint[i]).*int
        out2[:,i] .= int[i]
        out3[:,i] .= fiber(rint[i]*U)
    end
    out1[:,np] .= (σ*int[np]).+(λ*rint[np]).*int
    out2[:,np] .= int[np]
    for i ∈ 1:np-1
        out1[:,i+np] .= (σ*rint[i+1]).+(λ+(1-λ)*int[i+1]).*int
        out2[:,i+np] .= rint[i+1]
        out3[:,i+np] .= fiber(int[i+1]*L)
    end
    TensorField(int⊕(-1:step(int):1),Chain.(out1,out2,out3))
end

function convhull(p::PointCloud)
    n = length(p)
    out = Values{2,Int}[]
    for i ∈ 1:n
        pi = base(p)[i]
        for j ∈ 1:n
            i==j && continue
            pj = base(p)[j]
            pij = pi∧pj
            val = true
            for k ∈ 1:n
                k∈(i,j) && continue
                pk = base(p)[k]
                pijk = pij∧pk
                if Real(pijk) ≈ 0
                    if Real(abs((pi+pj)/2-pk)) < Real(abs(pi-pj))
                        val = false
                        break
                    end
                elseif Real(pijk) > 0
                    val = false
                    break
                end
            end
            val && push!(out,Values(i,j))
        end
    end
    return SimplexTopology(out,n)
end

function convhull(p::PointCloud,r)
    n = length(p)
    out = Values{2,Int}[]
    for i ∈ 1:n
        pi = base(p)[i]
        for j ∈ 1:n
            i==j && continue
            pj = base(p)[j]
            pij = pi∧pj
            avg = (pi+pj)/2
            val = true
            for k ∈ 1:n
                k∈(i,j) && continue
                pk = base(p)[k]
                dist = Real(abs(avg-pk))
                dist > r && continue
                pijk = pij∧pk
                if Real(pijk) ≈ 0
                    if dist < Real(abs(pi-pj))
                        val = false
                        break
                    end
                elseif Real(pijk) > 0
                    val = false
                    break
                end
            end
            val && push!(out,Values(i,j))
        end
    end
    return SimplexTopology(out,n)
end

#if !isdefined(Base, :get_extension)
using Requires
#end

#@static if !isdefined(Base, :get_extension)
function __init__()
    @require Makie="ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a" include("../ext/MakieExt.jl")
    @require UnicodePlots="b8865327-cd53-5732-bb35-84acbb429228" include("../ext/UnicodePlotsExt.jl")
    @require TetGen="c5d3f3f7-f850-59f6-8a2e-ffc6dc1317ea" include("../ext/TetGenExt.jl")
    @require MiniQhull="978d7f02-9e05-4691-894f-ae31a51d76ca" include("../ext/MiniQhullExt.jl")
    @require MATLAB="10e44e05-a98a-55b3-a45b-ba969058deb6" include("../ext/MATLABExt.jl")
end
#end

end # module
