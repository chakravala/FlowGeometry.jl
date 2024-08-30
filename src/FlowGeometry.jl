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

using Requires, Grassmann, Cartan, LinearAlgebra

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

rectangletriangle(i,m,p) = rectangletriangle(((i-1)%(2(m-1)))+1,((i-1)÷(2(m-1)))+1,m,p)
function rectangletriangle(i,j,m,p)
    k=(j-1)*m+(i÷2)+1;n=m+k
    isodd(i) ? Values(k,n,k+1) : Values(k,n-1,n)
end
#rectangletriangle(i,j,m,p) = (k=(j-1)*m+(i÷2)+1;n=m+k;Chain{p,1}(iseven(i) ? Values(k,n,n+1) : Values(k,k-1,n)))

rectangletriangles(p,m=51,JL=51) = [rectangletriangle(i,JL,p) for i ∈ 1:2*(m-1)*(JL-1)]
rectanglebounds(n=51,JL=51) = [1:JL:JL*n; JL*(n-1)+2:JL*n; JL*(n-1):-JL:JL; JL-1:-1:2]

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

function RakichPoints(P::CircularArc{T,m}=CircularArc{6,21}(),D=50,n=51,JL=51) where {T,m}
    t,c,x0 = T/100,1,0
    x = RakichPlate(P,D,n)
    y = RakichLine(0,D,JL,t/10)
    s = profile.(Ref(P),x,c,x0)
    κ = [k≠0 ? RakichNewton(D-k,JL,t/10) : 0.0 for k ∈ s]
    [RakichPoint(k,x,y,s,κ,JL) for k ∈ 1:n*JL]
end

function initrakich(P=CircularArc{6,61}(),D=50,n=101,JL=51)
    p = PointCloud(RakichPoints(P,D,n,JL))
    b = rectanglebounds(n,JL); push!(b,1)
    e = SimplexManifold(Values{2,Int}.([(b[i],b[i+1]) for i ∈ 1:length(b)-1]),length(p))
    return p,e,SimplexManifold(rectangletriangles(p,n),length(p))
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
function sphere(fac::Vector,r=1,p=points(fac))
    M = Manifold(fac)
    out = typeof(fac[1])[]
    for f ∈ fac
        p0,p1,p2 = p[value(f)]
        p3 = circlemid(p0+p1,r)
        p4 = circlemid(p1+p2,r)
        p5 = circlemid(p2+p0,r)
        v0,v1,v2 = value(f)
        v3,v4,v5 = length(p).+(1,2,3)
        push!(value(p),p3,p4,p5)
        push!(out,Chain{M,1}(v0, v3, v5))
        push!(out,Chain{M,1}(v3, v1, v4))
        push!(out,Chain{M,1}(v4, v2, v5))
        push!(out,Chain{M,1}(v3, v4, v5))
    end
    return out
end

# init

function __init__()
    @require Makie="ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a" begin
        function Makie.plot(N::Profile,args...)
            Makie.plot(interval(N),profile(N),args...)
        end
        function Makie.plot!(N::Profile,args...)
            Makie.plot!(interval(N),profile(N),args...)
        end
        function Makie.plot(N::Airfoil,args...)
            Makie.lines(N,args...)
            Makie.plot!(N.c,args...)
        end
        function Makie.plot!(N::Airfoil,args...)
            Makie.lines!(N,args...)
            Makie.plot!(N.c,args...)
        end
        function Makie.lines(N::Airfoil,args...)
            P = complex(N)
            Makie.lines([real.(P) imag.(P)],args...)
        end
        function Makie.lines!(N::Airfoil,args...)
            P = complex(N)
            Makie.lines!([real.(P) imag.(P)],args...)
        end
        function Makie.plot(N::DoubleArc,args...)
            U,L = upperlower(N)
            Makie.plot(real.(U),imag.(U),args...)
            length(U)==length(L) && Makie.plot!(real.(U),(imag.(U).+imag.(L))./2,args...)
            Makie.plot!(real.(L),imag.(L),args...)
        end
        function Makie.plot!(N::DoubleArc,args...)
            U,L = upperlower(N)
            Makie.plot!(real.(U),imag.(U),args...)
            length(U)==length(L) && Makie.plot!(real.(U),(imag.(U).+imag.(L))./2,args...)
            Makie.plot!(real.(L),imag.(L),args...)
        end
    end
    @require UnicodePlots="b8865327-cd53-5732-bb35-84acbb429228" begin
        function UnicodePlots.Plot(N::Profile,args...)
            UnicodePlots.lineplot(interval(N),profile(N),args...)
        end
        function Plot!(p,N::Profile,args...)
            UnicodePlots.lineplot!(p,interval(N),profile(N),args...)
        end
        function UnicodePlots.Plot(N::Airfoil,args...)
            Plot!(UnicodePlots.lineplot(N,args...),N.c,args...)
        end
        function Plot!(p,N::Airfoil,args...)
            Plot!(UnicodePlots.lineplot!(p,N,args...),N.c,args...)
        end
        function UnicodePlots.lineplot(N::Airfoil,args...)
            P = complex(N)
            UnicodePlots.lineplot(real.(P),imag.(P),args...)
        end
        function UnicodePlots.lineplot!(p,N::Airfoil,args...)
            P = complex(N)
            UnicodePlots.lineplot!(p,real.(P),imag.(P),args...)
        end
        function UnicodePlots.Plot(N::DoubleArc,args...)
            U,L = upperlower(N)
            p = UnicodePlots.lineplot(real.(U),imag.(U),args...)
            length(U)==length(L) && UnicodePlots.lineplot!(p,real.(U),(imag.(U).+imag.(L))./2,args...)
            UnicodePlots.lineplot!(p,real.(L),imag.(L),args...)
        end
        function Plot!(p,N::DoubleArc,args...)
            U,L = upperlower(N)
            UnicodePlots.lineplot!(p,real.(U),imag.(U),args...)
            length(U)==length(L) && UnicodePlots.lineplot!(p,real.(U),(imag.(U).+imag.(L))./2,args...)
            UnicodePlots.lineplot!(p,real.(L),imag.(L),args...)
        end
    end
    @require TetGen="c5d3f3f7-f850-59f6-8a2e-ffc6dc1317ea" begin
        spheremesh() = TetGen.tetrahedralize(cubesphere(),"vpq1.414a0.1";holes=[Point(0.0,0.0,0.0)])
        cubemesh() = TetGen.tetrahedralize(∂(MiniQhull.delaunay(cube(2))), "vpq1.414a0.1")
    end
    @require MiniQhull="978d7f02-9e05-4691-894f-ae31a51d76ca" begin
        function cubesphere(r=1,c=2)
            p = PointCloud(sphere(r))
            t = sphere(sphere(∂(MiniQhull.delaunay(p)),r),r)
            push!(p,cube(c*r)...)
            [t;∂(MiniQhull.delaunay(p,length(p)-7:length(p)))]
        end
    end
    @require MATLAB="10e44e05-a98a-55b3-a45b-ba969058deb6" begin
        decsg(args...) = MATLAB.mxcall(:decsg,1,args...)
        function decsg(N::Airfoil{pts}) where pts
            P = complex(N)[1:end-1]
            x1,x2,x3,x4 = -1.5,3.5,3.5,-1.5
            y1,y2,y3,y4 = 1.5,1.5,-1.5,-1.5
            R = [[3,4,x1,x2,x3,x4,y1,y2,y3,y4];zeros(2pts-8)]
            A = [[2,pts];[real.(P);imag.(P)]]
            decsg([R A],"R-A","RA")
        end
        export decsg
    end
end

end # module
