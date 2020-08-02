module FlowGeometry

#   This file is part of FlowGeometry.jl. It is licensed under the AGPL license
#   FlowGeometry Copyright (C) 2020 Michael Reed

using Requires, StaticArrays, Grassmann, Adapode

import Grassmann: points, edges, initedges
import Base: @pure

export Joukowski, NACA, joukowski, naca

include("profiles.jl")
include("airfoils.jl")

chord(N::A) where A<:Profile{p} where p = Grassmann.initpoints(interval(N))
chordedges(N::A) where A<:Profile{P} where P = Chain{ChainBundle(chord(N)),1}.(1:P-1,2:P)

function edgeslist(p::ChainBundle)
    Chain{p(2,3),1}.(1:length(p),[2:length(p);1])
end
function edgeslist!(p,r)
    l,n = length(p),length(r)
    push!(value(p),r...)
    Chain{p(2,3),1}.(l+1:l+n,l.+[2:n;1])
end

addbound(e,r=rectangle(xn,xm,yn,ym)) = [e;edgeslist!(points(e),r)]
airfoilbox(n) = addbound(airfoiledges(n),rectangle(xn,xm,yn,ym))
airfoiledges(n) = edgeslist(ChainBundle(points(n)))

function rectcirc(n,xn,xm,yn,ym,c::Chain{V}=Chain{SubManifold(ℝ^3),1}(1.0,0.0,0.0)) where V
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

# Rakich stretch mesh

export triangle, triangles, initrakich, arc, arcslope, FittedPoint

triangle(i,m,p) = triangle(((i-1)%(2(m-1)))+1,((i-1)÷(2(m-1)))+1,m,p)
triangle(i,j,m,p) = (k=(j-1)*m+(i÷2)+1;n=m+k;Chain{p,1}(isodd(i) ? SVector(k,n,k+1) : SVector(k,n-1,n)))
#triangle(i,j,m,p) = (k=(j-1)*m+(i÷2)+1;n=m+k;Chain{p,1}(iseven(i) ? SVector(k,n,n+1) : SVector(k,k-1,n)))

triangles(p,m=51,JL=51) = [triangle(i,JL,p) for i ∈ 1:2*(m-1)*(JL-1)]

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

function RakichPlate(x=0,c=1,n=21,D=50,JL=51)
    Δx = (c-x)/(n-1)
    r = RakichLine(x,D,((JL-n)÷2)+1,Δx)
    [-reverse(r);(x:Δx:c)[2:end-1];r.+1]
end

function RakichPoint(k,x,y,s,κ,JL=legnth(y))
    xk,yk = ((k-1)÷JL)+1,((k-1)%JL)+1
    Chain{SubManifold(ℝ^3),1}(1.0,x[xk],s[xk]≠0 ? Rakich(κ[xk],yk,s[xk],y[end],JL) : y[yk])
end

function RakichPoints(x0=0,c=1,t=0.06,D=50,n=51,m=21,JL=51)
    x = RakichPlate(x0,c,m,D,n)
    y = RakichLine(0,D,JL,t/10)
    s = arc.(x,t,c,x0)
    κ = [k≠0 ? RakichNewton(D-k,JL,t/10) : 0.0 for k ∈ s]
    [RakichPoint(k,x,y,s,κ,JL) for k ∈ 1:n*JL]
end

FittedPoint(k,JL=51) = Chain{SubManifold(ℝ^3),1}(1.0,(k-1)÷JL,(k-1)%JL)

RectangleBounds(n=51,JL=51) = [1:JL:JL*n; JL*(n-1)+2:JL*n; JL*(n-1):-JL:JL; JL-1:-1:2]

function initrakich(x0=0,c=1,t=0.06,D=50,n=101,m=61,JL=51)
    p = ChainBundle(RakichPoints(x0,c,t,D,n,m,JL))
    b,V = RectangleBounds(n,JL),p(2,3); push!(b,1)
    e = ChainBundle(Chain{V,1,Int}.([(b[i],b[i+1]) for i ∈ 1:length(b)-1]))
    return p,e,ChainBundle(triangles(p,n))
end

# points

export icosahedron, cube, box, sphere

square(x) = square(-x,x)
square(xn,xm) = rectangle(xn,xm,xn,xm)
rectangle(xn,xm,yn,ym) = Chain{SubManifold(ℝ^3),1}.(SVector{3,Float64}.(
    [(1.0,xn, yn), (1.0,xm, yn),
     (1.0,xm, ym), (1.0,xn, ym)]))

cube(x) = cube(-x,x)
cube(xn,xm) = box(xn,xm,xn,xm,xn,xm)
box(xn,xm,yn,ym,zn,zm) = Chain{SubManifold(ℝ^4),1}.(SVector{4,Float64}.(
    [(1.0,xn, yn, zn), (1.0,xm, yn, zn),
     (1.0,xm, ym, zn), (1.0,xn, ym, zn),
     (1.0,xn, yn, zm), (1.0,xm, yn, zm),
     (1.0,xm, ym, zm), (1.0,xn, ym, zm)]))

icosahedron(a=1,b=a*Irrational{:φ}()) =
    Chain{SubManifold(ℝ^4),1}.(SVector.(
     [(1,0,a,b),(1,b,0,a),(1,a,b,0),
      (1,0,a,-b),(1,-b,0,a),(1,a,-b,0),
      (1,0,-a,b),(1,b,0,-a),(1,-a,b,0),
      (1,0,-a,-b),(1,-b,0,-a),(1,-a,-b,0)]))

@generated function circlemid(x,r)
    v = Λ(SubManifold(ℝ^4)).v1
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

function cubesphere(r=1,c=2)
    p = ChainBundle(sphere(r))
    t = sphere(sphere(∂(MiniQhull.delaunay(p)),r),r)
    push!(value(p),cube(c*r)...)
    [t;∂(MiniQhull.delaunay(p,length(p)-7:length(p)))]
end

# init

function __init__()
    @require AbstractPlotting="537997a7-5e4e-5d89-9595-2241ea00577e" begin
        function AbstractPlotting.plot(N::Profile,args...)
            AbstractPlotting.plot(interval(N),profile(N),args...)
        end
        function AbstractPlotting.plot!(N::Profile,args...)
            AbstractPlotting.plot!(interval(N),profile(N),args...)
        end
        function AbstractPlotting.plot(N::Airfoil,args...)
            U,L = upperlower(N)
            AbstractPlotting.plot(real.(U),imag.(U),args...)
            length(U)==length(L) && AbstractPlotting.plot!(real.(U),(imag.(U).+imag.(L))./2,args...)
            AbstractPlotting.plot!(real.(L),imag.(L),args...)
        end
        function AbstractPlotting.plot!(N::Airfoil,args...)
            U,L = upperlower(N)
            AbstractPlotting.plot!(real.(U),imag.(U),args...)
            length(U)==length(L) && AbstractPlotting.plot!(real.(U),(imag.(U).+imag.(L))./2,args...)
            AbstractPlotting.plot!(real.(L),imag.(L),args...)
        end
        function AbstractPlotting.plot(N::NACA,args...)
            U,L = upperlower(N)
            AbstractPlotting.plot(real.(U),imag.(U),args...)
            AbstractPlotting.plot!(real.(L),imag.(L),args...)
            AbstractPlotting.plot!(N.c)
            #AbstractPlotting.plot!(N.t)
        end
        function AbstractPlotting.plot!(N::NACA,args...)
            U,L = upperlower(N)
            AbstractPlotting.plot!(real.(U),imag.(U),args...)
            AbstractPlotting.plot!(real.(L),imag.(L),args...)
            AbstractPlotting.plot!(N.c)
            #AbstractPlotting.plot!(N.t)
        end
        function AbstractPlotting.lines(N::Airfoil,args...)
            P = complex(N)
            AbstractPlotting.lines([real.(P) imag.(P)],args...)
        end
        function AbstractPlotting.lines!(N::Airfoil,args...)
            P = complex(N)
            AbstractPlotting.lines!([real.(P) imag.(P)],args...)
        end
    end
    @require TetGen="c5d3f3f7-f850-59f6-8a2e-ffc6dc1317ea" begin
        spheremesh() = TetGen.tetrahedralize(cubesphere(),"vpq1.414a0.1";holes=[Point(0.0,0.0,0.0)])
        cubemesh() = TetGen.tetrahedralize(∂(MiniQhull.delaunay(cube(2))), "vpq1.414a0.1")
    end
    @require MiniQhull="978d7f02-9e05-4691-894f-ae31a51d76ca" begin end
    @require MATLAB="10e44e05-a98a-55b3-a45b-ba969058deb6" begin
        decsg(args...) = MATLAB.mxcall(:decsg,1,args...)
        function decsg(N::Airfoil{pts}) where pts
            P = complex(N)[1:end-1]
            R = [[3,4,-1,2,2,-1,1,1,-1,-1];zeros(4pts-12)]
            A = [[2,2pts-2];[real.(P);imag.(P)]]
            decsg([R A],"R-A","RA")
        end
        export decsg
    end
end

end # module
