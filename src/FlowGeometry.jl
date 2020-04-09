module FlowGeometry

using Requires, StaticArrays

#=using Compose, Grassmann
import Grassmann: vector

@inline field_components(M,n=2) = ([[M[x,y][k] for x ∈ 1:length(M[:,1]), y ∈ 1:length(M[1,:])] for k ∈ 1:n]...,)

function __init__()
    @require PyPlot="d330b81b-6aea-500a-939a-2ce795aea3ee" begin
        export streamplot
        function PyPlot.streamplot(t::T,x=collect(-5:0.1:5),y=x) where T<:TensorAlgebra{V} where V
            PyPlot.streamplot(x,y,field_components([vector(t*(x*Λ(V).v1 - y*Λ(V).v2)) for x ∈ x, y ∈ y])...)
        end
    end
end=#

export NACA, naca, nacapoints

struct NACA{N,P} end

naca(n,p) = NACA{n,p}()

function nacapoints(::NACA{n,pts}) where {n,pts}
    ns=string(n,pad=4)
    Mi,Pi,Ti = parse.(Int,(ns[1],ns[2],ns[3:4]))
    a = SVector(0.2969,-0.1260,-0.3515,0.2843,-0.1036)
    M,P,T = Mi/100,Pi/10,Ti/100
    x = range(0,1,length=pts)
    # camber and gradient
    yc,dyc_dx = zeros(pts),zeros(pts)
    for i ∈ 1:pts
        if x[i] ≥ 0 && x[i] < P
            yc[i] = (M/P^2)*((2P*x[i])-x[i]^2)
            dyc_dx[i] = ((2M)/P^2)*(P-x[i])
        elseif x[i] ≥ P && x[i] ≤ 1
            yc[i] = (M/(1-P)^2)*(1-2P+2P*x[i]-x[i]^2)
            dyc_dx[i] = (2M/(1-P)^2)*(P-x[i])
        end
    end
    θ = atan.(dyc_dx)
    # thickness distribution
    yt = 5T*(a[1]*sqrt.(x).+a[2]*x+a[3]*x.^2+a[4]*x.^3+a[5]*x.^4)
    # upper surface points xu+im*yu
    U = (x.+im.*yc).+im.*cis.(θ).*yt; U[end] = 1
    # lower surface points xl+im*yl
    L = (x.+im.*yc).-im.*cis.(θ).*yt; L[end] = 1
    return U,L
end

function nacageom(N::NACA{n,pts}) where {n,pts}
    U,L = nacapoints(N)
    rect = [[3,4,-1,2,2,-1,1,1,-1,-1];zeros(4pts-12)]
    [rect [[2,2pts-2];real.(U);reverse(real.(L)[2:end-1]);imag.(U);reverse(imag.(L)[2:end-1])]]
end

function __init__()
    @require MATLAB="10e44e05-a98a-55b3-a45b-ba969058deb6" begin
        decsg(args...) = MATLAB.mxcall(:decsg,1,args...)
        decsg(N::NACA) = decsg(nacageom(N),"R-A","RA")
        export decsg
    end
end

end # module
