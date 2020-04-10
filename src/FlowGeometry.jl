module FlowGeometry

using Requires, StaticArrays
import Grassmann: points

export Joukowski, NACA, NACA5, joukowski, naca, naca5, nacasplit, nacapoints

abstract type Airfoil{P} end

struct Joukowski{R,f,g,b,p} <: Airfoil{p}
    Base.@pure Joukowski{R,f,g,b,p}() where {R,f,g,b,p} = new{R,f,g,b,p}()
end

joukowski(R,f,g,b,p) = Joukowski{R,f,g,b,p}()

struct NACA{N,P} <: Airfoil{P}
    Base.@pure NACA{N,P}() where {N,P} = new{N,P}()
end

naca(n,p) = NACA{n,p}()

function nacasplit(::NACA{n,p}) where {n,p}
    ns=string(n,pad=4)
    Mi,Pi,Ti = parse.(Int,(ns[1],ns[2],ns[3:4]))
    M,P,T = Mi/100,Pi/10,Ti/100
    x = range(0,1,length=p)
    yc = zeros(p) # camber
    for i ∈ 1:p
        if x[i] ≥ 0 && x[i] < P
            yc[i] = (M/P^2)*((2P*x[i])-x[i]^2)
        elseif x[i] ≥ P && x[i] ≤ 1
            yc[i] = (M/(1-P)^2)*(1-2P+2P*x[i]-x[i]^2)
        end
    end
    nacasplit(x,yc,P,M,T)
end

struct NACA5{N,P} <: Airfoil{P}
    Base.@pure NACA5{N,P}() where {N,P} = new{N,P}()
end

naca5(n,p) = NACA5{n,p}()

function nacasplit(::NACA5{n,p}) where {n,p}
    ns = string(n,pad=5)
    Ci,Pi,Ti = parse.(Int,(ns[1],ns[2:3],ns[4:5]))
    Cl,P,T = 3Ci/20,Pi/20,Ti/100
    x = range(0,1,length=p)
    yc = zeros(p) # camber
    for i ∈ 1:pts
        if x[i] ≥ 0 && x[i] < P
            yc[i] = k1*(x[i]^3-3M*x[i]^2+M^2*(3-M)*x[i])/6
        elseif x[i] ≥ P && x[i] ≤ 1
            yc[i] = k1*M^3*(1-x[i])/6
        end
    end
    nacasplit(x,yc,M,P,T)
end

function points(::Joukowski{R,f,g,b,p}) where {R,f,g,b,p}
    θ = range(0,2π,length=2p-2)
    z = R*cis.(θ).-(f-g*im)
    z.+b^2*conj.(z)./abs2.(z)
end

function points(N::Union{NACA,NACA5})
    U,L = nacasplit(N)
    [U;reverse(L)[2:end]]
end

function nacasplit(x,yc,P,M,T)
    p = length(x)
    dyc_dx = zeros(p) # gradient
    for i ∈ 1:p
        if x[i] ≥ 0 && x[i] < P
            dyc_dx[i] = ((2M)/P^2)*(P-x[i])
        elseif x[i] ≥ P && x[i] ≤ 1
            dyc_dx[i] = (2M/(1-P)^2)*(P-x[i])
        end
    end
    θ = atan.(dyc_dx)
    # thickness distribution
    a = SVector(0.2969,-0.1260,-0.3515,0.2843,-0.1036)
    yt = 5T*(a[1]*sqrt.(x).+a[2]*x+a[3]*x.^2+a[4]*x.^3+a[5]*x.^4)
    # upper surface points xu+im*yu
    U = (x.+im.*yc).+im.*cis.(θ).*yt; U[end] = 1
    # lower surface points xl+im*yl
    L = (x.+im.*yc).-im.*cis.(θ).*yt; L[end] = 1
    return U,L
end

function __init__()
    @require AbstractPlotting="537997a7-5e4e-5d89-9595-2241ea00577e" begin
        function AbstractPlotting.plot(N::Airfoil,args...)
            U,L = nacasplit(N)
            AbstractPlotting.plot(real.(U),imag.(U),args...)
            AbstractPlotting.plot!(real.(L),imag.(L),args...)
        end
        function AbstractPlotting.plot!(N::Airfoil,args...)
            U,L = nacasplit(N)
            AbstractPlotting.plot!(real.(U),imag.(U),args...)
            AbstractPlotting.plot!(real.(L),imag.(L),args...)
        end
        function AbstractPlotting.lines(N::Airfoil,args...)
            P = points(N)
            AbstractPlotting.lines([real.(P) imag.(P)],args...)
        end
        function AbstractPlotting.lines!(N::Airfoil,args...)
            P = points(N)
            AbstractPlotting.lines!([real.(P) imag.(P)],args...)
        end
    end
    @require MATLAB="10e44e05-a98a-55b3-a45b-ba969058deb6" begin
        decsg(args...) = MATLAB.mxcall(:decsg,1,args...)
        function decsg(N::Airfoil{pts}) where pts
            P = points(N)
            R = [[3,4,-1,2,2,-1,1,1,-1,-1];zeros(4pts-12)]
            A = [[2,2pts-2];[real.(P[1:end-1]);imag.(P[1:end-1])]]
            decsg([R A],"R-A","RA")
        end
        export decsg
    end
end

end # module
