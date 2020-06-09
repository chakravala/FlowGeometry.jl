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
    θ = range(0,2π,length=2p-1)
    z = R*cis.(θ).-(f-g*im)
    z.+b^2*inv.(z)
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
    yt = thickness(x,T)
    # upper surface points xu+im*yu
    U = (x.+im.*yc).+im.*cis.(θ).*yt; U[end] = 1
    # lower surface points xl+im*yl
    L = (x.+im.*yc).-im.*cis.(θ).*yt; L[end] = 1
    return U,L
end

const Thickness = SVector(
    SVector(0.2969,-0.1260,-0.3515,0.2843,-0.1036))

function thickness(x,T,::Val{n}=Val(1)) where n
    5T*(a[n][1]*sqrt.(x).+a[n][2]*x+a[n][3]*x.^2+a[n][4]*x.^3+a[n][5]*x.^4)
end

# circular arc

export triangle, triangles, initrakich, arc, arcslope, FittedPoint

triangle(i,m,p) = triangle(((i-1)%(2(m-1)))+1,((i-1)÷(2(m-1)))+1,m,p)
triangle(i,j,m,p) = (k=(j-1)*m+(i÷2)+1;n=m+k;Chain{p,1}(isodd(i) ? SVector(k,k+1,n) : SVector(k,n,n-1)))
#triangle(i,j,m,p) = (k=(j-1)*m+(i÷2)+1;n=m+k;Chain{p,1}(isodd(i) ? SVector(k,n+1,n) : SVector(k-1,k,n)))

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

FittedPoint(k,JL=51) = Chain{SubManifold(ℝ^3),1}(1,(k-1)÷JL,(k-1)%JL)

function arc(x,t=0.06,c=1,x0=0)
    (x<x0 || x>c) && return 0.0
    r = ((c/2)^2+t^2)/2t
    r*sin(acos((x-x0-c/2)/r))-r+t
    #ct = c^2-4t^2
    #((sqrt(-16(2(x-x0)-c)^2*t^2 + ct^2) - ct)*ct) / (8abs(ct)*t)
end

arcslope(x::Chain,t=0.06,c=1,x0=0) = arcslope(x[2],t,c,x0)
function arcslope(x,t=0.06,c=1,x0=0)
    (x<x0 || x>c) && return 0.0
    #r = ((c/2)^2+t^2)/2t # Algebra.df(arc(:x),:x)
    #((c-2(x-x0))*r) / (sqrt(4r^2-((2(x-x0) - c)^2))*abs(r))
    xc,t2 = 2(x-x0)-c,t^2
    (-4xc*t2) / (sqrt((c^2 + 4t2)^2 - 16xc^2*t2)*t)
end

function initrakich(x0=0,c=1,t=0.06,D=50,n=101,m=61,JL=51)
    p = ChainBundle(RakichPoints(x0,c,t,D,n,m,JL))
    return p,ChainBundle(triangles(p,n))
end

# init

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
