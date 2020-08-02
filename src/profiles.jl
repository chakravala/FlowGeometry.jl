
#   This file is part of FlowGeometry.jl. It is licensed under the AGPL license
#   FlowGeometry Copyright (C) 2020 Michael Reed

approx(x,y) = [x^i for i ∈ 0:length(y)-1]⋅y

export FlatPlate, ParabolicArc, CircularArc, ClarkY, Thickness, Modified, Camber
export profile, profileslope, interval

# profiles

abstract type Profile{P} end

@pure profile(p::T,x,c,x0=0) where T<:Profile = profile(p,(x-x0)/c)*c
@pure profileslope(p::T,x,c,x0=0) where T<:Profile = profileslope(p,(x-x0)/c)
profile(p::T) where T<:Profile = profile.(Ref(p),interval(p))
profileslope(p::T) where T<:Profile = profileslope.(Ref(p),interval(p))
@pure interval(::P,c=1,x0=0) where P<:Profile{p} where p = interval(p,c,x0)
@pure interval(p::Int,c=1,x0=0) = range(x0,c,length=p)
function points(N::A,c=1,x0=0) where A<:Profile{p} where p
    Chain{SubManifold(ℝ^3),1}.(1.0,interval(N,c,x0),profile(N))
end

# flat plate

struct FlatPlate{p} <: Profile{p}
    @pure FlatPlate{p}() where p = new{p}()
end

@pure profile(::FlatPlate,x,c=nothing,x0=nothing) = 0.0
@pure profileslope(::FlatPlate,x,c=nothing,x0=nothing) = 0.0
@pure profile(::FlatPlate{p}) where p = range(0,0,length=p)
@pure profileslope(f::FlatPlate) = profile(f)

# parabolic arc

struct ParabolicArc{t,p} <: Profile{p}
    @pure ParabolicArc{t,p}() where {t,p} = new{t,p}()
    @pure ParabolicArc{p}() where p = new{6,p}()
end

@pure profile(::ParabolicArc{t},x) where t = (x<0 || x>1) ? 0.0 : (t/25)*x*(1-x)
@pure profileslope(::ParabolicArc{t},x) where t = (x<0 || x>1) ? 0.0 : (1-2x)*t/400

# circular arc

struct CircularArc{t,p} <: Profile{p}
    @pure CircularArc{t,p}() where {t,p} = new{t,p}()
    @pure CircularArc{p}() where p = new{6,p}()
end

@pure function profile(::CircularArc{T},x) where T
    (x<0 || x>1) && return 0.0
    t = T/100; r = (t+1/4t)/2
    r*(sin(acos((x-1/2)/r))-1)+t
    #ct = -4t^2
    #((sqrt(-16(2x-1)^2*t^2 + ct^2) - ct)*ct) / (8abs(ct)*t)
end
@pure function profileslope(::CircularArc{T},x) where T
    (x<0 || x>1) && return 0.0
    #r = ((1/2)^2+t^2)/2t # Algebra.df(arc(:x),:x)
    #((1-2x)*r) / (sqrt(4r^2-((2x - 1)^2))*abs(r))
    t = T/100; t4 = 4t; xc = t4*(2x-1)
    -xc / sqrt((1 + t*t4)^2 - xc^2)
end
profileslope(P::CircularArc,x::Chain) = profileslope(P,x[2])
profileslope(P::CircularArc,x::Chain,c,x0=0) = profileslope(P,x[2],c,x0)

# NACA 4-series thickness

struct ClarkY{t,p} <: Profile{p}
    @pure ClarkY{t,p}() where {t,p} = new{t,p}()
end

@pure Y(x) = SVector(sqrt(x),x,x^2,x^3,x^4)
@pure dY(x) = SVector(0.5/sqrt(x),1.0,2x,3x^2,4x^3)
const _clarky = SVector(0.2969,-0.1260,-0.3515,0.2843,-0.1036)

@pure function thickness(x,a,t)
    (x<0 || x>1) ? 0.0 : 5t*(Y(x)⋅a)
end # thickness distribution
@pure function thickness_slope(x,a,t)
    (x<0 || x>1) ? 0.0 : 5t*(dY(x)⋅a)
end # thickness distribution slope

@pure profile(::ClarkY{t},x) where t = thickness(x,_clarky,t/100)
@pure profileslope(::ClarkY{t},x) where t = thickness_slope(x,_clarky,t/100)

# NACA custom thickness

struct Thickness{t,x,te,p} <: Profile{p}
    @pure Thickness{t,x,te,p}() where {t,x,te,p} = new{t,x,te,p}()
    @pure Thickness{t,x,p}() where {t,x,p} = Thickness{t,x,0.2,p}()
    @pure Thickness{t,p}() where {t,p} = Thickness{t,3,p}()
end

#[0.2,0.3,0.4,0.5,0.6] -> -[0.2,0.234,0.315,0.465,0.7]
@pure tailslope(x) = -approx(x,[31/200,151/300,-109/40,43/6,-5/2])
@pure riegel(x) = -(2.24-5.42x+12.3x^2)/(10(1-0.878x))

@pure function clarky(x=0.3,te=0.002,ts=tailslope(x),n=0.1,y=0.078)
    SMatrix{5,5}(Y(x)...,dY(x)...,Y(1)...,dY(1)...,Y(n)...)'\SVector(0.1,0,te,ts,y)
end

@pure @generated profile(::Thickness{t,m,te},x) where {t,m,te} = :(thickness(x,$(clarky(m/10,te/100)),$(t/100)))
@pure @generated profileslope(::Thickness{t,m,te},x) where {t,m,te} = :(thickness_slope(x,$(clarky(m/10,te/100)),$(t/100)))

# NACA mofified thickness

struct Modified{t,m,p} <: Profile{p}
    @pure Modified{t,m,p}() where {t,m,p} = new{t,m,p}()
    @pure Modified{t,p}() where {t,p} = new{t,63,p}()
end

@pure function decode(::Modified{t,m}) where {t,m}
    i = Int(m÷10)
    t/100,i≠9 ? i : 6sqrt(3),(m%10)/10
end

@pure D(x) = (x1 = 1-x; SVector(1.0,x1,x1^2,x1^3))
@pure dD(x) =  (x1 = x-1; SVector(0.0,-1.0,2x1,-3x1^2))
@pure ddD(x) =  SVector(0.0,0.0,2,-6(x-1))

@pure A(x) = SVector(sqrt(x),x,x^2,x^3)
@pure dA(x) =  SVector(0.5/sqrt(x),1.0,2x,3x^2)
@pure ddA(x) =  SVector(-1/(4x*sqrt(x)),0.0,2,6x)

@pure function modified(x,a,p,t=1)
    (x<0 || x>1) ? 0.0 : 5t*(x < p ? A(x)⋅a[1] : D(x)⋅a[2])
end # modified thickness
@pure function modified_slope(x,a,p,t=1)
    (x<0 || x>1) ? 0.0 : 5t*(x < p ? dA(x)⋅a[1] : dD(x)⋅a[2])
end # thickness distribution slope

@pure function profile(P::Modified,x)
    t,i,p = decode(P)
    modified(x,coef(t,i,p),p,t)
end
@pure function profileslope(P::Modified,x)
    t,i,p = decode(P)
    modified_slope(x,coef(t,i,p),p,t)
end

@pure radius(t,::Val{a0}=Val(_clarky[1])) where a0 = t^2*((a0/0.2)^2/2)

#transpose(Chain{V,1}(Chain{V,1}(D(x)),Chain{V,1}(dD(x)),Chain{V,1}(D(1)),Chain{V,1}(dD(1))))\Chain{V,1}(0.1,0.0,0.002,-tailangle(x))
#=function coef(t,i,x,te=0.002,ts=tailslope(x))
    di = SMatrix{4,4}(D(x)...,dD(x)...,D(1)...,dD(1)...)'\SVector(0.1,0,te,ts)
    ai = SMatrix{4,4}(A(x)...,dA(x)...,ddA(x)...,1,0,0,0)'\SVector(0.1,0,ddD(x)⋅di,(1/5t)*sqrt(2radius(0.1*i/6)))
    ai,di
end=#
@pure function coef(t,i,x,te=0.002,ts=tailslope(x))
    d0,d1,x1 = te,-ts,1-x
    d3 = (3d1 - (0.588/x1)) / 3x1^2
    d2 = (-1.5x1*d3) - (0.5d1 / x1)
    a0 = (0.2/t)*sqrt(2radius(t*i/6))
    a3 = 0.1/x^3 + (d1*x1-0.294)/(x*x1^2) - (3/8)*a0/(x^2.5)
    a2 = a0/2sqrt(x^3) - 0.1/x^2 - 2a3*x
    a1 = -(a0/2sqrt(x) + 2a2*x + 3a3*x^2)
    return SVector(a0,a1,a2,a3),SVector(d0,d1,d2,d3)
end

# NACA camber

struct Camber{n,p} <: Profile{p}
    @pure Camber{n,p}() where {n,p} = new{n,p}()
end

@pure function camber2(n) # decode
    ns = string(n,pad=2)
    M,P = parse.(Int,(ns[1],ns[2]))
    return M/100,P/10
end
@pure function camber3(n) # decode
    ns=string(n)
    C,P,R = parse.(Int,(ns[1],ns[2],ns[3]))
    m,p = C/2,P/20 # Cl = 3C/20
    r,k1,k1k2 = if !Bool(R)
        #p: [0.05,0.1,0.15,0.2,0.25]
        #m: p -> [0.058,0.126,0.2025,0.29,0.391]
        approx(p,[-1/250,359/300,7/10,10/3,0]),
        #k: p -> [361.4,51.64,15.957,6.643,3.230]
        approx(p,[284037/200,-3296847/100,4296829/15,-1087744,4544800/3]),0
    else
        #p: [0.1,0.15,0.2,0.25]
        #m: p -> [0.13,0.217,0.318,0.441]
        approx(p,[17/500,26/15,-2,32/3]),
        #k: p -> [51.99,15.793,6.52,3.191]
        approx(p,[72269/250,-583261/150,89864/5,83920/3]),
        #k/k: p-> [0.000764,0.00677,0.0303,0.1355]
        approx(p,[10763/50000,12081/25000,-87457/2500,10691/125])
    end
    return m*k1,r,k1k2
end

# derivation

#=yc(x) = SVector(1,x,x^2)
dyc(x) = SVector(0,1,2x)

@pure function clarkY(m,p)
    SMatrix{3,3}(yc(0)...,yc(p)...,dyc(p)...)'\SVector(0,m,0),
    SMatrix{3,3}(yc(1)...,yc(p)...,dyc(p)...)'\SVector(0,m,0)
end=#

@pure function camber2(::Camber{n},x) where n
    (x<0 || x>1) && return 0.0
    m,p = camber2(n)
    return m*(if x < p
        xp = x/p
        2xp-xp^2
    else
        (1+2p*(x-1)-x^2)/(1-p)^2
    end)
end
@pure function camber2slope(::Camber{n},x) where n
    (x<0 || x>1) && return 0.0
    m,p = camber2(n)
    2m*(p-x)/(x<p ? p^2 : (1-p)^2)
end
@pure function camber3(::Camber{n},x) where n
    (x<0 || x>1) && 0.0
    mk1,p,k1k2 = camber3(n)
    mk1*(p^3*(1-x)+(if x < p
        #x^3-3p*x^2+p^2*(3-p)*x
        (x-p)^3-(k1k2≠0 ? k1k2*x*(1-p)^3 : 0)
    else
        k1k2≠0 ? k1k2*((x-p)^3-x*(1-p)^3) : 0
    end))/6
end
@pure function camber3slope(::Camber{n},x) where n
    (x<0 || x>1) && 0.0
    mk1,p,k1k2 = camber3(n)
    mk1*(-p^3+(if x < p
        #3x^2-6p*x+p^2*(3-p)
        3(p-x)^2-(k1k2≠0 ? k1k2*(1-p)^3 : 0)
    else
        k1k2≠0 ? 3k1k2*((p-x)^2-(1-p)^3) : 0
    end))
end

@pure profile(P::Camber{n},x) where n = n < 100 ? camber2(P,x) : camber3(P,x)
@pure profileslope(P::Camber{n},x) where n = n < 100 ? camber2slope(P,x) : camber3slope(P,x)
profile(P::Camber{n}) where n = (x=interval(P); n < 100 ? camber2.(Ref(P),x) : camber3.(Ref(P),x))
profileslope(P::Camber{n}) where n = (x=interval(P); n < 100 ? camber2slope.(Ref(P),x) : camber3slope.(Ref(P),x))
