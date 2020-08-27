
#   This file is part of FlowGeometry.jl. It is licensed under the AGPL license
#   FlowGeometry Copyright (C) 2020 Michael Reed

approx(x,y::SVector{N}) where N = poly(x,Val(N))⋅y
approx(x,y::AbstractVector) where N = [x^i for i ∈ 0:length(y)-1]⋅y

@generated poly(x,::Val{N}) where N = Expr(:call,:SVector,[:(x^$i) for i ∈ 0:N-1]...)

export FlatPlate, ParabolicArc, CircularArc, ClarkY, Thickness, Modified
export profile, profileslope, profileangle, interval
export NACA4, NACA5, NACA6, NACA6A

# profiles

abstract type Profile{P} end

profile(p::T) where T<:Profile = p.(interval(p))
profileslope(p::T) where T<:Profile = profileslope.(Ref(p),interval(p))
profileangle(p::T) where T<:Profile = atan.(profileslope(p))
@pure profile(p::T,x) where T<:Profile = p(x)
@pure profile(p::T,x,c,x0=0) where T<:Profile = profile(p,(x-x0)/c)*c
@pure profileslope(p::T,x,c,x0=0) where T<:Profile = profileslope(p,(x-x0)/c)
@pure profileangle(p::T,x,c,x0=0) where T<:Profile = atan(profileslope(p,x,c,x0))
@pure interval(::P,c=1,x0=0) where P<:Profile{p} where p = interval(p,c,x0)
@pure interval(p::Int,c=1,x0=0) = range(x0,c,length=p)
function points(N::A,c=1,x0=0) where A<:Profile{p} where p
    Chain{SubManifold(ℝ^3),1}.(1.0,interval(N,c,x0),profile(N))
end

const AppendixI = [0,.005,.0125,.025,.05,.075,.1,.15,.2,.25,.3,.4,.5,.6,.7,.8,.9,.95,1]
const AppendixII = [0,.005,0.0075,.0125,.025,.05,.075,.1,.15,.2,.25,.3,.35,.4,.45,.5,.55,.6,.65,.7,.75,.8,.85,.9,.95,1]

Base.adjoint(n::P) where P<:Profile = Adjoint(n)
(n::LinearAlgebra.Adjoint{Any,<:Profile})(x) = profileslope(n.parent,x)
(n::LinearAlgebra.Adjoint{Any,<:Profile})(x,c,x0=0) = profileslope(n.parent,x,c,x0)

# flat plate

struct FlatPlate{p} <: Profile{p}
    @pure FlatPlate{p}() where p = new{p}()
end

@pure (::FlatPlate)(x) = 0.0
@pure (n::FlatPlate)(x,c,x0=0)  = profile(n,x,c,x0)
@pure profile(::FlatPlate,x,c,x0=nothing) = 0.0
@pure profileslope(::FlatPlate,x,c=nothing,x0=nothing) = 0.0
@pure profile(f::FlatPlate{p}) where p = range(0,0,length=p)
@pure profileslope(f::FlatPlate) = profile(f)

# parabolic arc

struct ParabolicArc{t,p} <: Profile{p}
    @pure ParabolicArc{t,p}() where {t,p} = new{t,p}()
    @pure ParabolicArc{p}() where p = new{6,p}()
end

(n::ParabolicArc)(x,c,x0=0)  = profile(n,x,c,x0)
@pure (::ParabolicArc{t})(x) where t = (x<0 || x>1) ? 0.0 : (t/25)*x*(1-x)
@pure profileslope(::ParabolicArc{t},x) where t = (x<0 || x>1) ? 0.0 : (1-2x)*t/400

# circular arc

struct CircularArc{t,p} <: Profile{p}
    @pure CircularArc{t,p}() where {t,p} = new{t,p}()
    @pure CircularArc{p}() where p = new{6,p}()
end

(n::CircularArc)(x,c,x0=0)  = profile(n,x,c,x0)
@pure function (::CircularArc{T})(x) where T
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

# NACA Clark Y
# http://naca.central.cranfield.ac.uk/reports/1933/naca-report-431.pdf

struct ClarkY{t,te,p} <: Profile{p}
    @pure ClarkY{t,te,p}() where {t,te,p} = new{t,te,p}()
    @pure ClarkY{t,p}() where {t,p} = ClarkY{t,0.21,p}()
end

@pure Y(x) = Chain{ℝ5,1}(sqrt(x),x,x^2,x^3,x^4)
@pure dY(x) = Chain{ℝ5,1}(0.5/sqrt(x),1.0,2x,3x^2,4x^3)

@pure function thickness(x,a,t)
    (x<0 || x>1) ? 0.0 : 5t*(Y(x)⋅a)[1]
end # thickness distribution
@pure function thickness_slope(x,a,t)
    (x<0 || x>1) ? 0.0 : 5t*(dY(x)⋅a)[1]
end # thickness distribution slope

@pure clarky(te=0.0021) = Chain{ℝ5,1}(0.2969,-0.1260,-0.3516+1e-16,0.2843,te-0.1036)
@pure radius(t,::Val{a0}=Val(clarky()[1])) where a0 = t^2*((a0/0.2)^2/2)

(n::ClarkY)(x,c,x0=0)  = profile(n,x,c,x0)
@pure (::ClarkY{t,te})(x) where {t,te} = thickness(x,clarky(te/100),t/100)
@pure profileslope(::ClarkY{t,te},x) where {t,te} = thickness_slope(x,clarky(te/100),t/100)

# NACA 4-series thickness
# https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19930091108.pdf

struct Thickness{t,x,te,p} <: Profile{p}
    @pure Thickness{t,x,te,p}() where {t,x,te,p} = new{t,x,te,p}()
    @pure Thickness{t,x,p}() where {t,x,p} = Thickness{t,x,t/100,p}()
    @pure Thickness{t,p}() where {t,p} = Thickness{t,3,p}()
end

#[0.2,0.3,0.4,0.5,0.6] -> -[0.2,0.234,0.315,0.465,0.7]
@pure tailslope(x) = -approx(x,SVector(31/200,151/300,-109/40,43/6,-5/2))
@pure riegel(x) = -(2.24-5.42x+12.3x^2)/(10(1-0.878x))

@pure function clarky(te,x,ts=tailslope(x),n=0.1,y=0.078) # te,x = 0.002,0.3
    transpose(Chain{ℝ5,1}(Y(x),dY(x),Y(1.0),dY(1.0),Y(n)))\Chain{ℝ5,1}(0.1,0.0,te,ts,y)
end

(n::Thickness)(x,c,x0=0)  = profile(n,x,c,x0)
@pure @generated (::Thickness{t,m,te})(x) where {t,m,te} = :(thickness(x,$(clarky(te/100,m/10)),$(t/100)))
@pure @generated profileslope(::Thickness{t,m,te},x) where {t,m,te} = :(thickness_slope(x,$(clarky(te/100,m/10)),$(t/100)))

# NACA mofified thickness
# http://naca.central.cranfield.ac.uk/reports/1935/naca-report-492.pdf

struct Modified{t,m,te,p} <: Profile{p}
    @pure Modified{t,m,te,p}() where {t,m,te,p} = new{t,m,te,p}()
    @pure Modified{t,m,p}() where {t,m,p} = new{t,m,0.2,p}()
    @pure Modified{t,p}() where {t,p} = Modified{t,63,p}()
end

@pure D(x) = (x1 = 1-x; Chain{ℝ4,1}(1.0,x1,x1^2,x1^3))
@pure dD(x) =  (x1 = x-1; Chain{ℝ4,1}(0.0,-1.0,2x1,-3x1^2))
@pure ddD(x) =  Chain{ℝ4,1}(0.0,0.0,2.0,-6(x-1))

@pure A(x) = Chain{ℝ4,1}(sqrt(x),x,x^2,x^3)
@pure dA(x) =  Chain{ℝ4,1}(0.5/sqrt(x),1.0,2x,3x^2)
@pure ddA(x) =  Chain{ℝ4,1}(-1/(4x*sqrt(x)),0.0,2.0,6x)

@pure function modified(t,i,x,te=0.002,ts=tailslope(x))
    d = transpose(Chain{ℝ4,1}(D(x),dD(x),D(1.0),dD(1.0)))\Chain{ℝ4,1}(0.1,0.0,te,ts)
    a = transpose(Chain{ℝ4,1}(A(x),dA(x),ddA(x),Chain{ℝ4,1}(1.0,1.0,1.0,1.0)))
    return a\Chain{ℝ4,1}(0.1,0.0,(ddD(x)⋅d)[1],(1/5t)*sqrt(2radius(0.2*i/6))),d
end

@pure function modified(::Modified{t,m,te}) where {t,m,te}
    i = Int(m÷10); t/100,i≠9 ? i : 6sqrt(3),(m%10)/10,te/100
end

(n::Modified)(x,c,x0=0)  = profile(n,x,c,x0)
@pure function (P::Modified)(x)
    t,i,p,te = modified(P); a = modified(t,i,p,te)
    (x<0 || x>1) ? 0.0 : 5t*(x < p ? A(x)⋅a[1] : D(x)⋅a[2])[1]
end
@pure function profileslope(P::Modified,x)
    t,i,p,te = modified(P); a = modified(t,i,p,te)
    (x<0 || x>1) ? 0.0 : 5t*(x < p ? dA(x)⋅a[1] : dD(x)⋅a[2])[1]
end

# NACA camber
# https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19930091610.pdf

struct NACA4{n,p} <: Profile{p}
    @pure NACA4{n,p}() where {n,p} = new{n,p}()
end

@pure C(x) = Chain{ℝ3,1}(1.0,x,x^2)
@pure dC(x) = Chain{ℝ3,1}(0.0,1.0,2x)

@pure function naca4(m,p)
    cm,y,dy = Chain{ℝ3,1}(m,0.0,0.0),C(p),dC(p)
    transpose(Chain{ℝ3,1}(y,dy,C(0.0)))\cm,
    transpose(Chain{ℝ3,1}(y,dy,C(1.0)))\cm
end

@pure function naca4(n) # decode
    ns = string(n,pad=2)
    M,P = parse.(Int,(ns[1],ns[2]))
    return M/100,P/10
end

(n::NACA4)(x,c,x0=0)  = profile(n,x,c,x0)
@pure function (::NACA4{n})(x) where n
    (x<0 || x>1) && (return 0.0)
    m,p = naca4(n)
    (C(x)⋅naca4(m,p)[x < p ? 1 : 2])[1]
end
@pure function profileslope(::NACA4{n},x) where n
    (x<0 || x>1) && (return 0.0)
    m,p = naca4(n)
    (dC(x)⋅naca4(m,p)[x < p ? 1 : 2])[1]
end

# NACA camber 5 series
# https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19930091610.pdf

struct NACA5{n,p} <: Profile{p}
    @pure NACA5{n,p}() where {n,p} = new{n,p}()
end

@pure function naca5(n) # decode
    ns=string(n)
    C,P,R = parse.(Int,(ns[1],ns[2],ns[3]))
    m,p = C/2,P/20 # Cl = 3C/20
    r,k1,k1k2 = if !Bool(R)
        #p: [0.05,0.1,0.15,0.2,0.25]
        #m: p -> [0.058,0.126,0.2025,0.29,0.391]
        approx(p,SVector(-1/250,359/300,7/10,10/3,0)),
        #k: p -> [361.4,51.64,15.957,6.643,3.230]
        approx(p,SVector(284037/200,-3296847/100,4296829/15,-1087744,4544800/3)),0
    else # reflexed
        #p: [0.1,0.15,0.2,0.25]
        #m: p -> [0.13,0.217,0.318,0.441]
        approx(p,SVector(17/500,26/15,-2,32/3)),
        #k: p -> [51.99,15.793,6.52,3.191]
        approx(p,SVector(72269/250,-583261/150,89864/5,83920/3)),
        #k/k: p-> [0.000764,0.00677,0.0303,0.1355]
        approx(p,SVector(10763/50000,12081/25000,-87457/2500,10691/125))
    end
    return m*k1,r,k1k2
end

(n::NACA5)(x,c,x0=0)  = profile(n,x,c,x0)
@pure function (::NACA5{n})(x) where n
    (x<0 || x>1) && (return 0.0)
    mk1,p,k1k2 = naca5(n)
    mk1*(p^3*(1-x)+(if x < p
        (x-p)^3-(k1k2≠0 ? k1k2*x*(1-p)^3 : 0)
    else
        k1k2≠0 ? k1k2*((x-p)^3-x*(1-p)^3) : 0
    end))/6
end
@pure function profileslope(::NACA5{n},x) where n
    (x<0 || x>1) && (return 0.0)
    mk1,p,k1k2 = naca5(n)
    mk1*(-p^3+(if x < p
        3(p-x)^2-(k1k2≠0 ? k1k2*(1-p)^3 : 0)
    else
        k1k2≠0 ? 3k1k2*((p-x)^2-(1-p)^3) : 0
    end))
end

# NACA camber 6 series
# https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19930091610.pdf

struct NACA6{c,n,p} <: Profile{p}
    a::SVector{n,Float64}
    Cl::SVector{n,Float64}
    @pure NACA6{p}(a,Cl) where p = new{10sum(Cl),length(Cl),p}(a,Cl)
    @pure NACA6{c,p}() where {c,p} = NACA6{p}(SVector(1.0),SVector(c/10))
end # angle of attack: αᵢ = Cla*h

@pure function naca6(x,a,g,h)
    x1 = 1 - x
    if a == 1.0
        -(x1*log(x1)+x*log(x))
    else
        a1,ax = 1-a,a-x
        x12,ax2 = x1^2,ax^2
        ((ax2*log(abs(ax))-x12*log(x1))+(x12-ax2)/2)/2a1-x*log(x)+g-h*x
    end
end
@pure function naca6(x,a,h)
    if a == 1.0
        log(1-x)-log(x)
    else
        x1,ax = 1-x,a-x
        (x1*log(x1)-ax*log(ax))/(1-a)-log(x)-1-h
    end
end

@pure function naca6(n::NACA6) # decode
    a1 = 1.0.-n.a; a12 = (a1.^2)./2
    g = .-((n.a.^2).*(log.(n.a)./2.0.-1/4).+1/4)./a1
    h = (a12.*log.(a1).-a12./2)./a1 .+ g
    #h = a1.*(log.(a1)./2.0.-1/4) .+ g
    return n.Cl./(2π.*(1.0.+n.a)),n.a,h,g
end

(n::NACA6)(x,c,x0=0)  = profile(n,x,c,x0)
@pure function (n::NACA6)(x)
    (x<1e-7 || x>1-1e-7) && (return 0.0)
    Cla,a,g,h = naca6(n)
    sum(Cla.*naca6.(x,a,g,h))
end
@pure function profileslope(n::NACA6,x)
    (x<1e-7 || x>1-1e-7) && (return 0.0)
    Cla,a,h = naca6(n)
    sum(Cla.*naca6.(x,a,h))
end

# NACA camber 6A series
# https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19930091610.pdf

struct NACA6A{c,p} <: Profile{p}
    @pure NACA6A{c,p}() where {c,p} = new{c,p}()
end

(n::NACA6A)(x,c,x0=0)  = profile(n,x,c,x0)
@pure function (::NACA6A{c})(x) where c
    (x<0 || x>1) && (return 0.0)
    (c/36π)*(if x < 0.87437
        sum.(naca6(x,0.8,SVector(0),SVector(0)))
    else # (0.0302164,-0.245209)/2π(1+a)
        0.34173943292855025-2.773248454778759(x-0.87437)
    end)
end
@pure function profileslope(::NACA6A{c},x) where c
    (x<0 || x>1) && (return 0.0)
    x < 0.87437 ? (c/36π)*sum.(naca6(x,0.8,SVector(0))) : -0.0245209c
end
