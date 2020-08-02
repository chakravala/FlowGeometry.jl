
#   This file is part of FlowGeometry.jl. It is licensed under the AGPL license
#   FlowGeometry Copyright (C) 2020 Michael Reed

approx(x,y) = [x^i for i ∈ 0:length(y)-1]⋅y

# profiles

abstract type Profile{P} end

interval(::Profile{p}) where p = interval(p)
interval(p::Int,c=1,x0=0) = range(x0,c,length=p)

# flat plate

struct FlatPlate{p} <: Profile{p}
    @pure FlatPlate{p}() where p = new{p}()
end

profile(::FlatPlate{p}) where p = range(0,0,length=p)
profileslope(f::FlatPlate) = profile(f)

# parabolic arc

struct ParabolicArc{n,p} <: Profile{p}
    @pure ParabolicArc{n,p}() where {n,p} = new{n,p}()
end

function parabolic(x,t=0.06,c=1,x0=0)
    x0c = x0+c
    (x<x0 || x>x0c) && return 0.0
    4t*(x-x0)*(x0c-x)
end

parabolicslope(x::Chain,t=0.06,c=1,x0=0) = parabolicslope(x[2],t,c,x0)
function parabolicslope(x,t=0.06,c=1,x0=0)
    (x<x0 || x>x0+c) && return 0.0
    t*(8(x0-x)+4c)
end

profile(::ParabolicArc{t,p}) where {t,p} = parabolic.(interval(p),t/100)
profileslope(::ParabolicArc{t,p}) where {t,p} = parabolicslope.(interval(p),t/100)

# circular arc

struct CircularArc{n,p} <: Profile{p}
    @pure CircularArc{n,p}() where {n,p} = new{n,p}()
    @pure CircularArc{p}() where p = new{6,p}()
end

function arc(x,t=0.06,c=1,x0=0)
    (x<x0 || x>x0+c) && return 0.0
    r = ((c/2)^2+t^2)/2t
    r*sin(acos((x-x0-c/2)/r))-r+t
    #ct = c^2-4t^2
    #((sqrt(-16(2(x-x0)-c)^2*t^2 + ct^2) - ct)*ct) / (8abs(ct)*t)
end

arcslope(x::Chain,t=0.06,c=1,x0=0) = arcslope(x[2],t,c,x0)
function arcslope(x,t=0.06,c=1,x0=0)
    (x<x0 || x>x0+c) && return 0.0
    #r = ((c/2)^2+t^2)/2t # Algebra.df(arc(:x),:x)
    #((c-2(x-x0))*r) / (sqrt(4r^2-((2(x-x0) - c)^2))*abs(r))
    xc,t2 = 2(x-x0)-c,t^2
    (-4xc*t2) / (sqrt((c^2 + 4t2)^2 - 16xc^2*t2)*t)
end

profile(::CircularArc{t,p}) where {t,p} = arc.(interval(p),t/100)
profileslope(::CircularArc{t,p}) where {t,p} = arcslope.(interval(p),t/100)

# NACA 4-series thickness

struct ClarkY{n,p} <: Profile{p}
    @pure ClarkY{n,p}() where {n,p} = new{n,p}()
end

Y(x) = SVector(sqrt(x),x,x^2,x^3,x^4)
dY(x) = SVector(0.5/sqrt(x),1.0,2x,3x^2,4x^3)

const _clarky = SVector(0.2969,-0.1260,-0.3515,0.2843,-0.1036)

thickness(x::T,a,t=1) where T<:AbstractRange = thickness.(x,Ref(a),t,x[end],x[1])
thickness(x::T,a,t=1) where T<:AbstractVector = thickness.(x,Ref(a),t,x[end],x[1])
function thickness(x,a,t=1,c=1,x0=0)
    (x<x0 || x>x0+c) ? (return 0.0) : 5t*(Y((x-x0)/c)⋅a)
end # thickness distribution

thickness_slope(x::T,a,t=1) where T<:AbstractRange = thickness_slope.(x,Ref(a),t,x[end],x[1])
thickness_slope(x::T,a,t=1) where T<:AbstractVector = thickness_slope.(x,Ref(a),t,x[end],x[1])
function thickness_slope(x,a,t=1,c=1,x0=0)
    (x<x0 || x>x0+c) ? (return 0.0) : 5t*(dY((x-x0)/c)⋅a)
end # thickness distribution slope

profile(::ClarkY{t,p}) where {t,p} = thickness(interval(p),_clarky,t/100)
profileslope(::ClarkY{t,p}) where {t,p} = thickness_slope(interval(p),_clarky,t/100)

# NACA custom thickness

struct Thickness{t,x,te,p} <: Profile{p}
    @pure Thickness{t,x,te,p}() where {t,x,te,p} = new{t,x,te,p}()
    @pure Thickness{t,x,p}() where {t,x,p} = Thickness{t,x,0.2,p}()
    @pure Thickness{t,p}() where {t,p} = Thickness{t,3,p}()
end

#[0.2,0.3,0.4,0.5,0.6] -> -[0.2,0.234,0.315,0.465,0.7]
tailslope(x) = -approx(x,[31/200,151/300,-109/40,43/6,-5/2])
riegel(x) = -(2.24-5.42x+12.3x^2)/(10(1-0.878x))

function clarky(x=0.3,te=0.002,ts=tailslope(x),n=0.1,y=0.078)
    SMatrix{5,5}(Y(x)...,dY(x)...,Y(1)...,dY(1)...,Y(n)...)'\SVector(0.1,0,te,ts,y)
end

profile(::Thickness{t,x,te,p}) where {t,x,te,p} = thickness(interval(p),clarky(x/10,te/100),t/100)
profileslope(::Thickness{t,x,te,p}) where {t,x,te,p} = thickness_slope(interval(p),clarky(x/10,te/100),t/100)

# NACA mofified thickness

struct Modified{t,m,p} <: Profile{p}
    @pure Modified{t,m,p}() where {t,m,p} = new{t,m,p}()
    @pure Modified{t,p}() where {t,m} = Modified{t,63,p}()
end

D(x) = (x1 = 1-x; SVector(1.0,x1,x1^2,x1^3))
dD(x) =  (x1 = x-1; SVector(0.0,-1.0,2x1,-3x1^2))
ddD(x) =  SVector(0.0,0.0,2,-6(x-1))

A(x) = SVector(sqrt(x),x,x^2,x^3)
dA(x) =  SVector(0.5/sqrt(x),1.0,2x,3x^2)
ddA(x) =  SVector(-1/(4x*sqrt(x)),0.0,2,6x)

modified(x::T,a,p,t=1) where T<:AbstractRange = modified.(x,Ref(a),p,t,x[end],x[1])
modified(x::T,a,p,t=1) where T<:AbstractVector = modified.(x,Ref(a),p,t,x[end],x[1])
function modified(x,a,p,t=1,c=1,x0=0)
    (x<x0 || x>x0+c) && (return 0.0)
    if x < p
        5t*(A((x-x0)/c)⋅a[1])
    else
        5t*(D((x-x0)/c)⋅a[2])
    end
end # modified thickness

modified_slope(x::T,a,p,t=1) where T<:AbstractRange = modified_slope.(x,Ref(a),p,t,x[end],x[1])
modified_slope(x::T,a,p,t=1) where T<:AbstractVector = modified_slope.(x,Ref(a),p,t,x[end],x[1])
function modified_slope(x,a,p,t=1,c=1,x0=0)
    (x<x0 || x>x0+c) && (return 0.0)
    if x < p
        5t*(dA((x-x0)/c)⋅a[1])
    else
        5t*(dD((x-x0)/c)⋅a[2])
    end
end # thickness distribution slope

function profile(::Modified{t,m,p}) where {t,m,p}
    i,x = Int(m÷10),(m%10)/10
    modified(interval(p),coef(t/100,i≠9 ? i : 6sqrt(3),x),x,t/100)
end
function profileslope(::Modified{t,m,p}) where {t,m,p}
    i,x = Int(m÷10),(m%10)/10
    modified_slope(interval(p),coef(t/100,i≠9 ? i : 6sqrt(3),x),x,t/100)
end

radius(t,::Val{a0}=Val(_clarky[1])) where a0 = t^2*((a0/0.2)^2/2)

#transpose(Chain{V,1}(Chain{V,1}(D(x)),Chain{V,1}(dD(x)),Chain{V,1}(D(1)),Chain{V,1}(dD(1))))\Chain{V,1}(0.1,0.0,0.002,-tailangle(x))
#=function coef(t,i,x,te=0.002,ts=tailslope(x))
    di = SMatrix{4,4}(D(x)...,dD(x)...,D(1)...,dD(1)...)'\SVector(0.1,0,te,ts)
    ai = SMatrix{4,4}(A(x)...,dA(x)...,ddA(x)...,1,0,0,0)'\SVector(0.1,0,ddD(x)⋅di,(1/5t)*sqrt(2radius(0.1*i/6)))
    ai,di
end=#
function coef(t,i,x,te=0.002,ts=tailslope(x))
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

function camber2(n) # decode
    ns = string(n,pad=2)
    M,P = parse.(Int,(ns[1],ns[2]))
    return M/100,P/10
end
function camber3(n) # decode
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

function clarkY(m,p)
    SMatrix{3,3}(yc(0)...,yc(p)...,dyc(p)...)'\SVector(0,m,0),
    SMatrix{3,3}(yc(1)...,yc(p)...,dyc(p)...)'\SVector(0,m,0)
end=#

function camber2(x,p,c=1,x0=0)
    (x<x0 || x>x0+c) && return 0.0
    xc = (x-x0)/c
    return if xc < p
        xcp = xc/p
        2xcp-xcp^2
    else
        (1+2p*(xc-1)-xc^2)/(1-p)^2
    end
end
function camber2slope(x,p,c=1,x0=0)
    (x<x0 || x>x0+c) && return 0.0
    xc = (x-x0)/c
    2(p-xc)/(xc<p ? p^2 : (1-p)^2)
end
function camber3(x,p,k,c=1,x0=0)
    (x<x0 || x>x0+c) && 0.0
    xc = (x-x0)/c
    (p^3*(1-xc)+(if xc < p
        #xc^3-3p*xc^2+p^2*(3-p)*xc
        (xc-p)^3-(k≠0 ? k*xc*(1-p)^3 : 0)
    else
        k≠0 ? k*((xc-p)^3-xc*(1-p)^3) : 0
    end))/6
end
function camber3slope(x,p,k,c=1,x0=0)
    (x<x0 || x>x0+c) && 0.0
    xc = (x-x0)/c
    -p^3+(if xc < p
        #3xc^2-6p*xc+p^2*(3-p)
        3(p-xc)^2-(k≠0 ? k*(1-p)^3 : 0)
    else
        k≠0 ? 3k*((p-xc)^2-(1-p)^3) : 0
    end)
end

function profile(::Camber{n,p}) where {n,p}
    x = interval(p)
    if n < 100
        m,r = camber2(n)
        return m.*camber2.(x,r,x[end],x[1])
    else
        mk1,r,k1k2 = camber3(n)
        return mk1*camber3.(x,r,k1k2,x[end],x[1])
    end
end
function profileslope(::Camber{n,p}) where {n,p}
    x = interval(p)
    if n < 100
        m,r = camber2(n)
        return m*camber2slope.(x,r,x[end],x[1])
    else
        mk1,r,k1k2 = camber3(n)
        return mk1.*camber3slope.(x,r,k1k2,x[end],x[1])
    end
end
