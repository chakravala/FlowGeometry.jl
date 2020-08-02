
#   This file is part of FlowGeometry.jl. It is licensed under the AGPL license
#   FlowGeometry Copyright (C) 2020 Michael Reed

abstract type Airfoil{p} end

# upper

struct UpperArc{A,p} <: Profile{p}
    a::A
    @pure UpperArc{A,p}(a) where {A,p} = new{A,p}(a)
end

@pure UpperArc(s::S) where S<:Airfoil{P} where P = UpperArc{S,length(upper(s))}(s)

interval(U::UpperArc)  = real.(upper(U.a))
profile(U::UpperArc) = imag.(upper(U.a))
profileslope(f::UpperArc,e=chordedges(f)) = interp(e,column(gradient(e,profile(f))))

# lower

struct LowerArc{A,p} <: Profile{p}
    a::A
    @pure LowerArc{A,p}(a) where {A,p} = new{A,p}(a)
end

@pure LowerArc(s::S) where S<:Airfoil{P} where P = LowerArc{S,length(lower(s))}(s)

interval(U::LowerArc) = real.(lower(U.a))
profile(U::LowerArc) = imag.(lower(U.a))
profileslope(f::LowerArc,e=chordedges(f)) = interp(e,column(gradient(e,profile(f))))

# Airfoils upper lower

upper(z,r) = (U = z.+r; U[end] = 1; U) # upper surface points xu+im*yu
lower(z,r) = (L = z.-r; L[end] = 1; L) # lower surface points xl+im*yl
upper(z::A) where A<:Profile = upper(interval(z),profile(z).*im)
lower(z::A) where A<:Profile = lower(interval(z),profile(z).*im)
upperlower(z::A) where A<:Profile = upperlower(interval(z),profile(z).*im)
upperlower(z::A) where A<:Airfoil = upper(z),lower(z)
upperlower(z,r) = upper(z,r),lower(z,r)

# Joukowski

struct Joukowski{R,f,g,b,p} <: Airfoil{p}
    @pure Joukowski{R,f,g,b,p}() where {R,f,g,b,p} = new{R,f,g,b,p}()
end

joukowski(R,f,g,b,p) = Joukowski{R,f,g,b,p}()

# SymmetricArc

struct SymmetricArc{S,P} <: Airfoil{P}
    s::S
    @pure SymmetricArc{S,P}(s) where {S,P} = new{S,P}(s)
end

@pure SymmetricArc(s::S) where S<:Profile{P} where P = SymmetricArc{S,2P-2}(s)

upper(n::SymmetricArc) = upper(n.s)
lower(n::SymmetricArc) = lower(n.s)
upperlower(n::SymmetricArc)= upperlower(n.s)

# DoubleArc

struct DoubleArc{U,L,P} <: Airfoil{P}
    u::U
    l::L
    @pure DoubleArc{U,L,P}(u,l) where {U,L,P} = new{U,L,P}(u,l)
end

@pure DoubleArc(u::U,l::L) where {U<:Profile{P},L<:Profile{Q}} where {P,Q} = DoubleArc{U,L,P+Q-2}(u,l)

upper(n::DoubleArc{U,L,P}) where {U,L,P} = upper(n.u)
lower(n::DoubleArc{U,L,P}) where {U,L,P} = lower(n.l)

# NACA

struct NACA{C,T,P} <: Airfoil{P} # NACA custom
    c::C
    t::T
    @pure NACA{C,T,P}(c,t) where {C,T,P} = new{C,T,P}(c,t)
end

@pure NACA{n,p}() where {n,p} = NACA(Camber{Int(n÷100),p}(),ClarkY{n%100,p}())
@pure NACA(c::C,t::T) where {C<:Profile{P},T<:Profile{P}} where P = NACA{C,T,2P-2}(c,t)

upper(n::NACA) = upper(getprofiles(n)...)
lower(n::NACA) = lower(getprofiles(n)...)
upperlower(n::NACA) = upperlower(getprofiles(n)...)
upper(x,yc,dyc_dx,yt) = upper(x.+im.*yc,im.*cis.(atan.(dyc_dx)).*yt)
lower(x,yc,dyc_dx,yt) = lower(x.+im.*yc,im.*cis.(atan.(dyc_dx)).*yt)
upperlower(x,yc,dyc_dx,yt) = upperlower(x.+im.*yc,im.*cis.(atan.(dyc_dx)).*yt)
getprofiles(n::NACA) = interval(n.c),profile(n.c),profileslope(n.c),profile(n.t)

# upper lower slope

function getprofileslopes(n::NACA{c,t,p}) where {c,t,p}
    range(0,1,length=p),profile(n.c),profileslope(n,c),profileslope(n.t)
end

#=upperlowerslope(n::NACA) = upperlowerslope(getprofiles(n)...)
function upperlowerslope(x,yc,dyc_dx,dyt_dx)
    ((1+dyc_dx.^2).*dyt_dx .- yt.*dyc_dx.*dyt_dxdx)/(1 + dyc_dx.^2)^(3/2)
    #upperlower(yc,cos.(θ).*dyt_dx)
end=#

# complex point sets

function Base.complex(::Joukowski{R,f,g,b,p}) where {R,f,g,b,p}
    θ = range(0,2π,length=2p-1)
    z = R*cis.(θ).-(f-g*im)
    z.+b^2*inv.(z)
end
function Base.complex(N::A) where A<:Airfoil
    U,L = upperlower(N)
    [U;reverse(L)[2:end]]
end
function points(N::A) where A<:Airfoil
    z = complex(N)[1:end-1]
    Chain{SubManifold(ℝ^3),1}.(1.0,real.(z),imag.(z))
end
function points(N::A) where A<:Profile{p} where p
    Chain{SubManifold(ℝ^3),1}.(1.0,interval(N),profile(N))
end

