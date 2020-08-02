
#   This file is part of FlowGeometry.jl. It is licensed under the AGPL license
#   FlowGeometry Copyright (C) 2020 Michael Reed

abstract type Airfoil{p} end

export Airfoil, UpperArc, LowerArc, SymmetricArc, DoubleArc, NACA, Joukowski
export upper, lower, upperlower

# Airfoils upper lower

for Side ∈ ("Upper","Lower")
    Arc,side = Symbol(Side,:Arc),Symbol(lowercase(Side))
    @eval begin
        struct $Arc{A,p} <: Profile{p}
            a::A
            @pure $Arc{A,p}(a) where {A,p} = new{A,p}(a)
        end
        @pure $Arc(s::S) where S<:Airfoil{P} where P = $Arc{S,length($side(s))}(s)
        $side(z::A,c=1,x0=0) where A<:Profile = $side(interval(z,c,x0),profile(z).*(c*im))
        interval(U::$Arc,c=1,x0=0)  = real.($side(U.a,c,x0))
        profile(U::$Arc) = imag.($side(U.a))
        profileslope(f::$Arc,e=initedges(f)) = interp(e,column(gradient(e,profile(f))))
    end
end

upper(z,r) = (U = z.+r; U[end] = 1; U) # upper surface points xu+im*yu
lower(z,r) = (L = z.-r; L[end] = 1; L) # lower surface points xl+im*yl
upperlower(z::A,c=1,x=0) where A<:Profile = upperlower(interval(z,c,x),profile(z).*(c*im))
upperlower(z::A,c=1,x=0) where A<:Airfoil = upper(z,c,x),lower(z,c,x)
upperlower(z,r) = upper(z,r),lower(z,r)

function Base.complex(N::A) where A<:Airfoil
    U,L = upperlower(N)
    [U;reverse(L)[2:end]]
end
function points(N::A) where A<:Airfoil
    z = complex(N)[1:end-1]
    Chain{SubManifold(ℝ^3),1}.(1.0,real.(z),imag.(z))
end

# SymmetricArc

struct SymmetricArc{S,P} <: Airfoil{P}
    s::S
    @pure SymmetricArc{S,P}(s) where {S,P} = new{S,P}(s)
end

@pure SymmetricArc(s::S) where S<:Profile{P} where P = SymmetricArc{S,2P-2}(s)

upper(n::SymmetricArc,c=1,x0=0) = upper(n.s,c,x0)
lower(n::SymmetricArc,c=1,x0=0) = lower(n.s,c,x0)
upperlower(n::SymmetricArc,c=1,x0=0)= upperlower(n.s,c,x0)

# DoubleArc

struct DoubleArc{U,L,P} <: Airfoil{P}
    u::U
    l::L
    @pure DoubleArc{U,L,P}(u,l) where {U,L,P} = new{U,L,P}(u,l)
end

@pure DoubleArc(u::U,l::L) where {U<:Profile{P},L<:Profile{Q}} where {P,Q} = DoubleArc{U,L,P+Q-2}(u,l)

upper(n::DoubleArc{U,L,P},c=1,x0=0) where {U,L,P} = upper(n.u,c,x0)
lower(n::DoubleArc{U,L,P},c=1,x0=0) where {U,L,P} = lower(n.l,c,x0)

# NACA

struct NACA{C,T,P} <: Airfoil{P} # NACA custom
    c::C
    t::T
    @pure NACA{C,T,P}(c,t) where {C,T,P} = new{C,T,P}(c,t)
end

@pure NACA{n,p}() where {n,p} = NACA(Camber{Int(n÷100),p}(),ClarkY{n%100,p}())
@pure NACA(c::C,t::T) where {C<:Profile{P},T<:Profile{P}} where P = NACA{C,T,2P-2}(c,t)

upper(n::NACA,c=1,x0=0) = upper(getprofiles(n,c,x0)...)
lower(n::NACA,c=1,x0=0) = lower(getprofiles(n,c,x0)...)
upperlower(n::NACA,c=1,x0=0) = upperlower(getprofiles(n,c,x0)...)
upper(x,yc,dyc_dx,yt) = upper(x.+im.*yc,im.*cis.(atan.(dyc_dx)).*yt)
lower(x,yc,dyc_dx,yt) = lower(x.+im.*yc,im.*cis.(atan.(dyc_dx)).*yt)
upperlower(x,yc,dyc_dx,yt) = upperlower(x.+im.*yc,im.*cis.(atan.(dyc_dx)).*yt)
getprofiles(n::NACA,c=1,x0=0) = interval(n.c,c,x0),c*profile(n.c),profileslope(n.c),c*profile(n.t)

# Joukowski

struct Joukowski{R,f,g,b,p} <: Airfoil{p}
    @pure Joukowski{R,f,g,b,p}() where {R,f,g,b,p} = new{R,f,g,b,p}()
end

joukowski(R,f,g,b,p) = Joukowski{R,f,g,b,p}()

function Base.complex(::Joukowski{R,f,g,b,p}) where {R,f,g,b,p}
    θ = range(0,2π,length=2p-1)
    z = R*cis.(θ).-(f-g*im)
    z.+b^2*inv.(z)
end
