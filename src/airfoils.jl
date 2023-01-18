
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

abstract type Airfoil{p} end

export Airfoil, UpperArc, LowerArc, SymmetricArc, DoubleArc, @NACA_str, Joukowski
export upper, lower, upperlower, American, British

struct Chord{p,c,x0}
    @pure Chord{p,c,x0}() where {p,c,x0} = new{p,c,x0}()
    @pure Chord{p,c}() where {p,c} = new{p,c,0}()
    @pure Chord{p}() where p = Chord{p,1}()
    @pure Chord(p=150,c=1,x0=0) = new{p,c,x0}()
end

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

@pure chord(p) = range(0,0,length=p)
@pure chord(::Profile{p}) where p = chord(p)
@pure chord(::Airfoil{p}) where p = chord(p)
upper(z,r) = (U = z.+r; U[end] = 1; U) # upper surface points xu+im*yu
lower(z,r) = (L = z.-r; L[end] = 1; L) # lower surface points xl+im*yl
upper(x,yc,dyc_dx,yt) = upper(x.+im.*yc,im.*cis.(atan.(dyc_dx)).*yt)
lower(x,yc,dyc_dx,yt) = lower(x.+im.*yc,im.*cis.(atan.(dyc_dx)).*yt)
upperlower(x,yc,dyc_dx,yt) = upperlower(x.+im.*yc,im.*cis.(atan.(dyc_dx)).*yt)
upperlower(z::A,c=1,x=0) where A<:Profile = upperlower(interval(z,c,x),profile(z).*(c*im))
upperlower(z::A,c=1,x=0) where A<:Airfoil = upper(z,c,x),lower(z,c,x)
upperlower(z,r) = upper(z,r),lower(z,r)

function Base.complex(N::A) where A<:Airfoil
    U,L = upperlower(N)
    [U;reverse(L)[2:end]]
end
function points(N::A) where A<:Airfoil
    z = complex(N)[1:end-1]
    Chain{Submanifold(ℝ^3),1}.(1.0,real.(z),imag.(z))
end

# SymmetricArc

"""
    SymmetricArc{U,L,P} <: Airfoil

Symmetric arc airfoil constructed from a profile without camber applied.
"""
struct SymmetricArc{S,P} <: Airfoil{P}
    s::S
    @pure SymmetricArc{S,P}(s) where {S,P} = new{S,P}(s)
end

@pure SymmetricArc(s::S) where S<:Profile{P} where P = SymmetricArc{S,2P-2}(s)

upper(n::SymmetricArc,c=1,x0=0) = upper(n.s,c,x0)
lower(n::SymmetricArc,c=1,x0=0) = lower(n.s,c,x0)
upperlower(n::SymmetricArc,c=1,x0=0)= upperlower(n.s,c,x0)

# DoubleArc

"""
    DoubleArc{U,L,P} <: Airfoil

Asymmetric double arc airfoil constructed from upper and lower profile.
"""
struct DoubleArc{U,L,P} <: Airfoil{P}
    u::U
    l::L
    @pure DoubleArc{U,L,P}(u,l) where {U,L,P} = new{U,L,P}(u,l)
end

@pure DoubleArc(u::U,l::L) where {U<:Profile{P},L<:Profile{Q}} where {P,Q} = DoubleArc{U,L,P+Q-2}(u,l)

upper(n::DoubleArc{U,L,P},c=1,x0=0) where {U,L,P} = upper(n.u,c,x0)
lower(n::DoubleArc{U,L,P},c=1,x0=0) where {U,L,P} = lower(n.l,c,x0)

# Rotated thickness by camber slope

"""
    British{C,T,P} <: Airfoil

Thickness measured perpendicular to the chord line, British convention.
"""
struct British{C,T,P} <: Airfoil{P}
    c::C
    t::T
    @pure British{C,T,P}(c,t) where {C,T,P} = new{C,T,P}(c,t)
end

@pure British{n,p}() where {n,p} = British(Camber{Int(n÷100),p}(),ClarkY{n%100,p}())
@pure British(c::C,t::T) where {C<:Profile{P},T<:Profile{P}} where P = British{C,T,2P-2}(c,t)

upper(n::British,c=1,x0=0) = upper(getprofiles(n,c,x0)...)
lower(n::British,c=1,x0=0) = lower(getprofiles(n,c,x0)...)
upperlower(n::British,c=1,x0=0) = upperlower(getprofiles(n,c,x0)...)
getprofiles(n::British,c=1,x0=0) = interval(n.c,c,x0),c*profile(n.c),chord(n.c),c*profile(n.t)

# Rotated thickness by camber slope

"""
    American{C,T,P} <: Airfoil

Thickness measured perpendicular to the camber line, American convention.
"""
struct American{C,T,P} <: Airfoil{P}
    c::C
    t::T
    @pure American{C,T,P}(c,t) where {C,T,P} = new{C,T,P}(c,t)
end

@pure American{n,p}() where {n,p} = American(Camber{Int(n÷100),p}(),ClarkY{n%100,p}())
@pure American(c::C,t::T) where {C<:Profile{P},T<:Profile{P}} where P = American{C,T,2P-2}(c,t)

upper(n::American,c=1,x0=0) = upper(getprofiles(n,c,x0)...)
lower(n::American,c=1,x0=0) = lower(getprofiles(n,c,x0)...)
upperlower(n::American,c=1,x0=0) = upperlower(getprofiles(n,c,x0)...)
getprofiles(n::American,c=1,x0=0) = interval(n.c,c,x0),c*profile(n.c),profileslope(n.c),c*profile(n.t)

# NACA specification selection

macro NACA_str(n)
    p = 150
    n5 = match(r"(\d{2}[01])(\d{2}(?>\.\d+)?)(?:-)?(?(?<=-)(\d{2}(?>\.\d+)?))?",n)
    n4 = match(r"(\d{2})(\d{2}(?>\.\d+)?)(?:-)?(?(?<=-)(\d{2}(?>\.\d+)?))?",n)
    n16 = match(r"1(\d)(?>-)(?:\((?=\d(?>\.\d+)?\)))?(\d(?>\.\d+)?)(?:(?<=\d)\))?(?:\((?=\d{2}(?>\.\d+)?\)))?(\d{2}(?>\.\d+)?)(?:(?<=\d)\))?",n)
    n6A = match(r"6(\d)(?>A)(?:\((?=\d(?>\.\d+)?\)))?(\d(?>\.\d+)?)(?:(?<=\d)\))?(?:\((?=\d{2}(?>\.\d+)?\)))?(\d{2}(?>\.\d+)?)(?:(?<=\d)\))?",n)
    if !isnothing(n5)
        British(NACA5{Meta.parse(n5[1]),p}(),isnothing(n5[3]) ? ClarkY{Meta.parse(n5[2]),p}() : Modified{Meta.parse(n5[2]),Meta.parse(n5[3]),p}())
    elseif !isnothing(n4)
        American(NACA4{Meta.parse(n4[1]),p}(),isnothing(n4[3]) ? ClarkY{Meta.parse(n4[2]),p}() : Modified{Meta.parse(n4[2]),Meta.parse(n4[3]),p}())
    elseif !isnothing(n16)
        American(NACA6{Meta.parse(n16[2]),p}(),Modified{Meta.parse(n16[1]),Meta.parse(n16[3]),p}())
    elseif !isnothing(n6A)
         American(NACA6A{Meta.parse(n6A[2]),p}(),Modified{Meta.parse(n6A[1]),Meta.parse(n6A[3]),p}())
    else
        throw(error("not valid"))
    end
end

# Joukowski

struct Joukowski{R,f,g,b,p} <: Airfoil{p}
    @pure Joukowski{R,f,g,b,p}() where {R,f,g,b,p} = new{R,f,g,b,p}()
end

joukowski(R,f,g,b,p) = Joukowski{R,f,g,b,p}()

function Base.complex(::Joukowski{R,f,g,b,p}) where {R,f,g,b,p}
    z = R*cis.(interval(2p-1,2π)).-(f-g*im)
    z .+ b^2*inv.(z)
end
