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

Makie.lines(N::Profile,args...) = Makie.lines(profile(N),args...)
Makie.lines!(N::Profile,args...) = Makie.lines!(profile(N),args...)
function Makie.lines(N::Airfoil,args...)
    display(Makie.lines(N,args...))
    Makie.lines!(N.c,args...)
end
function Makie.lines!(N::Airfoil,args...)
    Makie.lines!(N,args...)
    Makie.lines!(N.c,args...)
end
Makie.lines(N::Airfoil,args...) = Makie.lines(complex(N),args...)
Makie.lines!(N::Airfoil,args...) = Makie.lines!(complex(N),args...)
function Makie.lines(N::DoubleArc,args...)
    U,L = upperlower(N)
    display(Makie.lines(U,args...))
    length(U)==length(L) && Makie.lines!(fiber(real.(U)),fiber(imag.(U).+imag.(L))./2,args...)
    Makie.lines!(L,args...)
end
function Makie.lines!(N::DoubleArc,args...)
    U,L = upperlower(N)
    Makie.lines!(U,args...)
    length(U)==length(L) && Makie.lines!(fiber(real.(U)),fiber(imag.(U).+imag.(L))./2,args...)
    Makie.lines!(L,args...)
end

