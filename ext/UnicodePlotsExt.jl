module UnicodePlotsExt

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

using Grassmann, Cartan, FlowGeometry
isdefined(FlowGeometry, :Requires) ? (import FlowGeometry: UnicodePlots) : (using UnicodePlots)

UnicodePlots.lineplot(N::FlowGeometry.Profile,args...) = UnicodePlots.lineplot(profile(N),args...)
UnicodePlots.lineplot!(p,N::FlowGeometry.Profile,args...) = UnicodePlots.lineplot!(p,profile(N),args...)
function UnicodePlots.lineplot(N::Airfoil,args...)
    UnicodePlots.lineplot!(UnicodePlots.lineplot(complex(N),args...),N.c,args...)
end
function UnicodePlots.lineplot!(p,N::Airfoil,args...)
    UnicodePlots.lineplot!(UnicodePlots.lineplot!(p,complex(N),args...),N.c,args...)
end
function UnicodePlots.lineplot(N::DoubleArc,args...)
    U,L = upperlower(N)
    p = UnicodePlots.lineplot(U,args...)
    length(U)==length(L) && UnicodePlots.lineplot!(p,fiber(real.(U)),fiber(imag.(U).+imag.(L))./2,args...)
    UnicodePlots.lineplot!(p,L,args...)
end
function UnicodePlots.lineplot!(p,N::DoubleArc,args...)
    U,L = upperlower(N)
    UnicodePlots.lineplot!(p,U,args...)
    length(U)==length(L) && UnicodePlots.lineplot!(p,fiber(real.(U)),fiber(imag.(U).+imag.(L))./2,args...)
    UnicodePlots.lineplot!(p,L,args...)
end
Base.display(t::Airfoil) = (display(typeof(t)); display(UnicodePlots.lineplot(t)))
Base.display(t::FlowGeometry.Profile) = (display(typeof(t)); display(UnicodePlots.lineplot(t)))

end # module
