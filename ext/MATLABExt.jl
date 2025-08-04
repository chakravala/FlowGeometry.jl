module MATLABExt

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

using FlowGeometry
isdefined(FlowGeometry, :Requires) ? (import FlowGeometry: MATLAB) : (using MATLAB)

FlowGeometry.decsg(args...) = MATLAB.mxcall(:decsg,1,args...)

end # module
