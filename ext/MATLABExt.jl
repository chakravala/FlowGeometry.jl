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

decsg(args...) = MATLAB.mxcall(:decsg,1,args...)
function decsg(N::Airfoil{pts}) where pts
    P = fiber(complex(N))[1:end-1]
    x1,x2,x3,x4 = -1.5,3.5,3.5,-1.5
    y1,y2,y3,y4 = 1.5,1.5,-1.5,-1.5
    R = [[3,4,x1,x2,x3,x4,y1,y2,y3,y4];zeros(2pts-8)]
    A = [[2,pts];[real.(P);imag.(P)]]
    decsg([R A],"R-A","RA")
end
export decsg

