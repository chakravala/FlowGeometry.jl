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

spheremesh() = TetGen.tetrahedralize(cubesphere(),"vpq1.414a0.1";holes=[TetGen.Point(0.0,0.0,0.0)])
cubemesh(hmax=0.1) = TetGen.tetrahedralize(∂(MiniQhull.delaunay(cube(2))), "vpq1.414a$hmax")

