module ProjectiveGeometry

using Grassmann, PyPlot, LightGraphs, GraphPlot, Compose
import Grassmann: vector
import PyPlot: streamplot

export streamplot, graph

@inline field_components(M,n=2) = ([[M[x,y][k] for x ∈ 1:length(M[:,1]), y ∈ 1:length(M[1,:])] for k ∈ 1:n]...,)

function streamplot(t::T,x=collect(-5:0.1:5),y=x) where T<:TensorAlgebra{V} where V
    streamplot(x,y,field_components([vector(t*(x*Λ(V).v1 - y*Λ(V).v2)) for x ∈ x, y ∈ y])...)
end

graph(x,n="simplex.png",a=16cm,b=16cm) = draw(PNG(n,a,b), gplot(SimpleDiGraph(x),layout=circular_layout,nodelabel=collect(1:grade(vectorspace(x)))))

end # module
