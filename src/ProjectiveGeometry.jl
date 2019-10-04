module ProjectiveGeometry

using Requires, Compose, Grassmann
import Grassmann: vector

@inline field_components(M,n=2) = ([[M[x,y][k] for x ∈ 1:length(M[:,1]), y ∈ 1:length(M[1,:])] for k ∈ 1:n]...,)

function __init__()
    @require PyPlot="d330b81b-6aea-500a-939a-2ce795aea3ee" begin
        export streamplot
        function PyPlot.streamplot(t::T,x=collect(-5:0.1:5),y=x) where T<:TensorAlgebra{V} where V
            PyPlot.streamplot(x,y,field_components([vector(t*(x*Λ(V).v1 - y*Λ(V).v2)) for x ∈ x, y ∈ y])...)
        end
    end
    @require GraphPlot="a2cc645c-3eea-5389-862e-a155d0052231" begin
        using LightGraphs
        export graph
        graph(x,n="simplex.png",a=16cm,b=16cm) = draw(PNG(n,a,b), gplot(SimpleDiGraph(x),layout=circular_layout,nodelabel=collect(1:grade(vectorspace(x)))))
    end
end

end # module
