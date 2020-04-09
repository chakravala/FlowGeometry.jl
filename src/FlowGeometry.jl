module FlowGeometry

#=using Requires, Compose, Grassmann
import Grassmann: vector

@inline field_components(M,n=2) = ([[M[x,y][k] for x ∈ 1:length(M[:,1]), y ∈ 1:length(M[1,:])] for k ∈ 1:n]...,)

function __init__()
    @require PyPlot="d330b81b-6aea-500a-939a-2ce795aea3ee" begin
        export streamplot
        function PyPlot.streamplot(t::T,x=collect(-5:0.1:5),y=x) where T<:TensorAlgebra{V} where V
            PyPlot.streamplot(x,y,field_components([vector(t*(x*Λ(V).v1 - y*Λ(V).v2)) for x ∈ x, y ∈ y])...)
        end
    end
end=#

function nacas(n,pts)
    ns=string(n,pad=4)
    Mi,Pi,Ti = parse.(Int,(ns[1],ns[2],ns[3:4]))
    a = SVector(0.2969,-0.1260,-0.3515,0.2843,-0.1015) # open trailing edge
    #-0.1036 # closed trailing edge
    M,P,T = Mi/100,Pi/10,Ti/100
    x = range(0,1,length=pts)
    # camber and gradient
    yc,dyc_dx,θ = zeros(pts),zeros(pts),zeros(pts)
    for i ∈ 1:pts
        if x[i] ≥ 0 && x[i] < P
            yc[i] = (M/P^2)*((2P*x[i])-x[i]^2)
            dyc_dx[i] = ((2M)/P^2)*(P-x[i])
        elseif x[i] ≥ P && x[i] ≤ 1
            yc[i] = (M/(1-P)^2)*(1-2P+2P*x[i]-x[i]^2)
            dyc_dx[i] = (2M/(1-P)^2)*(P-x[i])
        end
        θ[i] = atan(dyc_dx[i])
    end
    # thickness distribution
    yt = zeros(pts)
    for i ∈ 1:pts
        term0 = a[1]*sqrt(x[i])
        term1 = a[2]*x[i]
        term2 = a[3]*x[i]^2
        term3 = a[4]*x[i]^3
        term4 = a[5]*x[i]^4
        yt[i] = 5T*(term0+term1+term2+term3+term4)
    end
    # upper surface points
    xu,yu = zeros(pts),zeros(pts)
    for i ∈ 1:pts
        xu[i] = x[i] - yt[i]*sin(θ[i])
        yu[i] = yc[i] + yt[i]*cos(θ[i])
    end
    # lower surface points
    xl,yl = zeros(pts),zeros(pts)
    for i ∈ 1:pts
        xl[i] = x[i] + yt[i]*sin(θ[i])
        yl[i] = yc[i] - yt[i]*cos(θ[i])
    end
    xu[end] = 1
    xl[end] = 1
    yu[end] = 0
    yl[end] = 0
    rect = [[3,4,-1,2,2,-1,1,1,-1,-1];zeros(4pts-12)]
    [rect [[2,2pts-2];xu;reverse(xl)[2:end-1];yu;reverse(yl)[2:end-1]]],"R-A","RA"
end

end # module
