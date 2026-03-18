using Pkg

println("Active project: ", Base.active_project())
Pkg.status()

import BifurcationKit
import OrdinaryDiffEq
import Plots
import Accessors

println("BifurcationKit loaded: ", BifurcationKit)
println("OrdinaryDiffEq loaded: ", OrdinaryDiffEq)
println("Plots loaded: ", Plots)
println("Accessors loaded: ", Accessors)
