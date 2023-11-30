using Pkg

PACKAGES = ["HTTP", "UUIDs", "JSON", "SparseArrays", "LinearAlgebra", "Printf", "JuliaFormatter", "Logging", "Test"]

for pkg in PACKAGES
    if !haskey(Pkg.installed(), Symbol(pkg))
        Pkg.add(pkg)
    end
end
