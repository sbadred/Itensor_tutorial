

"""
    create_one_electron_mpo(sites::Array{ITensors.Index{Int64}, 1},
                            one_electron_tensor::Array{Float64, 2})

Create an ITensor.MPO that represents the hamiltonian for the one-electron
part of the molecular hamiltonian.
"""
function create_one_electron_mpo(sites,
    one_electron_tensor::Array{Float64,2},
    tol::Float64)

    site_count = length(sites)
    mpo = ITensors.AutoMPO()
    for c in CartesianIndices(size(one_electron_tensor))
        t = one_electron_tensor[c]
        mpo += t, "Cdagup", c[1], "Cup", c[2]
        mpo += t, "Cdagdn", c[1], "Cdn", c[2]
    end

    return ITensors.MPO(mpo, sites; cutoff=tol)
end
"""
    create_two_electron_mpo(sites::Array{ITensors.Index{Int64},1},
                            two_electron_tensor::Array{Float64, 4})

Create an ITensor.MPO that represents the hamiltonian for the two-electron
part of the molecular hamiltonian.
"""
function create_two_electron_mpo(sites,
    two_electron_tensor::Array{Float64,4},
    tol::Float64)

    mpo = ITensors.AutoMPO()
    for c in CartesianIndices(size(two_electron_tensor))
        u = 1 / 2 * two_electron_tensor[c]
        mpo += u, "Cdagup", c[1], "Cdagup", c[3], "Cup", c[4], "Cup", c[2]
        mpo += u, "Cdagup", c[1], "Cdagdn", c[3], "Cdn", c[4], "Cup", c[2]
        mpo += u, "Cdagdn", c[1], "Cdagup", c[3], "Cup", c[4], "Cdn", c[2]
        mpo += u, "Cdagdn", c[1], "Cdagdn", c[3], "Cdn", c[4], "Cdn", c[2]
    end
   
    return ITensors.MPO(mpo, sites; cutoff=tol)
end


function create_mpo_chemicalH(one_electron_tensor::Array{Float64,2}, two_electron_tensor::Array{Float64,4}; tol=1e-14)
    mpo = ITensors.OpSum()
    for c in CartesianIndices(size(one_electron_tensor))
        t = one_electron_tensor[c]
        if norm(t) > tol
            ITensors.add!(mpo, t, "Cdagup", c[1], "Cup", c[2])
            ITensors.add!(mpo, t, "Cdagdn", c[1], "Cdn", c[2])
        end
    end
    for c in CartesianIndices(size(two_electron_tensor))
        u = 1 / 2 * two_electron_tensor[c]
        if abs(u) > tol
            if (c[1] ≠ c[3]) && (c[4] ≠ c[2])
                ITensors.add!(mpo, u, "Cdagup", c[1], "Cdagup", c[3], "Cup", c[4], "Cup", c[2])
                ITensors.add!(mpo, u, "Cdagup", c[1], "Cdagdn", c[3], "Cdn", c[4], "Cup", c[2])
            end
            ITensors.add!(mpo, u, "Cdagdn", c[1], "Cdagup", c[3], "Cup", c[4], "Cdn", c[2])
            ITensors.add!(mpo, u, "Cdagdn", c[1], "Cdagdn", c[3], "Cdn", c[4], "Cdn", c[2])
        end
    end
    return mpo
end

function ITensor_MPO(sites,h, v; tol=1e-12)
    mpo = create_mpo_chemicalH(h, v)
    return ITensors.MPO(mpo, sites; cutoff=tol)
end

function ITensors_MPO(sites, h, v; tol=1e-12)
    mpo1 = create_one_electron_mpo(sites, h, 1e-10)
    mpo2 = create_two_electron_mpo(sites, v, 1e-10)
    return mpo1+mpo2
end



function ITensors_MPS(K)
    sites = ITensors.siteinds("Electron", K, conserve_qns=false)
    return ITensors.randomMPS(sites)
end