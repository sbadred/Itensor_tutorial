
@testset "ITensor quantum chemical Hamiltonian" begin
    MOLECULE = "h4"
    δ = 1e-7
     data = "../data/fcidump_files/FCIDUMP." * MOLECULE
    Vnn, nsites,N, h, v = read_electron_integral_tensors(data)
    
    #ITensor
    #Define how many sweeps for algo
    nsweeps=2
    sweeps = ITensors.Sweeps(nsweeps)
    #cutoff two-site compression
    ITensors.cutoff!(sweeps, 1e-12)

    @info "Without particle number conservation" 
    sites = ITensors.siteinds("Electron", nsites , conserve_nf=false)

    @info "Creation of random mps"
    MPS_rank=16
    random_mps = ITensors.randomMPS(sites,MPS_rank)

    @info "Creation of mpo"
    Mpo =ITensors_MPO(sites,h, v,tol=δ)

    @timev energy_ITensors,Ψ=ITensors.dmrg(Mpo, random_mps,sweeps,eigsolve_krylovdim=2);


    @info "With particle number conservation"  
    sites = ITensors.siteinds("Electron", nsites , conserve_nf=true)
        
    #1: 0 particle #2-3: 1 particle #4 2particles
    state=[4,4,1,1]

    @info "Creation of random mps"
    random_mps = ITensors.randomMPS(sites,state,MPS_rank)
    println("Particle number: ", ITensors.flux(random_mps))

    @info "Creation of mpo"
    Mpo =ITensors_MPO(sites,h, v,tol=δ)

    @info "Ground_state search"
    ITensors.dmrg(Mpo, random_mps,sweeps,eigsolve_krylovdim=20);
 
end


