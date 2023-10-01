
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
    #Mpo =ITensors_MPO(sites,h, v,tol=δ)
    Mpo =ITensor_MPO(sites,h, v,tol=δ)

    @timev energy_ITensors,Ψ=ITensors.dmrg(Mpo, random_mps,sweeps,eigsolve_krylovdim=2);


    @info "With particle number conservation"  
    sites = ITensors.siteinds("Electron", nsites , conserve_nf=true)
        
    #create a Hartree fock state
    occ_state = [i in 1:N/2 ? 1 : 0 for i in 1:nsites]
    occ_to_state = Dict([0 => 1, 1 => 4])
    hartree_fock_state = [occ_to_state[n] for n in occ_state]
  

    @info "Creation of random mps"
    random_mps = ITensors.randomMPS(sites,hartree_fock_state,MPS_rank)
    println("Particle number: ", ITensors.flux(random_mps))

    @info "Creation of mpo"
    Mpo =ITensors_MPO(sites,h, v,tol=δ)

    @info "Ground_state search"
    ITensors.dmrg(Mpo, random_mps,sweeps,eigsolve_krylovdim=20);
 
end


