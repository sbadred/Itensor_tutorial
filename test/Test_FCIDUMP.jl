#""" Test generation of one and two electron integrals from pyscf"""
#using Test
include("../src/FCIDUMP.jl")

#h2
@testset "Generate FCIDUMP FILE" begin
 molecule="N 0 0 0; N 0 0 1.1"
 basis="sto-3g"
 name="FCIDUMP.N2"
 cd("../data/fcidump_files/")
 Generate_FCIDUMP(molecule, basis,name)

 Vnn,sites,N,h,v=read_electron_integral_tensors(name)
 cd("../../")
end


#h6
@testset "Generate FCIDUMP FILE" begin
    molecule="H 0 0 0; H 0 0 1.1; H 0 0 -1.1; H 0 0 2.1; H 0 0 -2.1; H 0 0 3.1"
    basis="sto-3g"
    name="FCIDUMP.h6"
    cd("data/fcidump_files/")
    Generate_FCIDUMP(molecule, basis,name)
   
    Vnn,sites,N,h,v=read_electron_integral_tensors(name)
    cd("../../")
end


#h8
@testset "Generate FCIDUMP FILE" begin
    molecule="H 0 0 0; H 0 0 1.1; H 0 0 -1.1; H 0 0 2.1; H 0 0 -2.1; H 0 0 3.1; H 0 0 -3.1; H 0 0 4.1"
    basis="sto-3g"
    name="FCIDUMP.h8"
    cd("data/fcidump_files/")
    Generate_FCIDUMP(molecule, basis,name)
   
    Vnn,sites,N,h,v=read_electron_integral_tensors(name)
    cd("../../")
end


#h10
@testset "Generate FCIDUMP FILE" begin
    molecule="H 0 0 0; H 0 0 1.1; H 0 0 -1.1; H 0 0 2.1; H 0 0 -2.1; H 0 0 3.1; H 0 0 -3.1; H 0 0 4.1; H 0 0 -4.1; H 0 0 5.1"
    basis="sto-3g"
    name="FCIDUMP.h10"
    cd("data/fcidump_files/")
    Generate_FCIDUMP(molecule, basis,name)
   
    Vnn,sites,N,h,v=read_electron_integral_tensors(name)
    cd("../../")
end


#h12
@testset "Generate FCIDUMP FILE" begin
    molecule="H 0 0 0; H 0 0 1.1; H 0 0 -1.1; H 0 0 2.1; H 0 0 -2.1; H 0 0 3.1; H 0 0 -3.1; H 0 0 4.1; H 0 0 -4.1; H 0 0 5.1
              ; H 0 0 -5.1; H 0 0 6.1"
    basis="sto-3g"
    name="FCIDUMP.h12"
    cd("data/fcidump_files/")
    Generate_FCIDUMP(molecule, basis,name)
   
    Vnn,sites,N,h,v=read_electron_integral_tensors(name)
    cd("../../")
end


#h14
@testset "Generate FCIDUMP FILE" begin
    molecule="H 0 0 0; H 0 0 1.1; H 0 0 -1.1; H 0 0 2.1; H 0 0 -2.1; H 0 0 3.1; H 0 0 -3.1; H 0 0 4.1; H 0 0 -4.1; H 0 0 5.1
              ; H 0 0 -5.1; H 0 0 6.1; H 0 0 -6.1; H 0 0 7.1"
    basis="sto-3g"
    name="FCIDUMP.h14"
    cd("data/fcidump_files/")
    Generate_FCIDUMP(molecule, basis,name)
   
    Vnn,sites,N,h,v=read_electron_integral_tensors(name)
    cd("../../")
end


#h16
@testset "Generate FCIDUMP FILE" begin
    molecule="H 0 0 0; H 0 0 1.1; H 0 0 -1.1; H 0 0 2.1; H 0 0 -2.1; H 0 0 3.1; H 0 0 -3.1; H 0 0 4.1; H 0 0 -4.1; H 0 0 5.1
              ; H 0 0 -5.1; H 0 0 6.1; H 0 0 -6.1; H 0 0 7.1; H 0 0 -7.1; H 0 0 8.1"
    basis="sto-3g"
    name="FCIDUMP.h16"
    cd("data/fcidump_files/")
    Generate_FCIDUMP(molecule, basis,name)
   
    Vnn,sites,N,h,v=read_electron_integral_tensors(name)
    cd("../../")
end


#h18
@testset "Generate FCIDUMP FILE" begin
    molecule="H 0 0 0; H 0 0 1.1; H 0 0 -1.1; H 0 0 2.1; H 0 0 -2.1; H 0 0 3.1; H 0 0 -3.1; H 0 0 4.1; H 0 0 -4.1; H 0 0 5.1
              ; H 0 0 -5.1; H 0 0 6.1; H 0 0 -6.1; H 0 0 7.1; H 0 0 -7.1; H 0 0 8.1; H 0 0 -8.1; H 0 0 9.1"
    basis="sto-3g"
    name="FCIDUMP.h18"
    cd("data/fcidump_files/")
    Generate_FCIDUMP(molecule, basis,name)
   
    Vnn,sites,N,h,v=read_electron_integral_tensors(name)
    cd("../../")
end


#h20
@testset "Generate FCIDUMP FILE" begin
    molecule="H 0 0 0; H 0 0 1.1; H 0 0 -1.1; H 0 0 2.1; H 0 0 -2.1; H 0 0 3.1; H 0 0 -3.1; H 0 0 4.1; H 0 0 -4.1; H 0 0 5.1
              ; H 0 0 -5.1; H 0 0 6.1; H 0 0 -6.1; H 0 0 7.1; H 0 0 -7.1; H 0 0 8.1; H 0 0 -8.1; H 0 0 9.1; H 0 0 -9.1; H 0 0 10.1"
    basis="sto-3g"
    name="FCIDUMP.h20"
    cd("data/fcidump_files/")
    Generate_FCIDUMP(molecule, basis,name)
   
    Vnn,sites,N,h,v=read_electron_integral_tensors(name)
    cd("../../")
end

@testset "Generate FCIDUMP FILE" begin
    molecule="O 0.0000 0.0000 0.1173; H 0.0000 0.7572 -0.4692; H 0.0000 -0.7572	 -0.4692    "
    basis="sto-3g"
    name="FCIDUMP.h2oo2"
    cd("data/fcidump_files/")
    Generate_FCIDUMP(molecule, basis,name)
   
    Vnn,sites,N,h,v=read_electron_integral_tensors(name)
    cd("../../")
end