@testitem "Consistency tests" begin
    using PDBTools
    using MolSimToolkit.PackmolInputCreator
    test_dir = PackmolInputCreator.PackmolInputCreatorDirectory*"/test"

    # system with water only, with constant density
    system = SolutionBoxUSC(
        solute_pdbfile = "$test_dir/data/poly_h.pdb",
        solvent_pdbfile = "$test_dir/data/water.pdb",
        cossolvent_pdbfile = "$test_dir/data/water.pdb",
        density_table = hcat(0:0.1:1, ones(11))
    )
    mw = 55.508250191225926
    for x in (0.0, 0.2, 0.5, 0.7, 1.0)
        @test convert_concentration(system, x, "x" => "vv") ≈ x
        @test convert_concentration(system, x, "x" => "mol/L"; density = 1.0) ≈ x * mw atol = 1e-3
        @test convert_concentration(system, x, "vv" => "x") ≈ x 
        @test convert_concentration(system, x, "vv" => "mol/L"; density = 1.0) ≈ x * mw
        @test convert_concentration(system, x * mw, "mol/L" => "x"; density = 1.0) ≈ x 
        @test convert_concentration(system, x * mw, "mol/L" => "vv"; density = 1.0) ≈ x
    end

    # system with ideal solution 
    system = SolutionBoxUSC(
        solute_pdbfile = "$test_dir/data/poly_h.pdb",
        solvent_pdbfile = "$test_dir/data/water.pdb",
        cossolvent_pdbfile = "$test_dir/data/water.pdb",
        density_table = hcat([0.0 + 0.05*i for i in 0:20], [1.0 + 0.05*i for i in 0:20])
    )
    @test system.solvent_molar_mass ≈ 18.01 atol = 0.01
    @test system.cossolvent_molar_mass ≈ 18.01 atol = 0.01
    @test density_pure_solvent(system) ≈ 1.0
    @test density_pure_cossolvent(system) ≈ 2.0

    # Concentration conversions in this ideal system, from the molar fraction, x
    ρc = density_pure_cossolvent(system) # g / mL
    ρw = density_pure_solvent(system) # g / mL
    M = system.solvent_molar_mass # g / mol
    ρ(x) = ρc*x + ρw*(1-x) # g / mL
    vv(x) = (x / ρc) / (x / ρc + (1-x) / ρw) 
    v(x) = M / (1000*ρ(x)) # L/mol: volume of 1 mol (c + w) of solution
    mx(x) = x / v(x) # molarity of cossolute

    # Concentration conversions
    for x in (0.0, 0.2, 0.5, 0.7, 1.0)
        @test convert_concentration(system, x, "x" => "vv") ≈ vv(x)
        @test convert_concentration(system, x, "x" => "mol/L"; density = ρ(x)) ≈ mx(x)
        @test convert_concentration(system, vv(x), "vv" => "x") ≈ x 
        @test convert_concentration(system, vv(x), "vv" => "mol/L"; density = ρ(x)) ≈ mx(x) 
        @test convert_concentration(system, mx(x), "mol/L" => "x"; density = ρ(x)) ≈ x 
        @test convert_concentration(system, mx(x), "mol/L" => "vv"; density = ρ(x)) ≈ vv(x) 
    end

end

@testitem "Write packmol input" begin

    using PDBTools
    using MolSimToolkit.PackmolInputCreator
    using DelimitedFiles
    test_dir = PackmolInputCreator.PackmolInputCreatorDirectory*"/test"

    # Ethanol-water mixture
    density_table = readdlm("$test_dir/data/water_ethanol.dat", comments=true, comment_char='#')
    system = SolutionBoxUSC(
        solute_pdbfile = "$test_dir/data/poly_h.pdb",
        solvent_pdbfile = "$test_dir/data/water.pdb",
        cossolvent_pdbfile = "$test_dir/data/ethanol.pdb",
        density_table = density_table,
    )
    Mc = 1000 * density_pure_cossolvent(system) / system.cossolvent_molar_mass # mol / L pure ethanol
    @test convert_concentration(system, 1.0, "x" => "mol/L"; density = density_pure_cossolvent(system)) ≈ Mc

    tmp_input_file = tempname()
    write_packmol_input(system; concentration = 0.5, cunit = "x", margin = 20.0, input = tmp_input_file)
    @test isfile(tmp_input_file)
    f1 = read("$test_dir/data/water_ethanol.inp", String)
    f2 = read(tmp_input_file, String)
    @test f1 == f2

end