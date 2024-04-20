using TestItemRunner
@run_package_tests

@testitem "Doctests" begin
    using Documenter: doctest
    doctest(MolSimToolkit)
end

@testitem "Aqua.test_all" begin
    import Aqua
    Aqua.test_all(
        MolSimToolkit;
        ambiguities = false
    )
end
