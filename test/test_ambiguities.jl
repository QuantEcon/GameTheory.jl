# ------------------------------ #
# Testing for method ambiguities #
# ------------------------------ #

@testset "Method ambiguities" begin
    ambiguities = Test.detect_ambiguities(GameTheory; recursive=true)
    for (m1, m2) in ambiguities
        @error "Ambiguous methods" m1 m2
    end
    @test isempty(ambiguities)
end

@testset "Exports are defined" begin
    undefined = [n for n in names(GameTheory) if !isdefined(GameTheory, n)]
    @test isempty(undefined)
end
