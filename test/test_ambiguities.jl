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
