@test_throws ArgumentError OrbitalIndex(-1, -2)
@test_throws ArgumentError OrbitalIndex(3, -1)
@test_throws ArgumentError OrbitalIndex(3, 5)

@testset ">> Operators" begin
    @test OrbitalIndex(0, 0) == OrbitalIndex(0, 0)
    @test OrbitalIndex(2, 1) == OrbitalIndex(2, 1)
    @test OrbitalIndex(2, 1) != OrbitalIndex(2, 0)
    @test OrbitalIndex(2, 1) ≠ OrbitalIndex(2, 0)

    @test OrbitalIndex(2, 1) > OrbitalIndex(2, 0)
    @test OrbitalIndex(2, 1) > OrbitalIndex(1, 0)
    @test OrbitalIndex(2, 1) < OrbitalIndex(3, 0)
end

@testset ">> Linear indexing for first numbers" begin
    index = 1
    for n in 0:100
        for l in 0:n
            orb = OrbitalIndex(n, l)
            @test convert(Integer, orb) == index
            @test convert(OrbitalIndex, convert(Integer, orb)) == orb
            index += 1
        end
    end
end

@testset ">> conversions for n=$n" for n in rand(0:1000000, 5)
    @testset ">> conversions for l=$l" for l in rand(0:n, 5)
        orb = OrbitalIndex(n, l)
        back  = convert(OrbitalIndex, convert(Integer, orb))
        @test typeof(back) === typeof(orb)
        @test back == orb
    end
end

@testset ">> String indexing" begin
    @test nl"2, 2" == OrbitalIndex(2, 2)
	@test nl"3, 2" == OrbitalIndex(3, 2)
end

@testset ">> Printing" begin
	buffer = IOBuffer()
	print(buffer, nl"12, 4")
	seekstart(buffer)
	@test readstring(buffer) == "n=12, l=4"
end

@testset ">> looping by 1" begin
    index = 3
    for u in nl"1,1":nl"5, 2"
        @test convert(Integer, u) == index
        index += 1
    end
    @test convert(Integer, nl"5, 3") == index
end

@testset ">> looping by (n, 0)" begin
    n = 1
    index = 0
    for u in nl"1,1":nl"2, 0":nl"5, 2"
        @test u.l == 1
        @test u.n == n
        index = convert(Integer, u)
        n += 2
    end
    @test index == convert(Integer, nl"5, 1")
end

@testset ">> looping by (n, l)" begin
    expected = [OrbitalIndex(1, 1), OrbitalIndex(4, 3), OrbitalIndex(7, 5)]
    actual = [u for u in nl"1,1":nl"3, 2":nl"9, 2"]
    @test expected == actual
end
