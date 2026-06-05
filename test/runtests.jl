using dn1, Test

@testset "RedkaMatrika" begin
    # Matrix:
    # [10  0  2]
    # [ 0  5  0]
    # [ 3  0  8]
    A = RedkaMatrika(
        [[10.0, 2.0], [5.0], [3.0, 8.0]],
        [[1, 3], [2], [1, 3]]
    )

    @testset "size" begin
        @test size(A) == (3, 3)
        @test size(A, 1) == 3
        @test size(A, 2) == 3
    end

    @testset "getindex" begin
        # Non-zero entries
        @test A[1,1] == 10.0
        @test A[1,3] == 2.0
        @test A[2,2] == 5.0
        @test A[3,1] == 3.0
        @test A[3,3] == 8.0

        # Zeros
        @test A[1,2] == 0.0
        @test A[2,1] == 0.0
        @test A[2,3] == 0.0
    end

    @testset "setindex!" begin
        # Modify existing entry
        A[1,1] = 20.0
        @test A[1,1] == 20.0

        # Add new non-zero entry
        A[2,3] = 7.0
        @test A[2,3] == 7.0

        # Remove entry by setting to zero
        A[3,1] = 0.0
        @test A[3,1] == 0.0

        # Other values unchanged
        @test A[3,3] == 8.0
    end

    @testset "Multiplication with vector" begin
        # Matrix
        # [20  0  2]
        # [ 0  5  7]
        # [ 0  0  8]
        x = [1.0, 2.0, 3.0]
        y = A * x
        @test y ≈ [
            20*1 + 2*3,
            5*2 + 7*3,
            8*3
        ]
    end

    @testset "SOR" begin
        A = RedkaMatrika(
            [
                [4.0, -1.0],
                [-1.0, 4.0, -1.0],
                [-1.0, 3.0]
            ],
            [
                [1, 2],
                [1, 2, 3],
                [2, 3]
            ]
        )

        b = [15.0, 10.0, 10.0]
        x0 = zeros(3)

        x, it = sor(A, b, x0, 1.25)

        @test it > 0
        @test x ≈ [5.0, 5.0, 5.0] atol=1e-8
    end
end # testset