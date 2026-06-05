using Test, FastGaussQuadrature, dn2

@testset "Gauss-Legendre implementation" begin
    ns = [1, 2, 3, 5, 7, 10, 17, 25, 50, 100]

    tol = 1e-10

    for n in ns
        x, w = gausslegendre(n)
        my_x, my_w = my_gausslegendre(n)

        @test isapprox(x, my_x; rtol=tol)
        @test isapprox(w, my_w; rtol=tol)
    end
end
