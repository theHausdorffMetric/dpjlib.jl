function v_asian_sample(X₀, K, r, σ, m, n)
    Δ = 1/365
    X = X₀ * exp((r-σ^2/2)*(m*Δ) + σ*√(m*Δ)*randn())
    x̂ = zero(X)
    for i in 1:n
        X *= exp((r-σ^2/2)*Δ + σ*√Δ*randn())
        x̂ += X
    end
    exp(-r*(m+n)*Δ)*max(x̂/n - K, 0)
end

function asian_main(F,K, sigma, m,n, nruns)
        xm = v_asian_sample(F,K,0.01,sigma,m,n)
        S = zero(typeof(xm))
     for i in 2:nruns
         x = v_asian_sample(F,K,0.01,sigma,m,n)
         delta = x - xm
         xm = xm + delta/i
         S = S + delta*(x-xm)
     end
     v = S/nruns
     sigma_xm = sqrt(v/nruns)
     (xm, sigma_xm)
end

