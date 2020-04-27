function v_asian_bm_sample(X₀, K, r, σ, m, n)
    Δ = 1/365
    X = X₀  + σ*√(m*Δ)*randn()
    x̂ = zero(X)
    for i in 1:n
        X += σ*√Δ*randn()
        x̂ += X
    end
    exp(-r*(m+n)*Δ)*max(x̂/n - K, 0)
end

function asian_bm(F, K, r, sigma, m, n, nruns)
        xm = v_asian_bm_sample(F,K,r,sigma,m,n)
        S = zero(typeof(xm))
     for i in 2:nruns
         x = v_asian_bm_sample(F,K,r,sigma,m,n)
         delta = x - xm
         xm = xm + delta/i
         S = S + delta*(x-xm)
     end
     v = S/nruns
     sigma_xm = sqrt(v/nruns)
     (xm, sigma_xm)
end

