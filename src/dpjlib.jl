module dpjlib

import Distributions
import Statistics


function GBSM_Option(FlagIsCall, S, X, T, r, b, σ)
 # b = r, give the BS stock option model
 # b = r - q, gives Merton 73, stock option with dividend
 # b = 0, give black 76
 # b = 0 && r = 0, gives Asay for margined futures options
 # b = r- r_f, gives garman kohlhagen currency options model
    d1 = (log(S / X) + (b + σ * σ * 0.5) * T) / (σ * sqrt(T))
    d2 = d1 - σ * sqrt(T)
    Price = 0.0
    if FlagIsCall
        Price = S * exp((b - r) * T) * cdf(Normal(), d1) -
                X * exp(-r * T) * cdf(Normal(), d2)
    else
        Price = -S * exp((b - r) * T) * cdf(Normal(), -d1) +
                X * exp(-r * T) * cdf(Normal(), -d2)
    end
    return(Price)
end

function GBSM_Delta(FlagIsCall, S, X, Time, r, b, σ)
  d1 = (log(S / X) + (b +  σ * σ / 2) * Time) / (σ * sqrt(Time))
    if FlagIsCall
            Price = exp((b - r) * Time) * cdf(Normal(),d1)
    else
            Price = exp((b - r) * Time) * (cdf(Normal(),d1) - 1)
    end
    return(Price)
 end

function mcgbsm(F,X,T,r,sigma,n)
        Y = exp.(randn(n)*sigma*sqrt(T).+(0-0.5*sigma^2)*T)*F
        mean(exp(-r*T) .* max.(Y.-X,0.0))
end

function mcgbsm(F,X,T,r,sigma,n)
        Y = exp.(randn(n)*sigma*sqrt(T).+(0-0.5*sigma^2)*T)*F
        mean(exp(-r*T) .* max.(Y.-X,0.0))
end

function v_asian_sample(X₀, K, r, σ, m, n)
    Δ = 1/365
    X = X₀ * exp((0-σ^2/2)*(m*Δ) + σ*√(m*Δ)*randn())
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

end # module
