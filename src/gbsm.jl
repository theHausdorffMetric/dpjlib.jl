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

