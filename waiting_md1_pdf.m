function [pw] = waiting_md1_pdf(T, lambda, x)

if (x == 0)
    pw = 1 - lambda * T;
else
    m = floor(x / T);
    pw = 0;
    pw = lambda * (1 - lambda * T) * exp(lambda * x);
    
    for k = 1 : m
        if (x - k * T > 0)
            pw = pw + (1 - lambda * T) * (-lambda) ^ k * (x - k * T) ^ (k - 1) / factorial(k) * exp(lambda * (x - k * T)) * (k + lambda * (x - k * T));
        end
    end
end