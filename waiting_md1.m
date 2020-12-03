function [pw] = waiting_md1(T, lambda, x)

m = floor(x / T);
pw = 0;
for k = 0 : m
    pw = pw + (1 - lambda * T) * (-lambda * (x - k * T)) ^ k / factorial(k) * exp(lambda * (x - k * T));
end