function [pi] = steady_md1(T, lambda, n)

lambda = lambda / T;

if (n == 0)
    pi = 1 - lambda;
else
    if (n == 1)
        pi = (1 - lambda) * (exp(lambda) - 1);
    else
        pi = (1 - lambda) * exp(n * lambda);
        for k = 1 : n - 1
            pik = exp(k * lambda) * (-1) ^ (n - k) * ((k * lambda) ^ (n - k) / factorial(n - k) + (k * lambda) ^ (n - k - 1) / factorial(n - k - 1));
            pi = pi + (1 - lambda) * pik;
        end
    end
end

