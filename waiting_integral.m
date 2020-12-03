function [total] = waiting_integral(x, D, lambda, beta)

total = 0;
if (beta == 0)
    total = waiting_md1(D, lambda, x) - (1 - lambda * D);
else
    % k = 0
    if (beta == -lambda)
        total = total + lambda * x;
    else
        total = total + (exp((beta + lambda) * x) - 1) * lambda / (beta + lambda);
    end
    % k > 0
    for k = 1 : floor(x / D)
        z = x - k * D;
        if (z > 0)
            if (beta == -lambda)
                total = total + exp(beta * k * D) * (-lambda) ^ k * (lambda * z ^ (k + 1) / factorial(k + 1) + z ^ k / factorial(k));
            else
                for j = 0 : k - 1
                    total = total + beta * exp(beta * k * D) * (lambda) ^ k * (-1) ^ (j + 1) * z ^ j * exp((lambda + beta) * z) / factorial(j) / (lambda + beta) ^ (k - j + 1);
                end
                total = total - exp(beta * k * D) * (-lambda) ^ (k + 1) * z ^ (k) * exp((lambda + beta) * z) / factorial(k) / (lambda + beta);
                total = total + beta * lambda ^ k / ((beta + lambda) ^ (k + 1)) * exp(beta * k * D);
            end
        end
    end
    total = total * (1 - lambda * D);
end


