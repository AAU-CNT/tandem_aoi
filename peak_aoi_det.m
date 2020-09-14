function [aoi] = peak_aoi_det(lambda, mu, D, delta)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function peak_aoi: computes the PAoI PMF for a tandem M/M/1-M/D/1 queue %
%                                                                         %
%              author - Federico Chiariotti <fchi@es.aau.dk>              %
%                                                                         %
% lambda (1*1): arrival rate                                              %
% mu: service rate of the M/M/1 system                                    %
% D: service time of the M/D/1 system                                     %
% delta (1*N): support vector of PAoI distribution                        %
%                                                                         %
% aoi (1*N): PAoI PMF                                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Auxiliary parameters
alpha = mu - lambda;
rho = lambda / mu;


aoi_a = zeros(1, length(delta));
aoi_b = zeros(1, length(delta));
aoi_c = zeros(1, length(delta));
aoi_d = zeros(1, length(delta));


for i = 1 : length(delta)
    if (delta(i) > 2 * D)
        x = delta(i) - 2 * D;
        % Case A
        aoi_a(i) = alpha * exp(-alpha * x) * waiting_integral(x, D, lambda, alpha);
        aoi_a(i) = aoi_a(i) - alpha * exp(-mu * x) * waiting_integral(x, D, lambda, mu);
        aoi_a(i) = aoi_a(i) - alpha * exp(-mu * D) * exp(-alpha * x) * waiting_integral(x, D, lambda, -lambda);
        aoi_a(i) = aoi_a(i) + alpha * exp(-mu * (delta(i) - D)) * waiting_integral(x, D, lambda, 0);
        % Consider zero queuing
        aoi_a(i) = aoi_a(i) + (1 - lambda * D) * alpha * exp(-alpha * x) * (1 - exp(-mu * D)) * (1 - exp(-lambda * x));
        
        % Case C
        aoi_c(i) = aoi_c(i) + waiting_integral(x, D, lambda, mu) * alpha * exp(-mu * x);
        aoi_c(i) = aoi_c(i) - waiting_integral(x, D, lambda, alpha) * mu * exp(-mu * x) * exp(-lambda * D);
        aoi_c(i) = aoi_c(i) + exp(-mu * (x + D)) * lambda * waiting_integral(x, D, lambda, 0);
        % Consider zero queuing
        aoi_c(i) = aoi_c(i) + (1 - lambda * D) * (alpha * exp(-mu * x) * (1 - exp(-lambda * D)) - lambda * exp(-alpha * x - lambda * (x + D)) * (1 - exp(-alpha * D)));
        
        % Cases B and D
        for k = 0 : floor(delta(i) - 2 * D)
            z = delta(i) - (k + 2) * D;
            if (z >= 0)
                % Case B
                for j = 0 : k + 1
                    aoi_b(i) = aoi_b(i) - alpha * mu * (-lambda) ^ j * exp(-alpha * (delta(i) - D) - lambda * (k + 1) * D) * z ^ j / factorial(j) / lambda * (1 - lambda * D);
                    if (j <= k)
                        % Case D
                        aoi_d(i) = aoi_d(i) + exp(-mu * (delta(i) - D)) * mu * lambda * (1 - lambda * D) * exp(alpha * D * (k + 1)) * lambda ^ k / (mu ^ (k - j + 1)) * exp(mu * z) * z ^ j * (-1) ^ j / factorial(j);
                        aoi_d(i) = aoi_d(i) + exp(-mu * (delta(i) - D)) * mu * lambda * (1 - lambda * D) * (-lambda) ^ (j - 1) * exp(lambda * z) * z ^ j / factorial(j);
                    end
                end
                % Lower limit of the integral
                aoi_b(i) = aoi_b(i) + alpha * mu * exp(-mu * (delta(i) - D)) / lambda * (1 - lambda * D);
                aoi_d(i) = aoi_d(i) - exp(-mu * (delta(i) - D)) * mu * lambda * (1 - lambda * D) * exp(alpha * D * (k + 1)) * rho ^ k / mu;
                aoi_d(i) = aoi_d(i) + exp(-mu * (delta(i) - D)) * mu * (1 - lambda * D);
            end
        end
        
    end
end

% General case
aoi = aoi_a + aoi_b + aoi_c + aoi_d;

end

