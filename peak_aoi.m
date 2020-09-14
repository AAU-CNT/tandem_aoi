function [aoi] = peak_aoi(lambda, mu, delta)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function peak_aoi: computes the PAoI PMF for a tandem M/M/1-M/M/1 queue %
%                                                                         %
%              author - Federico Chiariotti <fchi@es.aau.dk>              %
%                                                                         %
% lambda (1*1): arrival rate                                              %
% mu (1*2): service rates (do not set mu(2)-mu(1)=lambda, undef. case)    %
% delta (1*N): support vector of PAoI distribution                        %
%                                                                         %
% aoi (1*N): PAoI PMF                                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Auxiliary parameters
alpha = mu - lambda;
rho = lambda ./ mu;

% Check for special case (same service rate)
if (mu(1) == mu(2))
    alpha = alpha(1);
    mu = mu(1);
    
    % Case A
    aoi_a = (alpha ^ 2) / (lambda ^ 2) * mu .* exp(-mu .* delta) .* (1 + (lambda .* delta - 1) .* exp(lambda .* delta));
    aoi_a = aoi_a + exp(-mu .* delta) .* (mu .* (exp(lambda .* delta) - 1) + lambda .* (exp(-alpha .* delta) - exp(lambda .* delta)));
    aoi_a = aoi_a + mu ^ 2 / (lambda ^ 2) .* exp(-mu .* delta) .* (lambda .* delta - exp(lambda .* delta) + 1) * alpha;
    
    % Case B
    aoi_b = - (mu / lambda) ^ 2 / 2 * alpha .* exp(-mu .* delta) .* (lambda .* delta .* (lambda .* delta + 2) - 2 .* exp(lambda .* delta) + 2);
    aoi_b = aoi_b + mu ^ 2 / lambda .* exp(-mu .* delta) .* (1 + lambda .* delta - exp(lambda .* delta));
    aoi_b = aoi_b + mu / alpha .* exp(-2 * mu .* delta) .* (alpha .* exp((mu + lambda) .* delta) - mu .* exp(mu .* delta) + lambda .* exp(lambda .* delta));
    
    % Case C
    aoi_c = -mu ^ 2 .* exp(-mu .* delta) .* ((delta .^ 2) / 2 * alpha + delta .* (alpha - lambda) / lambda); 
    aoi_c = aoi_c + mu .* exp(-alpha .* delta) .* (1 - exp(-lambda .* delta)) * ((alpha / lambda) ^ 2);
    aoi_c = aoi_c - mu .* exp(-mu .* delta) * lambda / alpha .* (1 - exp(-alpha .* delta));
    
    % Case D
    aoi_d = mu ^ 2 * lambda .* exp(-mu .* delta) .* (2 .* cosh(alpha .* delta) - 2 - (alpha .* delta) .^ 2) / (alpha ^ 2);
else
    exp_am1 = exp(-alpha(1) .* delta) - exp(-mu(1) .* delta);
    exp_am2 = exp(-alpha(2) .* delta) - exp(-mu(2) .* delta);
    exp_m12 = exp(-mu(1) .* delta) - exp(-mu(2) .* delta);
    exp_a12 = exp(-alpha(1) .* delta) - exp(-alpha(2) .* delta);
    exp_m1a2 = exp(-mu(1) .* delta) - exp(-alpha(2) .* delta);
    exp_a1m2 = exp(-alpha(1) .* delta) - exp(-mu(2) .* delta);

    % Case A
    aoi_a = alpha(1) * alpha(2) * mu(2) / (mu(2) - mu(1)) / (mu(1) - alpha(2)) .* exp_m1a2;
    aoi_a = aoi_a + alpha(2) * mu(1) * mu(2) / (mu(2) - mu(1)) / (mu(2) - alpha(1)) .* exp_m12;
    aoi_a = aoi_a - lambda .* exp(-mu(1) .* delta) .* (1 - exp(-alpha(2) .* delta));
    aoi_a = aoi_a + alpha(1) * mu(1) * alpha(2) / (mu(2) - mu(1)) / (mu(2) - alpha(1)) .* exp_am1; 

    % Case B
    aoi_b = mu(1) * exp_am1;
    aoi_b = aoi_b - lambda / (alpha(2)) * mu(1) .* exp(-mu(1) .* delta) .* (1 - exp(-alpha(2) .* delta));
    aoi_b = aoi_b + (mu(1) * mu(2) * alpha(1)) / (lambda * (mu(2) - mu(1))) .* (exp(-alpha(1) .* delta) - (1 + lambda .* delta) .* exp(-mu(1) .* delta));
    aoi_b = aoi_b + (mu(1) * mu(2) * alpha(2)) / ((mu(2) - mu(1)) * (mu(2) - alpha(1))) * (lambda .* exp_m12 / (mu(2) - mu(1)));
    aoi_b = aoi_b - (mu(1) * mu(2) * alpha(2)) / ((mu(2) - mu(1)) * (mu(2) - alpha(1))) .* exp_am1;

    % Case C
    mult = mu(2) / (mu(2) - mu(1));
    aoi_c = alpha(2) * mu(1) .* delta .* exp(-mu(2) .* delta) + lambda .* exp(-mu(1) .* delta) .* (1 - exp(-alpha(2) .* delta));
    aoi_c = aoi_c + alpha(1) * alpha(2) / (alpha(2) - mu(1)) .* exp_m1a2;
    aoi_c = aoi_c - mu(1) * alpha(2) / (mu(2) - mu(1)) .* exp_m12 - alpha(1) * alpha(2) / lambda * exp_am2;
    aoi_c = aoi_c - lambda * alpha(2) / alpha(1) .* exp(-mu(2) .* delta) .* (1 - exp(-alpha(1) .* delta));
    aoi_c = aoi_c * mult;

    % Case D
    mult = mu(1) * mu(2) * lambda;
    aoi_d = (exp(-lambda * delta) + exp(-(mu(2) + alpha(1)) * delta)) ./ (alpha(1) * alpha(2));
    aoi_d = aoi_d + delta .* (exp(-mu(2) * delta) - exp(-mu(1) * delta)) / (mu(2) - mu(1)) ;
    aoi_d = aoi_d - (exp(-mu(2) .* delta) + exp(-mu(1) .* delta)) ./ (alpha(1) * alpha(2));
    aoi_d = aoi_d * mult;
end

% General case
aoi = aoi_a + aoi_b + aoi_c + aoi_d;

end

