function [aoi] = peak_aoi(lambda, mu, delta)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     function peak_aoi: computes the Peak AoI PMF for a tandem queue     %
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
aoi_b = aoi_b - lambda * mu(1) / (alpha(2)) .* exp(-mu(1) .* delta) .* (1 - exp(-alpha(2) .* delta));
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

% General case
aoi = aoi_a + aoi_b + aoi_c + aoi_d;

end

