function [t] = system_time(lambda, mu, delta)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    function peak_aoi: compute the system time PMF for a tandem queue    %
%                                                                         %
%              author - Federico Chiariotti <fchi@es.aau.dk>              %
%                                                                         %
% lambda (1*1): arrival rate                                              %
% mu (1*2): service rates (do not set mu(2)-mu(1)=lambda, undef. case)    %
% delta (1*N): support vector of system time distribution                 %
%                                                                         %
% t (1*N): system time PMF                                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Auxiliary parameters

% Check special case (same service rate)
if (mu(1) == mu(2))
    % System time is an Erlang distribution
    t = lambda ^ 2 .* delta .* exp(-lambda .* delta);
else
    alpha = mu - lambda;
    rho = lambda ./ mu;
    exp_am1 = exp(-alpha(1) .* delta) - exp(-mu(1) .* delta);
    exp_am2 = exp(-alpha(2) .* delta) - exp(-mu(2) .* delta);
    exp_m12 = exp(-mu(1) .* delta) - exp(-mu(2) .* delta);
    exp_a12 = exp(-alpha(1) .* delta) - exp(-alpha(2) .* delta);
    exp_m1a2 = exp(-mu(1) .* delta) - exp(-alpha(2) .* delta);
    exp_a1m2 = exp(-alpha(1) .* delta) - exp(-mu(2) .* delta);

    % Case A
    mult = mu(2) * alpha(2) / mu(1);
    t_a = mult * alpha(1) / (mu(2) - mu(1)) .* exp_a12;
    t_a = t_a + mult * lambda / mu(2) .* exp(-alpha(1) .* delta) .* (1 - exp(-mu(2) .* delta));
    t_a = t_a - mult * mu(1) / (mu(2) - alpha(1)) .* exp_a1m2;

    % Case B
    t_b = lambda .* exp(-alpha(1) .* delta) .* (1 - exp(-mu(2) .* delta)) + alpha(1) * mu(2) / (mu(2) - mu(1)) .* exp_am1;
    t_b = t_b - lambda * alpha(2) * mu(2) / (mu(2) - mu(1)) / (mu(2) - alpha(1)) .* exp_a1m2;

    % Case C
    t_c = alpha(2) / mu(1) .* exp(-alpha(2) .* delta) .* (alpha(1) - mu(1) .* exp(-lambda .* delta) + lambda .* exp(-mu(1) .* delta));

    % Case D
    t_d = (mu(2) * (mu(1) - lambda) .* exp(-mu(1) .* delta) - mu(1) * (mu(2) - lambda) .* exp(-mu(2) .* delta)) / (mu(2) - mu(1));
    t_d = t_d + lambda .* exp(-(mu(2) + alpha(1)) .* delta);

    % General case
    t = t_a + t_b + t_c + t_d;
end

end

