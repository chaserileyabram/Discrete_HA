function u = utility1(p, risk_aver,c)
	% first derivative of utility function
    if p.quad_b > 0
        u = 1 - p.quad_b .* c;
    elseif p.exp_a > 0
        u = exp(-p.exp_a .* c);
    else
        u = c.^(-risk_aver);
    end
end