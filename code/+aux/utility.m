function u = utility(p, risk_aver,c)
	if isequal(size(risk_aver),size(c))
		u = zeros(size(risk_aver));
		Ilog_u = (risk_aver == 1);

		u(Ilog_u) = log(c(Ilog_u));
		u(~Ilog_u) = (c(~Ilog_u) .^(1-risk_aver(~Ilog_u) )-1)./(1-risk_aver(~Ilog_u) );
	elseif numel(risk_aver) == 1
        if p.quad_b > 0
            u = c - (p.quad_b/2) .* c.^2;
        elseif p.exp_a > 0
            u = -(1/p.exp_a) .* exp(-p.exp_a .* c);
        elseif risk_aver == 1
			u = log(c);
		else
			u = (c .^(1-risk_aver) - 1) ./ (1 - risk_aver);
        end
	end
end