function c = u1inv(p, risk_aver,u)
    if p.quad_b > 0
        c = (1 - u) ./ p.quad_b;
    elseif p.exp_a > 0
        c = -log(u) ./ p.exp_a;
    else
        c = u .^ (-1 ./ risk_aver);
    end
end