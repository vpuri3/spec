function [ue,ve,pe] = walsh_ex(x,y,nu,t)
    %exact solution of Walsh problem (Fig 1 '92 paper)
    lam = 25;
    ue = exp(-  lam*nu*t)   * (sin(5*y)+     cos(3*x).*cos(4*y));
    ve = exp(-  lam*nu*t)   * (cos(5*x)+0.75*sin(3*x).*sin(4*y));   
    pe = exp(-2*lam*nu*t)/4 * (cos(7*x)+0.75*cos(7*y)); % not sure about 0.75
end
