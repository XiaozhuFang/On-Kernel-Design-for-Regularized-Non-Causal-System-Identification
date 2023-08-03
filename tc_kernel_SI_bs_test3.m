function K=tc_kernel_SI_bs_test3(x,le,ri)
sigma=x(7);
m=le+ri+1;
K_si=zeros(m,m);
for s=-le:ri
    for t= s:ri
        if t>s && s>1 
            K_si(le+1+s,le+1+t)=K1(t,s,x)+K2(t,s,x)+K3(t,s,x)+K4(t,s,x)+K5(t,s,x);
        end
        if t>s && s==1
            K_si(le+1+s,le+1+t)=K1(t,s,x)+K2(t,s,x)+K3(t,s,x)+K4(t,s,x);  
        end
        if t>1 && 1>s
                 K_si(le+1+s,le+1+t)=K1(t,s,x)+K2(t,s,x)+K3(t,s,x); 
        end
        if t==1 && s<1
            K_si(le+1+s,le+1+t)=K1(t,s,x)+K2(t,s,x);
        end
        if t<1 && t>s
            K_si(le+1+s,le+1+t)=K1(t,s,x);
        end
    end
end
K_diag=zeros(m,m);
for t=-le:ri
        if t>1
            K_diag(le+1+t,le+1+t)=K1(t,t,x)+K5(t,t,x)+K6(t,t,x);
        end
        if t==1
            K_diag(le+1+t,le+1+t)=K1(t,t,x)+K6(t,t,x);
        end
        if t<1
            K_diag(le+1+t,le+1+t)=K1(t,t,x);
        end
end

K=sigma^2*( triu(K_si,1)+triu(K_si,1)'+K_diag);
end

function value=K1(t,s,x)
Ac=x(1);
Aa=x(2);
sqcc=x(3);
sqc0=x(4);
sqca=x(5);
sq_lambda=x(6);
sigma=x(7);
g2_u2=sq_lambda^2/(Ac^2);
g2s=sq_lambda^2*Aa;
s_u=Aa/Ac;
g2_u=sq_lambda^2/Ac;
g2s2=sq_lambda^2*Aa^2;
value= sqca^2*(sq_lambda^(2*max(1,t+1))*Aa^(max(2-s-t,2+t-s)))/(1-g2s2);
end

    function value=K2(t,s,x)
        Ac=x(1);
        Aa=x(2);
        sqcc=x(3);
        sqc0=x(4);
        sqca=x(5);
        sq_lambda=x(6);
        sigma=x(7);
        g2_u2=sq_lambda^2/(Ac^2);
        g2s=sq_lambda^2*Aa;
        s_u=Aa/Ac;
        g2_u=sq_lambda^2/Ac;
        g2s2=sq_lambda^2*Aa^2;
        value = sqc0*sqca*Aa^(t-s)*sq_lambda^(2*t);
    end
    function value=K3(t,s,x)
        Ac=x(1);
        Aa=x(2);
        sqcc=x(3);
        sqc0=x(4);
        sqca=x(5);
        sq_lambda=x(6);
        sigma=x(7);
        g2_u2=sq_lambda^2/(Ac^2);
        g2s=sq_lambda^2*Aa;
        s_u=Aa/Ac;
        g2_u=sq_lambda^2/Ac;
        g2s2=sq_lambda^2*Aa^2;
        value= sqcc*sqca*Aa^(-s)*((g2s)^(max(s+1,1))*Ac^(t-max(s+1,1))-g2s^t)*Ac/(Ac-g2s);
    end
    function value=K4(t,s,x)
        Ac=x(1);
        Aa=x(2);
        sqcc=x(3);
        sqc0=x(4);
        sqca=x(5);
        sq_lambda=x(6);
        sigma=x(7);
        g2_u2=sq_lambda^2/(Ac^2);
        g2s=sq_lambda^2*Aa;
        s_u=Aa/Ac;
        g2_u=sq_lambda^2/Ac;
        g2s2=sq_lambda^2*Aa^2;
        value= sqc0*sqcc*Ac^(t-s)*sq_lambda^(2*s);
    end
    function value=K5(t,s,x)
        Ac=x(1);
        Aa=x(2);
        sqcc=x(3);
        sqc0=x(4);
        sqca=x(5);
        gamma=x(6);
        sigma=x(7);
        g2_u2=gamma^2/(Ac^2);
        g2s=gamma^2*Aa;
        s_u=Aa/Ac;
        g2_u=gamma^2/Ac;
        g2s2=gamma^2*Aa^2;
        value = sqcc^2* (Ac^(t+s)*g2_u2-gamma^(2*s)*Ac^(t-s))/(1-g2_u2);
    end
    function value=K6(t,s,x)
        Ac=x(1);
        Aa=x(2);
        sqcc=x(3);
        sqc0=x(4);
        sqca=x(5);
        gamma=x(6);
        sigma=x(7);
        g2_u2=gamma^2/(Ac^2);
        g2s=gamma^2*Aa;
        s_u=Aa/Ac;
        g2_u=gamma^2/Ac;
        g2s2=gamma^2*Aa^2;
        value = sqc0^2*gamma^(2*t);
    end

