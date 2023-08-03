% the upper bound of hyper is [1-tol 1-tol  inf inf inf 1-tol 1-tol inf inf 4.99]
% the lower bound of hyper is [-1+tol -1+tol -inf -inf -inf tol tol  -inf -inf 0.99]
c_sign = hyper(10);
sign_u = sign(2.5-c_sign);
sign_s = sign(abs(2.5-c_sign)-1);

switch kernel        
          case {'NCSI-FO'}
            a_c = hyper(1);
            a_a = hyper(2);
            sqrt_lambda_c = hyper(6);
            sqrt_lambda_a = hyper(7);
            sigma_c = exp(hyper(8));
            sigma_a = exp(hyper(9));
            c_c = sign_s*sqrt(exp(hyper(3)));
            c_0 = sqrt(exp(hyper(4)));
            c_a = sign_u*sqrt(exp(hyper(5)));
            Pi1 = tc_kernel_SI_bs([a_c;a_a;c_c;c_0;c_a;sqrt_lambda_c;sigma_c],ind(di),ind(di+1));
            Pi2 = tc_kernel_SI_bs([a_a;a_c;c_a;c_0;c_c;sqrt_lambda_a;sigma_a],ind(di),ind(di+1));
            Pi0 = tc_kernel_SI_delta([a_c;a_a;c_c;c_0;c_a],ind(di),ind(di+1),1);
            Pi = Pi0+Pi1+flip(flip(Pi2,1),2);
        case {'NCSI-DC'}
            rho_c = hyper(6);
            rho_a = hyper(7);
            sqrt_lambda_c = hyper(1);
            sqrt_lambda_a = hyper(2);
            c_c = sign_s*sqrt(exp(hyper(3)));
            c_0 = sqrt(exp(hyper(4)));
            c_a = sign_u*sqrt(exp(hyper(5)));
            Pi1 = tc_kernel_SI_bs([rho_c*sqrt_lambda_c;rho_a*sqrt_lambda_a;c_c;c_0;c_a;sqrt_lambda_c;sqrt(1-rho_c^2)],ind(di),ind(di+1));
            Pi2 = tc_kernel_SI_bs([rho_a*sqrt_lambda_a;rho_c*sqrt_lambda_c;c_a;c_0;c_c;sqrt_lambda_a;sqrt(1-rho_a^2)],ind(di),ind(di+1));
            Pi0 = tc_kernel_SI_delta([rho_c*sqrt_lambda_c;rho_a*sqrt_lambda_a;c_c;c_0;c_a],ind(di),ind(di+1),1);
            Pi = Pi0+Pi1+flip(flip(Pi2,1),2);
        case {'NCSI-TC'}
            lambda_c = hyper(1);
            lambda_a = hyper(2);
            c_c = sign_s*sqrt(exp(hyper(3)));
            c_0 = sqrt(exp(hyper(4)));
            c_a = sign_u*sqrt(exp(hyper(5)));
            Pi1 = tc_kernel_SI_bs([lambda_c;lambda_a;c_c;c_0;c_a;sqrt(lambda_c);sqrt(1-lambda_c)],ind(di),ind(di+1));
            Pi2 = tc_kernel_SI_bs([lambda_a;lambda_c;c_a;c_0;c_c;sqrt(lambda_a);sqrt(1-lambda_a)],ind(di),ind(di+1));
            Pi0 = tc_kernel_SI_delta([lambda_c;lambda_a;c_c;c_0;c_a],ind(di),ind(di+1),1);
            Pi = Pi0+Pi1+flip(flip(Pi2,1),2);

end
