function K_si=tc_kernel_SI_bs(x,le,ri)
As=x(1);
Au=x(2);
sqcs=x(3);
sqc0=x(4);
sqcu=x(5);
sq_lambda=x(6);
sigma=x(7);
m=le+ri+1;
K_si=zeros(m,m);
g2=sq_lambda^2;
g2u2=sq_lambda^2*Au^2;
s =(-le:ri)';
t= -le:ri;
mat_t_s=triu(t-s,1); 
mat_tps=((t+s)+abs(t+s))/2; 
max_t_s2= triu((2*max(1,t+1)-t-s),1);
% mat_t_max_s = (abs(triu(t-max(1,s+1),1))+triu(t-max(1,s+1),1))/2;
            part1= sqcs^2*(g2*As.^(mat_tps)-As^2*sq_lambda.^(2*s).*(As.^mat_t_s))/(As^2-g2);
            K_si(le+3:le+1+ri,le+3:le+1+ri) =  K_si(le+3:le+1+ri,le+3:le+1+ri)+part1(le+3:le+1+ri,le+3:le+1+ri); 
            part2=  sqc0*sqcs*As.^mat_t_s.*sq_lambda.^(2*s);
            K_si(le+2:le+1+ri,le+2:le+1+ri) =  K_si(le+2:le+1+ri,le+2:le+1+ri)+part2(le+2:le+1+ri,le+2:le+1+ri); 
%             for ss=-le:ri
%                 for tt= max(ss+1,2):ri
% 
%                     smax=max(ss+1,1);
%                     vec = (lams.^(tt- (smax:(tt-1)))).*((lamu).^((smax:(tt-1))-ss)).*((gamma^2).^(smax:(tt-1)));
%                      part3(le+1+ss,le+1+tt)=   sqcu*sqcs*sum(vec);
%                      part3(le+1+ss,le+1+tt) =  sqcu*sqcs*lams^(-ss)*((g2s)^(max(ss+1,1))*lamu^(tt-max(ss+1,1))-g2s^tt)*lamu/(lamu-g2s);
%                 end
%             end

            vec=zeros(1,ri-1);
            mat_su= As.^(1:ri)'* (Au*(sq_lambda^2)).^(ri:-1:1);
            for i = 1: ri-1
                vec(i)=sum( sum(triu(mat_su,i)-triu(mat_su,i+1)));
            end
            vec = flip(vec);
            part3_m= zeros(le+1+ri);
            for ss=-le:ri
                for tt= max(ss+2,2):ri
                    part3_m(le+1+ss,le+1+tt) =  sqcu*sqcs * vec(tt-max(ss+1,1))*sq_lambda^(2*max(ss,0));
                end
            end
            part3= part3_m.* Au.^(max(le:-1:-ri,0))';
            
            K_si(:,le+3:le+1+ri) =  K_si(:,le+3:le+1+ri)+part3(:,le+3:le+1+ri); 
            part4 = sqc0*sqcu*Au.^mat_t_s.*sq_lambda.^(2*t);
            K_si(:,le+2:le+1+ri) =  K_si(:,le+2:le+1+ri)+part4(:,le+2:le+1+ri); 
            part5 = sqcu^2*sq_lambda.^(2*max(1,t+1)).*Au.^(max_t_s2)/(1-g2u2);
            K_si =  K_si+part5;
            K_si= triu(K_si,1)+triu(K_si,1)';
            
            part1_diag= diag(sqcs^2*(As.^(2*t).*g2-As^2*g2.^t)/(As^2-g2));
            K_si(le+3:le+1+ri,le+3:le+1+ri) =  K_si(le+3:le+1+ri,le+3:le+1+ri)+part1_diag(le+3:le+1+ri,le+3:le+1+ri); 
            part2_diag= diag(sqc0^2*sq_lambda.^(2*t));
            K_si(le+2:le+1+ri,le+2:le+1+ri) =  K_si(le+2:le+1+ri,le+2:le+1+ri)+part2_diag(le+2:le+1+ri,le+2:le+1+ri); 
            part3_diag= diag(sqcu^2*sq_lambda.^(2*max(1,t+1)).*Au.^(2*max(1-t,1))/(1-g2u2));
            K_si =  K_si+part3_diag; 
            
K_si=sigma^2*K_si;
end
