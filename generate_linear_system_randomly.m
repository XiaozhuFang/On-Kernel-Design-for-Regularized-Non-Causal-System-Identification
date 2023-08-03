function [system, system_new]= generate_linear_system_randomly(order,poleub, f_times, option)
% Randomly generate a stable non-minimum phase lienar system
%
%   system = generate_linear_system_randomly(order,poleub ,Num, option)
%   input: 
%   order: the order of the system; 
%   poleub: the upper bound of the pole's magnitude option:
%     -1 replace the two poles with two real poles (or two complex pairs)
%     at [0.8 0.9] and [1.1, 1.2]; 
%     -2 replace the two poles with two real
%     poles at [0.8 0.9] and [1.1, 1.2];
%
%   Author: XiaozhuFang 20220618.  
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('............generating system start...........\n')
pole_max = 2; % pole check
value = 1;
while value > 1e8 || pole_max>poleub
    bw = inf;
    % step 1: generating a random system
    while isinf(bw)||isnan(bw) % finite bandwidth of generated system

        mc= rss(order,1,1);
        bw = bandwidth(mc);
        if mc.d==0 || isinf(bw)||isnan(bw)
            bw=inf;
            continue
        end
        f = f_times*bw*2*pi; %100 times of bandwidth
        try
        md = c2d(mc,1/f,'zoh');
        catch
            bw=inf;
            continue
        end
        if ~all((abs(zero(md))<0.96)+(abs(zero(md))>1.04))
            bw=inf;
            continue
        end
    end
    system = idpoly(md);
    if system.B(1)==0
        pole_max=2;
        continue;
    end
    % step 2: Relocating the domiant zero near the unit circle
    if option <=2
    system_new= replacepole(system, order, option);
    end
    % step 3: reject some systems
    pole_max = max(abs(pole(system))); % reject unstable system

    Polepart = system.F;
        Zeropart = system.B;
        Polepart(end+1:numel(Zeropart))=0;
        Zeropart(end+1:numel(Polepart))=0;
        [r,~,~]= residue(Polepart,Zeropart);

    value = max(abs([system.f system.b]));   % reject system having large coefficient for numerical issue

    if sum(isnan(r))>0
        pole_max=2;
        continue
    end
    if option <=2
        pole_max = max(abs(pole(system_new))); % reject unstable system

        Polepart = system_new.F;
        Zeropart = system_new.B;
        Polepart(end+1:numel(Zeropart))=0;
        Zeropart(end+1:numel(Polepart))=0;
        [r,~,~]= residue(Polepart,Zeropart);

        value_new = max(abs([system_new.f system_new.b]));   % reject system having large coefficient for numerical issue
        value =max(value, value_new);
        if sum(isnan(r))>0
            pole_max=2;
            continue
        end
    end
end
end

function system = replacepole(system, order, option)
root_z = [];
ang = 1.7; % this angle constraint the angle of compelx pole in [-1,1]
if option==1
    root_zero =(1.1+0.1*rand()).*exp(1j*rand()*ang) ;
    root_z = [root_z root_zero conj(root_zero)];
    root_zero = (0.9-0.1*rand()).*exp(1j*rand()*ang);
    root_z = [root_z root_zero conj(root_zero)];
elseif option == 2
    root_zero = 1.1+0.1*rand() ;
    root_z = [root_z root_zero];
    root_zero = 0.9-0.1*rand();
    root_z = [root_z root_zero];
end
rootsys=roots(system.B);
%             roots(polyz*system_temp.B(1))
root_z= [ root_z  rootsys'];
root_z=root_z(1:order);
if imag(root_z(end))~=0 && real(root_z(end))~= real(root_z(end-1))
    root_z(end)=abs(root_z(end));
end
polyz=poly(root_z);
system.B =polyz*system.B(1);
% we relocate zeros, but keep the original poles
end