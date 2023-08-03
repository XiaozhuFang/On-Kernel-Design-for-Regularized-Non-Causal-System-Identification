function K_si=tc_kernel_SI_delta(x,le,ri,blk)
lams=x(1);
lamu=x(2);
sqcs=x(3);
sqc0=x(4);
sqcu=x(5);

m=le+ri+1;
vu = flip(sqcu*lamu.^(1:le));
v0 = sqc0;
vc = sqcs*lams.^(1:ri);
%sum( [vu,v0,vc])
if blk==0
    K_si=blkdiag(vu'*vu,v0'*v0,vc'*vc);
else
    K_si=[vu,v0,vc]'*[vu,v0,vc];    
end
end