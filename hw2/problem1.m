a = 1.0;
eps=1;
while true
    if((a==1.0+eps))
        disp(eps)
        break
    end
    eps=eps/2;
end