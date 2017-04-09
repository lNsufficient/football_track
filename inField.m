function answer = inField(corners, X)
%inField Checks if a pixel is inside the field
%   Detailed explanation goes here

    c = corners;
    k1 = (c(1,2)-c(2,2))/(c(1,1)-c(2,1));
    m1 = c(1,2)-k1*c(1,1);
    l1 = k1*X(1) + m1;
    
    if l1>X(2);
        answer = 0;
        return
    end
    
    l2 = (c(2,2) + c(3,2))/2;
    
    if l2 > X(2)
        answer = 0;
        return
    end
    k3 = (c(3,2)-c(4,2))/(c(3,1)-c(4,1));
    m3 = c(3,2)-k3*c(3,1);
    l3 = k3*X(1) + m3;
    if l3 > X(2)
       answer = 0;
       return
    end
    
    l4 = (c(1,2) + c(4,2))/2;
    
    if l4 < X(2)
       answer = 0;
       return
    end
    
    answer = 1;
end

