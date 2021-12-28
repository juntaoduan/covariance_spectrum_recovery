function res = earthmover(u_pos,u_mas,v_pos,v_mas)
    cap = max(u_pos(end),v_pos(end));
    u_pos = [u_pos cap];
    v_pos = [v_pos cap];
    u_mas = [u_mas 0];
    v_mas = [v_mas 0];
    j=1;
    i=1;
    last = 0;
    G=0;
    dist = 0;
    while(i<=size(u_pos,2)&&j<=size(v_pos,2))
        if(v_pos(j)<u_pos(i))
            dist = dist + abs(G)*(v_pos(j)-last);
            G = G-v_mas(j);
            last = v_pos(j);
            j=j+1;
        else
            dist = dist + abs(G)*(u_pos(i)-last);
            G = G+u_mas(i);
            last = u_pos(i);
            i=i+1;
        end
        
    end
    res = dist;
end
