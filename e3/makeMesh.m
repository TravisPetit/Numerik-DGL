function [F, P] = makeMesh(C, H, G, B)
    l = size(G, 2);
    m = size(B, 2);

    P = [G.'; B.'];

    F = zeros(m+l, 5);
    n = 0;
    for kk=1:l
        le = C(1, kk);  ri = C(2, kk);  do = C(3, kk);  up = C(4, kk);
        
        % internal quads, anchored top-right
        if le > 0 && do > 0
            dl = C(1, do);
            if dl > 0
                n = n+1;
                F(n, 1) = kk;  F(n, 2) = le;  F(n, 3) = dl;  F(n, 4) = do;  F(n, 5) = NaN;
            end
        end
        
        % edge quads
        if le > 0
            ld = C(3, le);
            if do < 0 && ld < 0
                n = n+1;
                F(n, 1) = kk;  F(n, 2) = le;  F(n, 3) = l-ld;  F(n, 4) = l-do;  F(n, 5) = NaN;
            end
        end
        if do > 0
            dr = C(2, do);
            if ri < 0 && dr < 0
                n = n+1;
                F(n, 1) = kk;  F(n, 2) = do;  F(n, 3) = l-dr;  F(n, 4) = l-ri;  F(n, 5) = NaN;
            end
        end
        if ri > 0
            ru = C(4, ri);
            if up < 0 && ru < 0
                n = n+1;
                F(n, 1) = kk;  F(n, 2) = ri;  F(n, 3) = l-ru;  F(n, 4) = l-up;  F(n, 5) = NaN;
            end
        end
        if up > 0
            ul = C(1, up);
            if le < 0 && ul < 0
                n = n+1;
                F(n, 1) = kk;  F(n, 2) = up;  F(n, 3) = l-ul;  F(n, 4) = l-le;  F(n, 5) = NaN;
            end
        end

        % edge triangles
        if le < 0 && do < 0
            n = n+1;
            F(n, 1) = kk;  F(n, 2) = l-le;  F(n, 3) = l-do;  F(n, 4) = NaN;  F(n, 5) = NaN;
        end
        if do < 0 && ri < 0
            n = n+1;
            F(n, 1) = kk;  F(n, 2) = l-do;  F(n, 3) = l-ri;  F(n, 4) = NaN;  F(n, 5) = NaN;
        end
        if ri < 0 && up < 0
            n = n+1;
            F(n, 1) = kk;  F(n, 2) = l-ri;  F(n, 3) = l-up;  F(n, 4) = NaN;  F(n, 5) = NaN;
        end
        if up < 0 && le < 0
            n = n+1;
            F(n, 1) = kk;  F(n, 2) = l-up;  F(n, 3) = l-le;  F(n, 4) = NaN;  F(n, 5) = NaN;
        end

        % edge pentangles
        if le > 0 && do > 0
            ld = C(3, le);
            dl = C(1, do);
            if ld < 0
                n = n+1;
                F(n, 1) = kk;  F(n, 2) = le;  F(n, 3) = l-ld;  F(n, 4) = l-dl;  F(n, 5) = do;
            end
        end
        if do > 0 && ri > 0
            dr = C(2, do);
            rd = C(3, ri);
            if dr < 0
                n = n+1;
                F(n, 1) = kk;  F(n, 2) = do;  F(n, 3) = l-dr;  F(n, 4) = l-rd;  F(n, 5) = ri;
            end
        end
        if ri > 0 && up > 0
            ru = C(4, ri);
            ur = C(2, up);
            if ru < 0
                n = n+1;
                F(n, 1) = kk;  F(n, 2) = ri;  F(n, 3) = l-ru;  F(n, 4) = l-ur;  F(n, 5) = up;
            end
        end
        if up > 0 && le > 0
            ul = C(1, up);
            lu = C(4, le);
            if ul < 0
                n = n+1;
                F(n, 1) = kk;  F(n, 2) = up;  F(n, 3) = l-ul;  F(n, 4) = l-lu;  F(n, 5) = le;
            end
        end
    end
    F = F(1:n, :);
end