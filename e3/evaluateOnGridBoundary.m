function gv = evaluateOnGridBoundary(gh, B)
[~, nb] = size(B);
gv = zeros(nb, 1);

for i = 1:nb
    gv(i) = gh(B(:,i));
end

end