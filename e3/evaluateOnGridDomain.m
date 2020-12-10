%% evaluateOnGridDomain
function fv = evaluateOnGridDomain(fh, G)
[~, nd] = size(G);
fv = zeros(nd, 1);

for i = 1:nd
    fv(i) = fh(G(:,i));
end

end