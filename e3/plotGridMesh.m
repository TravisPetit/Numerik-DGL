function plotGridMesh(C, H, G, B)
    [F, P] = makeMesh(C, H, G, B);
    patch('Faces', F, 'Vertices', P, 'FaceColor', [.5 .5 .5]);
    axis('equal');
end