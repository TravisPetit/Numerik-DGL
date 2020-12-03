function plotMeshFunction(F, P, u, v)
    c = [u(:); v(:)];
    Pz = [P, c];
    patch('Faces', F, 'Vertices', Pz, 'FaceVertexCData', c, 'FaceColor', 'interp', 'EdgeColor', 'none');
    axis('equal');
    view(3);
end