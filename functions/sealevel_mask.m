function flags = sealevel_mask(md, gl_buffer)
    phi = md.mask.ocean_levelset(:);
    x = md.mesh.x(:);
    y = md.mesh.y(:);
    tri = md.mesh.elements;

    floating_ocean = phi < 0;

    %construct mesh edges
    edges = [
        tri(:,[1 2]);
        tri(:,[2 3]);
        tri(:,[3 1])
        ];

    edges = unique(sort(edges,2), 'rows');

    %identify edges crossed by the grounding line zero contour
    phi1 = phi(edges(:,1));
    phi2 = phi(edges(:,2));

    crosses_gl = (phi1 < 0 & phi2 >= 0) | (phi2 < 0 & phi1 >= 0);

    gl_edges = edges(crosses_gl, :);

    %interpolate grounding-line crossing points
    p1 = phi(gl_edges(:,1));
    p2 = phi(gl_edges(:,2));

    alpha = -p1 ./ (p2-p1);

    xgl = x(gl_edges(:,1)) + alpha .* (x(gl_edges(:,2))-x(gl_edges(:,1)));
    ygl = y(gl_edges(:,1)) + alpha .* (y(gl_edges(:,2))-y(gl_edges(:,1)));

    distance_to_gl = inf(md.mesh.numberofvertices,1);

    for i = 1:md.mesh.numberofvertices
        distance_to_gl(i) = min(hypot(x(i)-xgl, y(i)-ygl));
    end

    flags = floating_ocean & distance_to_gl >= gl_buffer;

end
