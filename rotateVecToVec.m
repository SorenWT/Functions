function R = rotateVecToVec(a, b)
% R = rotateVecToVec(a, b)
% Returns a 3x3 rotation matrix R such that R*a is aligned with b.
% a and b must be 3-element vectors.

    % Ensure column vectors
    a = a(:); 
    b = b(:);

    % Normalize
    a = a / norm(a);
    b = b / norm(b);

    v = cross(a, b);
    c = dot(a, b);

    % If vectors are already aligned
    if norm(v) < 1e-12
        if c > 0
            R = eye(3);
        else
            % 180-degree rotation: pick orthogonal axis
            % Find vector orthogonal to a
            [~, idx] = min(abs(a));
            tmp = zeros(3,1);
            tmp(idx) = 1;
            v = cross(a, tmp);
            v = v / norm(v);
            R = -eye(3) + 2*(v*v');
        end
        return
    end

    % Skew-symmetric cross-product matrix
    vx = [   0   -v(3)  v(2);
           v(3)    0   -v(1);
          -v(2)  v(1)    0  ];

    % Rodrigues' formula
    R = eye(3) + vx + vx^2 * ((1 - c) / (norm(v)^2));
end