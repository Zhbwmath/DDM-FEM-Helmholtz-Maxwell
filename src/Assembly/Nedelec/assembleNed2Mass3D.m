function M = assembleNed2Mass3D(node, elem)
% ASSEMBLENED2MASS3D  Assemble the conforming 3D NE_2 vector mass matrix.

M = assembleNed2Core3D(node, elem, 'mass');
end
