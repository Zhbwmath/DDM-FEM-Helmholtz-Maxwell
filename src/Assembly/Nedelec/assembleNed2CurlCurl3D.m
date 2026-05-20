function A = assembleNed2CurlCurl3D(node, elem)
% ASSEMBLENED2CURLCURL3D  Assemble the conforming 3D NE_2 curl-curl matrix.

A = assembleNed2Core3D(node, elem, 'curl');
end
