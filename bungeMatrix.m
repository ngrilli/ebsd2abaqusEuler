% Nicolo Grilli
% University of Oxford
% AWE project 2019

% rotation matrix from Euler angles
% as needed in the crystal plasticity code
% use angles in radians

function [ret] = bungeMatrix(phi1,Phi,phi2)
	matrix = zeros(3,3);
	matrix(1,1) = cos(phi1) * cos(phi2) - sin(phi1) * sin(phi2) * cos(Phi);
	matrix(1,2) = - cos(phi1) * sin(phi2) - sin(phi1) * cos(phi2) * cos(Phi);
	matrix(1,3) = sin(phi1) * sin(Phi);
	matrix(2,1) = sin(phi1) * cos(phi2) + cos(phi1) * sin(phi2) * cos(Phi);
	matrix(2,2) = - sin(phi1) * sin(phi2) + cos(phi1) * cos(phi2) * cos(Phi);
	matrix(2,3) = - cos(phi1) * sin(Phi);
	matrix(3,1) = sin(phi2) * sin(Phi);
	matrix(3,2) = cos(phi2) * sin(Phi);
	matrix(3,3) = cos(Phi);
	ret = matrix;
end
