function phi = gyroid_phi_at_point(p, L)
% Evaluates the standard gyroid implicit function at point p for period L.
% p : [x y z] in mm
% L : cell size (mm)

k = 2*pi / L;
phi = sin(k*p(1))*cos(k*p(2)) + ...
      sin(k*p(2))*cos(k*p(3)) + ...
      sin(k*p(3))*cos(k*p(1));
end
