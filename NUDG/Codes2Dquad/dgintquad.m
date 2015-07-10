%Domain integrate the scalar-field 'field', return a scalar value
%"result"

%input: 'field': scalar field to be integrated (NpxK)
%           'V': Local Vandermonde matrix      (NpxNp)
%           'J': Jacobian scaling defined on each element. (NpxK)
%output: 'result': some scalar value (1x1)

function result =  dgintquad(field,V,J)
    Po=2;  %rectangle.
    %Po=sqrt(2); %scaling factor from legendre basis (triangles)
    %get modes of our field on each triangle
    modes = V\field;
    %modes = invV*field; %or switch to do this, could be faster.
    modezero = modes(1,:); %get mode 0 on each triangle.
    
    result = Po*(modezero*J(1,:)'); %dot product over all triangles to get integral
end

%note: only works for constant-jacobians
%think can generalize by V\(integrand.*J) at line 13.