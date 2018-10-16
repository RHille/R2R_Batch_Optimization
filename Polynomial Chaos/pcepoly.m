function r = pcepoly(name,varargin)

name = lower(name);
switch(name)
    case 'legendre'
        r = LegendrePoly(varargin{:});
    case 'hermite'
        r = HermitePoly(varargin{:});
    case 'laguerre'
        r = LaguerrePoly(varargin{:});
    case 'jacobi'
        r = JacobiPoly(varargin{:});
    otherwise
        error('Invalid argument.');
end