function r = pcepoly2D(name,varargin)

name = lower(name);
switch(name)
    case 'legendre'
        r = LegendrePoly2D(varargin{:});
    case 'hermite'
        r = HermitePoly2D(varargin{:});
    case 'laguerre'
        r = LaguerrePoly2D(varargin{:});
    case 'jacobi'
        r = JacobiPoly2D(varargin{:});
    otherwise
        error('Invalid argument.');
end