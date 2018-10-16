function r = pcepoly2DIP(name,varargin)

name = lower(name);
switch(name)
    case 'legendre'
        r = Legendre2DIP(varargin{:});
    case 'hermite'
        r = Hermite2DIP(varargin{:});
    case 'laguerre'
        r = Laguerre2DIP(varargin{:});
    case 'jacobi'
        r = Jacobi2DIP(varargin{:});
    otherwise
        error('Invalid argument.');
end