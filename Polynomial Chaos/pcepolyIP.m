function r = pcepolyIP(name,varargin)

name = lower(name);
switch(name)
    case 'legendre'
        r = LegendreIP(varargin{:});
    case 'hermite'
        r = HermiteIP(varargin{:});
    case 'laguerre'
        r = LaguerreIP(varargin{:});
    case 'jacobi'
        r = JacobiIP(varargin{:});
    otherwise
        error('Invalid argument.');
end