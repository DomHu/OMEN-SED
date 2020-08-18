function [b,fval,exitflag,output] = fzero_vec(FunFcn,x1, x2, options,varargin)
% Vectorized single-variable nonlinear zero finding based on fzero
%
%  NB: unlike fzero, there is no error and exit if a solution fails (so check exitflag ..)
%  
%   X = fzero_vec(FUN,X1, X2), where X1,X2 is
%   finite interval where the sign of FUN(X1) differs from the sign of 
%   FUN(X2). An error occurs if this is not true.  Calling FZERO with a
%   finite interval guarantees FZERO will return a value near a point where
%   FUN changes sign.
%
%   exitflag(i):
%     1  FZERO found a zero X(i).
%    -5  FZERO may have converged to a singular point.
%    -6  FZERO can not detect a change in sign of the function.

if nargin > 3 && isfield(options,'Display') && ~isempty(options.Display)
    switch options.Display
        case {'notify','notify-detailed'}
            trace = 1;
        case {'none', 'off'}
            trace = 0;
        case {'iter','iter-detailed'}
            trace = 3;
        case {'final','final-detailed'}
            trace = 2;
        otherwise
            trace = 1;
    end
else
    trace = 0;
end

if nargin > 3 && isfield(options,'TolX') && ~isempty(options.TolX)
    tol = options.TolX;
else
    tol = eps;
end

% Initialization
fcount = 0;
iter = 0;

% Interval input
a = x1;
b = x2;

fa = FunFcn(a,varargin{:});
fb = FunFcn(b,varargin{:});

npts = max(length(fa),length(fb));

% expand any scalar-valued inputs
a = a.*ones(1,npts);
b = b.*ones(1,npts);
fa = fa.*ones(1,npts);
fb = fb.*ones(1,npts);

savea = a; saveb = b;
savefa = fa; savefb = fb;



fval = NaN(1,npts);

cmplt = false(1,npts);
exitflag = ones(1,npts);

fcount = fcount + 2;

% Check supplied interval for sign change or root at endpoints
for i=1:npts
    if ( fa(i) == 0 )
        b(i) = a(i);
        
        
        fval(i) = fa(i);
        cmplt(i) = true;
       
    elseif ( fb(i) == 0)
        % b = b;
        
        fval(i) = fb(i);
        cmplt(i) = true;
     
    elseif (fa(i) > 0) == (fb(i) > 0)
        cmplt(i) = true;
        exitflag(i) = -6;
    end
    
end

fc = fb;
%  allocate temporary variables
c = NaN(1,npts);
d = NaN(1,npts);
e = NaN(1,npts);

if trace > 1
    fprintf('\nfzero_vec npts %g\n',npts);
end

% Main loop
while any(~cmplt)
    for i=find(~cmplt)
        if fb(i) == 0 || a(i) == b(i)
            fval(i)= fb(i);
            cmplt(i) = true;
           
        else
            % Insure that b is the best result so far, a is the previous
            % value of b, and c is on the opposite side of the zero from b.
            if (fb(i) > 0) == (fc(i) > 0)
                c(i) = a(i);  fc(i) = fa(i);
                d(i) = b(i) - a(i);  e(i) = d(i);
            end
            if abs(fc(i)) < abs(fb(i))
                a(i) = b(i);    b(i) = c(i);    c(i) = a(i);
                fa(i) = fb(i);  fb(i) = fc(i);  fc(i) = fa(i);
            end
            
            % Convergence test and possible exit
            m = 0.5*(c(i) - b(i));
            toler = 2.0*tol*max(abs(b(i)),1.0);
            if (abs(m) <= toler) || (fb(i) == 0.0)
                fval(i) = fb(i);
                cmplt(i) = true;
           
            else
                
                if trace > 2
                    fprintf('i=%4g  fcount %5.0f   x=%13.6g f(x)=%13.6g\n',i,fcount, b(i), fb(i));
                end
                % Choose bisection or interpolation
                if (abs(e(i)) < toler) || (abs(fa(i)) <= abs(fb(i)))
                    % Bisection
                    d(i) = m;  e(i) = m;
                    % procedure='bisection';
                else
                    % Interpolation
                    s = fb(i)/fa(i);
                    if (a(i) == c(i))
                        % Linear interpolation
                        p = 2.0*m*s;
                        q = 1.0 - s;
                    else
                        % Inverse quadratic interpolation
                        q = fa(i)/fc(i);
                        r = fb(i)/fc(i);
                        p = s*(2.0*m*q*(q - r) - (b(i) - a(i))*(r - 1.0));
                        q = (q - 1.0)*(r - 1.0)*(s - 1.0);
                    end;
                    if p > 0, q = -q; else p = -p; end;
                    % Is interpolated point acceptable
                    if (2.0*p < 3.0*m*q - abs(toler*q)) && (p < abs(0.5*e(i)*q))
                        e(i) = d(i);  d(i) = p/q;
                        %procedure='interpolation';
                    else
                        d(i) = m;  e(i) = m;
                        % procedure='bisection';
                    end;
                end % Interpolation
                
                % Next point
                a(i) = b(i);
                fa(i) = fb(i);
                if abs(d(i)) > toler, b(i) = b(i) + d(i);
                elseif b(i) > c(i), b(i) = b(i) - toler;
                else b(i) = b(i) + toler;
                end
            end
           
            
        end
    end
    if any(~cmplt)
        fb = FunFcn(b,varargin{:});
        fcount = fcount + 1;
        iter = iter + 1;
    end
end % Main loop

%fval = fb; % b is the best value


for i=1:npts
    switch exitflag(i)
        case 1
            if abs(fval(i)) <= max(abs(savefa(i)),abs(savefb(i)))
                if trace > 1
                    fprintf('i=%4g zero found x=%13.6g f(x)=%13.6g\n',i, b(i), fb(i));
                end
            else
                exitflag(i) = -5;                
                if trace > 1
                    fprintf('i=%4g current value x(i)=%g near singular point %g %g \n',i,fb(i),savea(i),saveb(i));
                end              
            end
        case -6
            if trace > 1
                fprintf('i=%4g no sign change at points %g %g \n',i,savea(i),saveb(i));
            end
        otherwise
            error('unexpected exitflag(%g) = %g',i,exitflag(i));
    end
end


output.nfail        = sum(exitflag ~=  1);
output.nsingular    = sum(exitflag == -5);
output.nsignerr     = sum(exitflag == -6);
output.fcount       = fcount;
output.iter         = iter;

if trace > 1
    fprintf('fzero_vec all complete nfail %g nsingular %g nsignerr %g func-count %g\n\n', ... 
        output.nfail, output.nsingular, output.nsignerr, fcount);
end

end

