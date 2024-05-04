classdef mypolynomial

properties (Access = public)
    poly
end % properties


methods (Access = public)

function p = mypolynomial(a)
    if (nargin > 1)
        print_usage();
    end

    if (nargin == 0)
        p.poly = 0;
    else
        % Is this a copy constructor?
        if ( isa(a, "mypolynomial") )
            p = a;
        elseif (isreal(a) && isvector(a))
            p.poly = a(:).';  % force row vector
            disp("I am called")
        else
            error("@mypolynomial: A must be a real vector");
        end
    end
end

function disp (p)
    a = p.poly;
    first = true;
    for i = 1:length(a);
        if (a(i) != 0)
            if (first)
                first = false;
            elseif (a(i) > 0 || isnan (a(i)))
                printf (" +");
            end
            if (a(i) < 0)
                printf (" -");
            end
            if (i == 1)
                printf (" %.5g", abs (a(i)));
            elseif (abs (a(i)) != 1)
                printf (" %.5g *", abs (a(i)));
            end
            if (i > 1)
                printf (" X");
            end
            if (i > 2)
                printf (" ^ %d", i - 1);
            end
        end
    end
  
    if (first)
        printf(" 0");
    end
    printf("\n");
  
end


end % methods

end % classdef