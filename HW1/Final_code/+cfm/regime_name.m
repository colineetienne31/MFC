function name = regime_name(Minf)
%REGIME_NAME Return the flow regime used by this assignment.
if abs(Minf-0.3)<1e-8
    name='subsonic';
elseif abs(Minf-0.91)<1e-8
    name='transonic';
elseif abs(Minf-2.2)<1e-8
    name='supersonic';
else
    if Minf<0.8
        name='subsonic';
    elseif Minf<1.2
        name='transonic';
    else
        name='supersonic';
    end
end
end
