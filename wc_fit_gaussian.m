function out=wc_fit_gaussian(x,bl,xmax,ymax,sx)
        out = bl + ymax*exp(-(x-xmax).^2/(2*sx^2));
        out=cumsum(out);
end