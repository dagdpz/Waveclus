unit = 1;

for unit = 1 : size(indexspikemax,1)

for i = 1 : length(indexspikemax)
    
    tempshift(i,:) = conv(spikeall(i,:),-fliplr(spikeall(indexspikemax(unit),:)));
    
end

end
[~,tempshift] = min(spikeallconv,[],2);
tempshift = tempshift - 64;

figure;plot(spikeallconv')

int_factor = 2;
np = 127;


ints=1/int_factor:1/int_factor:np;
    intspikes = spline(1:np,tempshift,ints);
    
    figure;plot(intspikes')
    
    
    
    
    ints2 = [];
        for i = 1 : handles.par.scales - 1 
            divfac = 2^(i + 1);
            divfac2 = 2^i;
            offfac = size(handles.spikecomp,2)/(2^i);
            intsB = round((handles.par.w_pre/divfac2-handles.par.w_pre/divfac:handles.par.w_pre/divfac2+handles.par.w_post/divfac-1)+offfac);
            ints2 = cat(2,ints2,intsB);
        end
        intsB = round(handles.par.w_pre/divfac2-handles.par.w_pre/divfac:handles.par.w_pre/divfac2+handles.par.w_post/divfac-1);
        ints2 = cat(2,ints2,intsB);