
for ii = 1:size(VSDI.roi.labels,1)
    roiname  =VSDI.roi.labels{ii,1};
    VSDI.roi.labels{ii,2} = roiname_ipsicontra(roiname, VSDI.info.Sside); 
end

