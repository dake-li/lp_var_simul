function ismpl = smpl(calvec,nfirst,nlast,nper);
 i1 = nfirst(1) + (nfirst(2)-1)/nper;
 i2 = nlast(1) + (nlast(2)-1)/nper;
 ismpl = (calvec > i1-.0001).*(calvec < i2+.0001);
end

