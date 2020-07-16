function [istart,iend] = smpl_HO(smpl_par)

nfirst = smpl_par.nfirst;
nlast  = smpl_par.nlast;
calvec = smpl_par.calvec;
nper   = smpl_par.nper;

i1 = nfirst(1) + (nfirst(2)-1)/nper;
i2 = nlast(1) + (nlast(2)-1)/nper;
%  ismpl = (calvec > i1-.0001).*(calvec < i2+.0001);
istart = find(calvec>i1-.0001, 1 );
iend = find(calvec<i2+.0001, 1, 'last' );

end

