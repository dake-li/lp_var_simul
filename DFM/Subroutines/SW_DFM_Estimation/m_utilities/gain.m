function fgain = gain(h,w)

% Calculate Gain of filter h at Frequency w
im = sqrt(-1);
z=exp(-w*im);
h1=h(1);
z1=1;
for i = 2:size(h,1);
  z1=z1*z;
  h1=h1+h(i)*z1;
end;
g2=h1*h1';
fgain = sqrt(g2);  

end